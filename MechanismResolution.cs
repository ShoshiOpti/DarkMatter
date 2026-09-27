using System.Globalization;
using Result = System.Collections.Generic.Dictionary<string, object?>;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

/// <summary>
/// Resolution audit for the saved charged/zero-contact control. Equivalence is
/// never inferred from a nonsignificant p value or from an unspecified tolerance.
/// </summary>
public static class MechanismResolution
{
    sealed record Point(string Galaxy, int Index, double Radius, double Observed, double Sigma, bool Outer,
        string Model, double Predicted, bool Failed);

    public static Result Run(string dataRoot, string outputRoot, AnalysisStudyProtocol? protocol = null)
    {
        protocol ??= new(); protocol.Validate();
        if (protocol.TestAlpha >= .5) throw new ArgumentException("Mechanism interval inclusion requires TestAlpha below one half.");
        string input = Path.Combine(dataRoot, "publication", "statistics", "charged", "predicted_rows.csv");
        string output = Path.Combine(outputRoot, "research", "mechanism_resolution");
        Directory.CreateDirectory(output);
        var points = Csv.Read(input).Where(r => r["model"] is "charged" or "free_control").Select(Parse).ToArray();
        var byKey = new Dictionary<(string, int, string), Point>();
        foreach (var point in points)
            if (!byKey.TryAdd((point.Galaxy, point.Index, point.Model), point))
                throw new InvalidDataException("Duplicate charged/control prediction row.");
        var a = points.Where(p => p.Model == "charged").OrderBy(p => p.Galaxy, StringComparer.Ordinal).ThenBy(p => p.Index).ToArray();
        var b = new List<Point>();
        foreach (var left in a)
        {
            if (!byKey.TryGetValue((left.Galaxy, left.Index, "free_control"), out var right))
                throw new InvalidDataException("Charged and control observations are not paired.");
            if (left.Radius != right.Radius || left.Observed != right.Observed || left.Sigma != right.Sigma || left.Outer != right.Outer)
                throw new InvalidDataException("Charged and control observation coordinates or split differ.");
            b.Add(right);
        }
        if (a.Length != points.Count(p => p.Model == "free_control") || a.Length != 3034 ||
            a.Count(p => p.Outer) != 659 || a.Select(p => p.Galaxy).Distinct().Count() != 131)
            throw new InvalidDataException("The retained charged/control resolution study requires the frozen 131-galaxy, 3034-observation, 659-outer sample.");
        var failed = points.GroupBy(p => p.Galaxy).Where(g => g.Any(p => p.Failed || p.Outer && (!double.IsFinite(p.Predicted) || p.Predicted < 0)))
            .Select(g => g.Key).ToHashSet(StringComparer.Ordinal);
        var pointRows = new List<Result>();
        var galaxyInputs = new List<PairedGalaxyObservation>();
        var galaxyRows = new List<Result>();
        foreach (var group in a.Where(p => p.Outer).GroupBy(p => p.Galaxy))
        {
            bool valid = !failed.Contains(group.Key);
            var stable = new List<double>(); var old = new List<double>(); var speed = new List<double>(); var errorBounds = new List<double>();
            foreach (var left in group)
            {
                var right = byKey[(left.Galaxy, left.Index, "free_control")];
                double delta = valid ? StableSquaredLossDifference(left.Predicted, right.Predicted, left.Observed, left.Sigma) : double.NaN;
                double legacy = valid ? Math.Pow((left.Predicted - left.Observed) / left.Sigma, 2) - Math.Pow((right.Predicted - left.Observed) / left.Sigma, 2) : double.NaN;
                double dv = valid ? left.Predicted - right.Predicted : double.NaN;
                double? budget = valid && protocol.SolverSpeedErrorBoundKms is double epsilon
                    ? SquaredLossErrorBound(left.Predicted, right.Predicted, left.Observed, left.Sigma, epsilon) : null;
                if (valid) { stable.Add(delta); old.Add(legacy); speed.Add(Math.Abs(dv)); if (budget is double bound) errorBounds.Add(bound); }
                pointRows.Add(new()
                {
                    ["galaxy"] = left.Galaxy, ["index"] = left.Index, ["radius_kpc"] = left.Radius,
                    ["observed_kms"] = left.Observed, ["sigma_kms"] = left.Sigma,
                    ["charged_kms"] = double.IsFinite(left.Predicted) ? left.Predicted : null,
                    ["free_control_kms"] = double.IsFinite(right.Predicted) ? right.Predicted : null,
                    ["speed_difference_kms"] = valid ? dv : null, ["absolute_speed_difference_kms"] = valid ? Math.Abs(dv) : null,
                    ["stable_squared_loss_difference"] = valid ? delta : null,
                    ["subtracted_squared_loss_difference"] = valid ? legacy : null,
                    ["supplied_solver_loss_error_bound"] = budget, ["comparison_valid"] = valid
                });
            }
            double? effect = valid ? StatisticalEngine.Sum(stable) / stable.Count : null;
            double? boundMean = errorBounds.Count > 0 ? StatisticalEngine.Sum(errorBounds) / errorBounds.Count : null;
            galaxyRows.Add(new()
            {
                ["galaxy"] = group.Key, ["outer_points"] = group.Count(), ["comparison_valid"] = valid,
                ["stable_squared_loss_difference"] = effect,
                ["subtracted_squared_loss_difference"] = valid ? StatisticalEngine.Sum(old) / old.Count : null,
                ["maximum_absolute_speed_difference_kms"] = valid ? speed.Max() : null,
                ["mean_absolute_speed_difference_kms"] = valid ? StatisticalEngine.Sum(speed) / speed.Count : null,
                ["supplied_solver_loss_error_bound"] = boundMean
            });
            if (valid) galaxyInputs.Add(new(group.Key, group.Count(), effect!.Value));
        }
        bool complete = failed.Count == 0;
        var options = new PairedBootstrapOptions
        {
            Draws = protocol.BootstrapDraws, ConfidenceLevel = protocol.ConfidenceLevel,
            PairsSeed = StatisticalEngine.DeriveSeed(protocol.Seed, 9000000001),
            WildSeed = StatisticalEngine.DeriveSeed(protocol.Seed, 9000000002)
        };
        PairedBootstrapResult? inference = complete ? StatisticalEngine.Run(galaxyInputs, options) : null;
        PairedBootstrapResult? inclusionInference = complete && protocol.MechanismLossMargin is not null
            ? StatisticalEngine.Run(galaxyInputs, options with { ConfidenceLevel = 1 - 2 * protocol.TestAlpha }) : null;
        double[] differences = pointRows.Where(r => r["comparison_valid"] is true).Select(r => (double)r["absolute_speed_difference_kms"]!).Order().ToArray();
        double? maximumSpeed = differences.Length > 0 ? differences[^1] : null;
        double? lossErrorBudget = complete && protocol.SolverSpeedErrorBoundKms is not null
            ? galaxyRows.Average(r => (double)r["supplied_solver_loss_error_bound"]!) : null;
        double? maxSpeedIncludingError = complete && maximumSpeed is double maximum && protocol.SolverSpeedErrorBoundKms is double error ? maximum + 2 * error : null;
        var interval = inclusionInference?.BootstrapTInterval;
        StatisticalInterval? expandedInterval = interval is not null && lossErrorBudget is double budgetLoss
            ? new(interval.Lower - budgetLoss, interval.Upper + budgetLoss) : null;
        bool? speedWithin = protocol.MechanismSpeedToleranceKms is double speedMargin && maxSpeedIncludingError is double upperSpeed ? upperSpeed <= speedMargin : null;
        bool? lossWithin = protocol.MechanismLossMargin is double lossMargin && expandedInterval is not null
            ? expandedInterval.Lower > -lossMargin && expandedInterval.Upper < lossMargin : null;
        string Status(double? margin, bool? outcome) => !complete ? "undefined_full_sample_due_to_failed_pair"
            : margin is null ? "no_margin_specified"
            : protocol.SolverSpeedErrorBoundKms is null ? "solver_error_bound_not_specified"
            : outcome is null ? "undefined_interval_or_budget"
            : outcome.Value ? "within_explicitly_configured_margin" : "not_within_explicitly_configured_margin";
        Csv.Write(Path.Combine(output, "pointwise_differences.csv"), pointRows);
        Csv.Write(Path.Combine(output, "galaxy_contrasts.csv"), galaxyRows);
        var report = new Result
        {
            ["schema_version"] = 1, ["status"] = "completed", ["protocol"] = protocol,
            ["source"] = Path.GetRelativePath(dataRoot, input).Replace('\\', '/'), ["source_sha256"] = Data.FileSha256(input),
            ["total_galaxies"] = 131, ["outer_points"] = 659, ["failed_galaxies"] = failed.Order(StringComparer.Ordinal).ToArray(),
            ["full_sample_continuous_defined"] = complete,
            ["scope"] = "Resolution of frozen charged/control predictions. This neither estimates solver error nor demonstrates physical equivalence or population formation.",
            ["stable_contrast_formula"] = "(v_charged-v_free)*(v_charged+v_free-2*v_observed)/sigma^2, averaged within galaxy and then equally across galaxies",
            ["maximum_absolute_speed_difference_kms"] = maximumSpeed,
            ["median_absolute_speed_difference_kms"] = differences.Length > 0 ? StatisticalEngine.Quantile(differences, .5) : null,
            ["p95_absolute_speed_difference_kms"] = differences.Length > 0 ? StatisticalEngine.Quantile(differences, .95) : null,
            ["stable_effect_inference"] = inference,
            ["solver_budget"] = new Result
            {
                ["per_prediction_speed_error_bound_kms"] = protocol.SolverSpeedErrorBoundKms,
                ["provenance"] = "User-supplied bound, not an internally estimated or independently certified error. Absent bounds remain absent.",
                ["effect_error_bound"] = lossErrorBudget, ["maximum_speed_difference_including_error_kms"] = maxSpeedIncludingError,
                ["formula"] = "Per outer row: [2*|vA-vobs|*epsilon + epsilon^2 + 2*|vB-vobs|*epsilon + epsilon^2]/sigma^2; equal-galaxy averaging bounds the mean effect. Bounds are added adversarially, without assuming error independence."
            },
            ["speed_margin_assessment"] = new Result
            {
                ["status"] = Status(protocol.MechanismSpeedToleranceKms, speedWithin),
                ["configured_margin_kms"] = protocol.MechanismSpeedToleranceKms, ["within_margin_including_supplied_error"] = speedWithin
            },
            ["loss_margin_assessment"] = new Result
            {
                ["status"] = Status(protocol.MechanismLossMargin, lossWithin), ["configured_symmetric_margin"] = protocol.MechanismLossMargin,
                ["confidence_level"] = 1 - 2 * protocol.TestAlpha, ["bootstrap_t_interval"] = interval,
                ["interval_expanded_by_supplied_solver_error"] = expandedInterval, ["within_margin_including_supplied_error"] = lossWithin,
                ["interpretation"] = "An optional two-one-sided-test-style interval-inclusion diagnostic using a marginal pairs-bootstrap-t interval; not an inversion of the wild test. Conditional on the selected margin, claimed solver bound, independent galaxies and bootstrap calibration."
            },
            ["no_equivalence_from_nonsignificance"] = true,
            ["published_reproduction_unchanged"] = true
        };
        Data.SaveJson(Path.Combine(output, "report.json"), report);
        Console.WriteLine($"Mechanism resolution: maximum saved speed difference {maximumSpeed:G6} km/s; speed assessment {Status(protocol.MechanismSpeedToleranceKms, speedWithin)}.");
        return report;
    }

    static Point Parse(Row row)
    {
        bool Flag(string value) => value.ToLowerInvariant() switch
        {
            "true" or "1" or "1.0" => true, "false" or "0" or "0.0" => false,
            _ => throw new InvalidDataException("Invalid prediction flag.")
        };
        int index = int.Parse(row["index"], CultureInfo.InvariantCulture);
        double radius = Csv.Number(row, "radius"), observed = Csv.Number(row, "observed"), sigma = Csv.Number(row, "sigma");
        if (index < 0 || !double.IsFinite(radius) || radius <= 0 || !double.IsFinite(observed) || !double.IsFinite(sigma) || sigma <= 0)
            throw new InvalidDataException("Invalid charged/control observational coordinate.");
        return new(row["galaxy"], index, radius, observed, sigma, Flag(row["is_outer"]), row["model"],
            Csv.Number(row, "predicted"), row.ContainsKey("fit_failed") && Flag(row["fit_failed"]));
    }

    static double StableSquaredLossDifference(double a, double b, double observed, double sigma) =>
        ((a - b) / sigma) * (((a - observed) + (b - observed)) / sigma);
    static double SquaredLossErrorBound(double a, double b, double observed, double sigma, double epsilon) =>
        (2 * Math.Abs((a - observed) / sigma) + 2 * Math.Abs((b - observed) / sigma)) * (epsilon / sigma) + 2 * Math.Pow(epsilon / sigma, 2);

    public static Result SelfCheck()
    {
        int count = 0;
        void Check(bool ok, string label) { if (!ok) throw new InvalidDataException("Mechanism resolution self-check: " + label); count++; }
        Check(StableSquaredLossDifference(4, 3, 1, 2) == 1.25, "analytic squared-loss difference");
        Check(StableSquaredLossDifference(3, 4, 1, 2) == -1.25, "swapped control reverses sign");
        Check(StableSquaredLossDifference(4, 4, 1, 2) == 0, "identical predictions give exactly zero");
        Check(StableSquaredLossDifference(40, 30, 10, 20) == 1.25, "unit scaling invariant");
        double baseline = StableSquaredLossDifference(20, 21, 19, 2), bound = SquaredLossErrorBound(20, 21, 19, 2, .1);
        foreach (double ea in new[] { -.1, 0, .1 })
        foreach (double eb in new[] { -.1, 0, .1 })
            Check(Math.Abs(StableSquaredLossDifference(20 + ea, 21 + eb, 19, 2) - baseline) <= bound + 1e-13, "perturbed prediction lies inside supplied error bound");
        Check(SquaredLossErrorBound(20, 21, 19, 2, 0) == 0, "zero supplied error produces zero budget");
        return new() { ["status"] = "pass", ["checks"] = count };
    }
}
