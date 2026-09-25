using MathNet.Numerics.Statistics;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using System.Globalization;
using System.Text.Json;
using Row = System.Collections.Generic.Dictionary<string, string>;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class Statistics
{
    static readonly JsonSerializerOptions JsonOptions = new() { WriteIndented = true };
    static readonly (string Left, string Right)[] Pairs = [("fixed_core", "nfw"), ("variable_core", "nfw"), ("variable_core", "fixed_core")];
    public static readonly Dictionary<string, string> FailureConvention = new()
    {
        ["invalid_prediction_loss"] = "positive_infinity",
        ["eligible_denominator"] = "all paired eligible galaxies; no prediction-based exclusions",
        ["one_invalid"] = "finite prediction wins",
        ["joint_invalid"] = "joint failure, neither model wins; not a finite tie",
        ["paired_difference"] = "omit only undefined joint-failure differences; report their count; opposite infinite central values give undefined median",
        ["ecdf"] = "finite loss steps divided by all eligible galaxies; failure mass remains at positive infinity",
        ["json_nonfinite"] = "positive_infinity, negative_infinity, undefined are explicit strings; absent/not-applicable values are null"
    };

    public static Result Run(string dataRoot, string outputRoot)
    {
        string numerics = Directory.Exists(Path.Combine(dataRoot, "numerics")) ? Path.Combine(dataRoot, "numerics") : dataRoot;
        string destination = Path.Combine(outputRoot, "significance");
        Directory.CreateDirectory(destination);
        var frames = new Dictionary<string, List<Row>>();
        var hashes = new Dictionary<string, string>();
        List<Row> Read(string relative)
        {
            if (!frames.TryGetValue(relative, out var rows))
            {
                string path = Path.Combine(numerics, relative);
                hashes[relative] = Data.FileSha256(path);
                frames[relative] = rows = Csv.Read(path);
            }
            return rows;
        }

        var primary = new List<Result>();
        var sensitivity = new List<Result>();
        var descriptive = new List<Result>();
        var paired = new List<Result>();
        var audits = new Result();
        var auditRows = new List<Result>();
        foreach (string epsilon in new[] { "0.001", "0.01" })
        {
            string basePath = $"regime_baseline/results/guard_{epsilon}_merged_fits.csv";
            var baseline = Read(basePath);
            var predictions = Read($"regime_baseline/results/guard_{epsilon}_predictions.csv");
            foreach (string model in new[] { "core", "nfw" })
            {
                var (report, rows) = AuditPredictions(baseline, predictions, model);
                audits[basePath + ":" + model] = report;
                auditRows.AddRange(rows.Select(r => WithPrefix(r, basePath, model)));
            }
            foreach (bool wide in new[] { false, true })
            {
                string familyPath = $"population_family/guard_{(wide ? "wide_" : "")}{epsilon}/free_family_fits.csv";
                var free = Read(familyPath);
                var freePredictions = Read(familyPath.Replace("_fits.csv", "_predictions.csv"));
                var (report, rows) = AuditPredictions(free, freePredictions);
                audits[familyPath] = report;
                auditRows.AddRange(rows.Select(r => WithPrefix(r, familyPath, "variable_core")));
                foreach (string model in new[] { "core", "nfw" })
                    audits[familyPath + ":paired:" + model] = CheckPairedObservations(freePredictions, predictions.Where(r => Text(r, "model") == model));
                var models = new Dictionary<string, List<Row>>
                {
                    ["fixed_core"] = baseline.Where(r => Text(r, "model") == "core").ToList(),
                    ["nfw"] = baseline.Where(r => Text(r, "model") == "nfw").ToList(),
                    ["variable_core"] = free
                };
                string scenario = $"guard_{epsilon}" + (wide ? "_wide_eta" : "");
                foreach (var (left, right) in Pairs)
                {
                    if (wide && left == "fixed_core" && right == "nfw")
                        continue;
                    var (result, losses) = Compare(Subset(models[left], true), Subset(models[right], true), left, right, scenario);
                    if (scenario == "guard_0.001")
                    {
                        result["tail_and_near_tie_diagnostics"] = TailDiagnostics(losses);
                        primary.Add(result);
                    }
                    else
                    {
                        result["multiplicity_status"] = "exploratory sensitivity; unadjusted";
                        sensitivity.Add(result);
                    }
                    paired.AddRange(losses);
                }
                foreach (var (left, right) in Pairs)
                {
                    if (wide && left == "fixed_core" && right == "nfw")
                        continue;
                    var (result, losses) = Compare(Subset(models[left], false), Subset(models[right], false), left, right, scenario, false, "full");
                    descriptive.Add(result);
                    paired.AddRange(losses);
                }
            }
        }
        double[] adjusted = Holm(primary.Select(r => (double)r["p_two_sided"]!).ToArray());
        for (int i = 0; i < primary.Count; i++)
        {
            primary[i]["p_holm_three_comparisons"] = adjusted[i];
            primary[i]["z_holm_two_sided"] = NormalMagnitude(adjusted[i]);
        }
        foreach (string directory in new[] { "uncertainty", "uncertainty_wide_pulls", "uncertainty_distance_floor" })
        {
            var frame = Read($"observations/{directory}/penalized_fits.csv");
            foreach (string configuration in frame.Select(r => Text(r, "configuration")).Distinct().Order(StringComparer.Ordinal))
                foreach (bool holdout in new[] { true, false })
                {
                    var left = Subset(frame, holdout, "core", configuration);
                    var right = Subset(frame, holdout, "nfw", configuration);
                    if (left.Count == 0 && right.Count == 0)
                        continue;
                    var (result, rows) = Compare(left, right, "fixed_core", "nfw", directory + ":" + configuration, holdout, holdout ? "holdout" : "full");
                    if (holdout)
                    {
                        result["multiplicity_status"] = "exploratory sensitivity; unadjusted";
                        sensitivity.Add(result);
                    }
                    else
                        descriptive.Add(result);
                    paired.AddRange(rows);
                }
        }
        var results = new Result
        {
            ["analysis_status"] = "exploratory conditional comparison of previously inspected saved model fits",
            ["statistical_unit"] = "one galaxy, never an individual radius",
            ["primary_family"] = "three paired outer-holdout comparisons at regime guard 0.001 and original eta range",
            ["null"] = "independent fair win/loss signs conditional on non-ties and non-joint failures",
            ["median_and_win_ci_assumption"] = "independent identically distributed galaxies from the specified selected population; exchangeability alone is insufficient; sample randomness not established",
            ["z_definition"] = "Phi^{-1}(1-p/2), a two-sided normal-equivalent magnitude; direction is separate",
            ["multiplicity_scope"] = "Holm only across the three declared main comparisons, not historical choices or sensitivity analyses",
            ["not_established"] = new[] { "absolute goodness of fit", "probability that a theory is true", "significance of full self-consistent baryon-field theory", "marginalization over nuisance parameters", "inter-galaxy independence", "a prospective independent validation sample" },
            ["full_fit_comparisons"] = "descriptive only; parameter flexibility and optimization act on these same points",
            ["failure_convention"] = FailureConvention,
            ["primary"] = primary,
            ["sensitivities"] = sensitivity,
            ["descriptive_full_fit"] = descriptive,
            ["prediction_reconstruction"] = audits,
            ["input_sha256"] = hashes,
            ["software_versions"] = new Result { ["runtime"] = System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription, ["language"] = "C#", ["mathnet_numerics"] = typeof(SpecialFunctions).Assembly.GetName().Version?.ToString() }
        };
        WriteJson(Path.Combine(destination, "results.json"), results);
        Csv.Write(Path.Combine(destination, "paired_galaxy_losses.csv"), paired.Select(SafeRow));
        Csv.Write(Path.Combine(destination, "prediction_reconstruction.csv"), auditRows.Select(SafeRow));
        string[] keys = ["left_model", "right_model", "eligible_galaxies", "left_wins", "right_wins", "finite_ties", "joint_failures", "p_two_sided", "z_two_sided", "p_holm_three_comparisons", "z_holm_two_sided"];
        var table = primary.Select(r =>
        {
            var row = keys.ToDictionary(k => k, k => r[k]);
            var ci = (double[])r["left_win_fraction_ci95"]!;
            var median = (Result)r["median_difference_interval"]!;
            var mci = (double[])median["ci95"]!;
            row["win_ci95_lower"] = ci[0];
            row["win_ci95_upper"] = ci[1];
            row["median_delta"] = median["median"];
            row["median_ci95_lower"] = mci[0];
            row["median_ci95_upper"] = mci[1];
            return SafeRow(row);
        }).ToList();
        Csv.Write(Path.Combine(destination, "primary_summary.csv"), table);
        var additional = AdditionalAnalyses(numerics, outputRoot);
        var validation = Validate(numerics, destination, results);
        validation["additional_analyses"] = additional;
        validation["self_checks"] = SelfChecks();
        validation["prediction_audits"] = audits.Count;
        validation["input_files"] = hashes.Count;
        WriteJson(Path.Combine(destination, "validation.json"), validation);
        Console.WriteLine($"Statistics: {primary.Count} primary, {sensitivity.Count} sensitivity, {descriptive.Count} descriptive comparisons; {audits.Count} audits passed.");
        return results;
    }

    static Result WithPrefix(Result row, string input, string model)
    {
        var result = new Result { ["input"] = input, ["model"] = model };
        foreach (var item in row)
            result[item.Key] = item.Value;
        return result;
    }
    static string Text(Row row, string key) => row.GetValueOrDefault(key, "");
    static bool Flag(Row row, string key) => string.Equals(Text(row, key), "true", StringComparison.OrdinalIgnoreCase);
    static double Number(Row row, string key) => TryNumber(Text(row, key), out double n) ? n : double.NaN;
    static int Count(Row row, string key)
    {
        double n = Number(row, key);
        if (!double.IsFinite(n) || n < 0 || n > int.MaxValue || n != Math.Truncate(n))
            throw new InvalidDataException("Invalid integer count: " + key);
        return (int)n;
    }
    public static double PredictiveLoss(Row row)
    {
        double value = Number(row, "holdout_chi2_per_row");
        return Flag(row, "holdout_admissible") && double.IsFinite(value) && value >= 0 ? value : double.PositiveInfinity;
    }
    static List<Row> Subset(IEnumerable<Row> rows, bool holdout, string? model = null, string? configuration = null) =>
        rows.Where(r => Flag(r, "holdout") == holdout && (model is null || Text(r, "model") == model) && (configuration is null || Text(r, "configuration") == configuration)).ToList();

    static (Result Summary, List<Result> Rows) Compare(IEnumerable<Row> leftInput, IEnumerable<Row> rightInput, string leftName, string rightName, string name, bool inferential = true, string split = "holdout")
    {
        List<Row> Normalize(IEnumerable<Row> source) => source.Select(r =>
        {
            var row = new Row(r);
            foreach (string key in new[] { "n", "n_train", "n_test" })
                if (row.ContainsKey(key))
                    Count(row, key);
            if (split == "full")
            {
                row["n_test"] = row["n_train"];
                row["holdout_chi2_per_row"] = row["chi2_data_per_row"];
                row["holdout_admissible"] = row.GetValueOrDefault("train_admissible", "True");
            }
            return row;
        }).ToList();
        var left = Normalize(leftInput);
        var right = Normalize(rightInput);
        if (left.Select(r => Text(r, "galaxy")).Distinct().Count() != left.Count || right.Select(r => Text(r, "galaxy")).Distinct().Count() != right.Count)
            throw new InvalidDataException("Paired summaries require one row per galaxy.");
        var a = left.ToDictionary(r => Text(r, "galaxy"));
        var b = right.ToDictionary(r => Text(r, "galaxy"));
        if (!a.Keys.ToHashSet().SetEquals(b.Keys))
            throw new InvalidDataException("Paired summaries require identical eligible galaxy sets.");
        int points = 0, finite = 0, leftWins = 0, rightWins = 0, ties = 0, jointFailures = 0, leftFailures = 0, rightFailures = 0;
        var differences = new List<double>();
        var rows = new List<Result>();
        foreach (string galaxy in a.Keys.Order(StringComparer.Ordinal))
        {
            int count = Count(a[galaxy], "n_test");
            if (count <= 0 || count != Count(b[galaxy], "n_test"))
                throw new InvalidDataException("Paired predictions require equal positive point counts.");
            points += count;
            double x = PredictiveLoss(a[galaxy]), y = PredictiveLoss(b[galaxy]);
            double delta = x - y;
            bool fx = double.IsFinite(x), fy = double.IsFinite(y);
            string outcome;
            if (!fx && !fy)
            {
                jointFailures++;
                outcome = "joint_failure";
            }
            else
            {
                if (fx && fy)
                {
                    finite++;
                    if (x == y)
                        ties++;
                }
                else if (!fx)
                    leftFailures++;
                else
                    rightFailures++;
                if (x < y)
                    leftWins++;
                else if (y < x)
                    rightWins++;
                differences.Add(delta);
                outcome = x < y ? "left_win" : y < x ? "right_win" : "finite_tie";
            }
            rows.Add(new Result { ["analysis"] = name, ["split"] = split, ["left_model"] = leftName, ["right_model"] = rightName, ["galaxy"] = galaxy, ["n_points"] = count, ["left_loss"] = x, ["right_loss"] = y, ["delta_left_minus_right"] = delta, ["outcome"] = outcome });
        }
        var result = new Result
        {
            ["analysis"] = name,
            ["split"] = split,
            ["left_model"] = leftName,
            ["right_model"] = rightName,
            ["eligible_galaxies"] = a.Count,
            ["points"] = points,
            ["jointly_admissible"] = finite,
            ["left_wins"] = leftWins,
            ["right_wins"] = rightWins,
            ["finite_ties"] = ties,
            ["joint_failures"] = jointFailures,
            ["left_only_failures"] = leftFailures,
            ["right_only_failures"] = rightFailures,
            ["paired_difference_galaxies"] = differences.Count,
            ["median_per_galaxy_delta_chi2_per_point"] = Median(differences),
            ["left_median_loss"] = Median(left.Select(PredictiveLoss)),
            ["right_median_loss"] = Median(right.Select(PredictiveLoss))
        };
        if (inferential)
        {
            foreach (var item in ExactSign(leftWins, rightWins))
                result[item.Key] = item.Value;
            result["median_difference_interval"] = MedianInterval(differences);
            result["direction"] = leftWins > rightWins ? leftName : rightWins > leftWins ? rightName : "no win imbalance";
        }
        else
            result["inference"] = "descriptive in-sample comparison; no significance claim";
        return (result, rows);
    }

    static (Result Report, List<Result> Rows) AuditPredictions(IEnumerable<Row> fitInput, IEnumerable<Row> predictionInput, string? model = null)
    {
        var fits = fitInput.Where(r => model is null || Text(r, "model") == model).ToList();
        var predictions = predictionInput.Where(r => model is null || Text(r, "model") == model).ToList();
        string trainingKey = predictions.Count > 0 && predictions[0].ContainsKey("training_row") ? "training_row" : "is_training";
        var grouped = predictions.GroupBy(r => (Galaxy: Text(r, "galaxy"), Holdout: Flag(r, "holdout"))).ToDictionary(g => g.Key, g => g.ToList());
        var seen = new HashSet<(string Galaxy, bool Holdout)>();
        var reports = new List<Result>();
        double maxAbs = 0, maxRel = 0;
        foreach (var fit in fits)
        {
            var key = (Galaxy: Text(fit, "galaxy"), Holdout: Flag(fit, "holdout"));
            if (!seen.Add(key))
                throw new InvalidDataException("Duplicate fit key in prediction audit.");
            if (!grouped.TryGetValue(key, out var points) || points.Count != Count(fit, "n"))
                throw new InvalidDataException("Incorrect complete prediction count: " + key);
            if (points.Select(r => Number(r, "r_catalogue_kpc")).Distinct().Count() != points.Count)
                throw new InvalidDataException("Duplicate radius within galaxy/split.");
            foreach (string partition in key.Holdout ? new[] { "train", "test" } : new[] { "train" })
            {
                var picked = points.Where(r => Flag(r, trainingKey) == (partition == "train")).ToList();
                int count = Count(fit, partition == "train" ? "n_train" : "n_test");
                if (picked.Count != count || count <= 0)
                    throw new InvalidDataException("Incorrect train/test prediction count: " + key);
                double sum = 0;
                bool finite = true;
                foreach (var point in picked)
                {
                    double err = Number(point, "e_Vobs_kms"), obs = Number(point, "Vobs_kms"), pred = Number(point, "Vpred_catalogue_kms");
                    if (!double.IsFinite(err) || err <= 0 || !double.IsFinite(obs))
                        throw new InvalidDataException("Invalid observed data in prediction audit.");
                    finite &= double.IsFinite(pred) && pred >= 0;
                    double residual = (pred - obs) / err;
                    sum += residual * residual;
                }
                double reconstructed = finite ? sum / count : double.PositiveInfinity;
                double saved = partition == "train" ? Number(fit, "chi2_data_per_row") : PredictiveLoss(fit);
                double discrepancy = 0, relative = 0;
                if (!double.IsFinite(reconstructed) || !double.IsFinite(saved))
                {
                    if (reconstructed != saved)
                        throw new InvalidDataException("Admissibility mismatch in prediction audit.");
                }
                else
                {
                    discrepancy = Math.Abs(reconstructed - saved);
                    relative = discrepancy / Math.Max(1, Math.Abs(saved));
                    if (discrepancy > 1e-10 + 1e-10 * Math.Abs(saved))
                        throw new InvalidDataException($"Stored loss mismatch: {key}, {partition}, {reconstructed:R}, {saved:R}");
                }
                maxAbs = Math.Max(maxAbs, discrepancy);
                maxRel = Math.Max(maxRel, relative);
                reports.Add(new Result { ["galaxy"] = key.Galaxy, ["holdout"] = key.Holdout, ["partition"] = partition, ["n_points"] = count, ["saved_loss"] = saved, ["reconstructed_loss"] = reconstructed, ["absolute_difference"] = discrepancy });
            }
        }
        if (!seen.SetEquals(grouped.Keys))
            throw new InvalidDataException("Extra or missing prediction galaxy/split.");
        return (new Result { ["fit_rows"] = fits.Count, ["partition_checks"] = reports.Count, ["prediction_rows"] = predictions.Count, ["max_absolute_difference"] = maxAbs, ["max_relative_difference_scale_at_least_one"] = maxRel, ["passed"] = true }, reports);
    }

    static Result CheckPairedObservations(IEnumerable<Row> left, IEnumerable<Row> right)
    {
        static (string Galaxy, double Radius, double Observed, double Error, bool Training)[] Canonical(IEnumerable<Row> source) => source.Where(r => Flag(r, "holdout"))
            .Select(r => (Text(r, "galaxy"), Number(r, "r_catalogue_kpc"), Number(r, "Vobs_kms"), Number(r, "e_Vobs_kms"), Flag(r, r.ContainsKey("training_row") ? "training_row" : "is_training")))
            .OrderBy(r => r.Item1, StringComparer.Ordinal).ThenBy(r => r.Item2).ToArray();
        var a = Canonical(left);
        var b = Canonical(right);
        if (!a.SequenceEqual(b))
            throw new InvalidDataException("Compared models have unequal observations or partitions.");
        return new Result { ["matched_rows"] = a.Length, ["heldout_rows"] = a.Count(r => !r.Training), ["galaxies"] = a.Select(r => r.Galaxy).Distinct().Count(), ["passed"] = true };
    }

    static Result TailDiagnostics(IEnumerable<Result> source)
    {
        var rows = source.Where(r => double.IsFinite((double)r["left_loss"]!) && double.IsFinite((double)r["right_loss"]!)).ToList();
        if (rows.Count == 0)
            return new Result { ["jointly_finite"] = 0 };
        double[] x = rows.Select(r => (double)r["left_loss"]!).ToArray(), y = rows.Select(r => (double)r["right_loss"]!).ToArray();
        double[] weights = rows.Select(r => Convert.ToDouble(r["n_points"], CultureInfo.InvariantCulture)).ToArray();
        double[] delta = x.Zip(y, (a, b) => a - b).ToArray();
        var cuts = new List<Result>();
        foreach (string mode in new[] { "absolute", "relative_scale_at_least_one" })
            foreach (double level in mode == "absolute" ? new[] { 1e-6, 1e-4, .01 } : new[] { 1e-6, 1e-4 })
            {
                int near = 0, left = 0, right = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    double threshold = mode == "absolute" ? level : level * Math.Max(1, Math.Max(Math.Abs(x[i]), Math.Abs(y[i])));
                    if (Math.Abs(delta[i]) <= threshold)
                        near++;
                    if (delta[i] < -threshold)
                        left++;
                    if (delta[i] > threshold)
                        right++;
                }
                cuts.Add(new Result { ["mode"] = mode, ["tolerance"] = level, ["near_ties"] = near, ["left_lower_beyond_tolerance"] = left, ["right_lower_beyond_tolerance"] = right });
            }
        double[] q = [.5, .75, .9, .95, .99, 1], dq = [0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1];
        return new Result
        {
            ["jointly_finite"] = rows.Count,
            ["mean_delta_galaxy_equal"] = delta.Average(),
            ["point_weighted_mean_delta"] = delta.Zip(weights, (d, w) => d * w).Sum() / weights.Sum(),
            ["left_mean_loss"] = x.Average(),
            ["right_mean_loss"] = y.Average(),
            ["quantile_method"] = "numpy linear interpolation",
            ["quantile_levels"] = q,
            ["left_loss_quantiles"] = Quantiles(x, q),
            ["right_loss_quantiles"] = Quantiles(y, q),
            ["delta_quantiles"] = Quantiles(delta, dq),
            ["delta_quantile_levels"] = dq,
            ["left_improvements_over_one"] = delta.Count(d => d < -1),
            ["left_regressions_over_one"] = delta.Count(d => d > 1),
            ["left_improvements_over_ten"] = delta.Count(d => d < -10),
            ["left_regressions_over_ten"] = delta.Count(d => d > 10),
            ["near_tie_diagnostics"] = cuts,
            ["near_tie_scope"] = "descriptive sign stability; thresholded comparisons are not extra significance tests",
            ["largest_regressions"] = rows.OrderByDescending(r => (double)r["delta_left_minus_right"]!).Take(5).ToList(),
            ["largest_improvements"] = rows.OrderBy(r => (double)r["delta_left_minus_right"]!).Take(5).ToList()
        };
    }

    public static double Median(IEnumerable<double> source)
    {
        double[] a = source.Order().ToArray();
        if (a.Length == 0)
            return double.NaN;
        int m = a.Length / 2;
        if (a.Length % 2 == 1)
            return a[m];
        if (a[m - 1] == a[m])
            return a[m];
        return (a[m - 1] + a[m]) / 2;
    }
    static double[] Quantiles(double[] values, double[] levels)
    {
        var sorted = values.Order().ToArray();
        return levels.Select(level => SortedArrayStatistics.QuantileCustom(sorted, level, QuantileDefinition.R7)).ToArray();
    }
    public static Result ExactSign(int wins, int losses)
    {
        if (wins < 0 || losses < 0)
            throw new ArgumentOutOfRangeException(nameof(wins));
        int informative = checked(wins + losses);
        if (informative == 0)
            return new Result
            {
                ["n_informative"] = 0,
                ["left_win_fraction"] = null,
                ["left_win_fraction_ci95"] = null,
                ["p_two_sided"] = null,
                ["z_two_sided"] = null
            };
        double probability = Math.Min(1, 2 * BinomialCdf(Math.Min(wins, losses), informative));
        // Clopper-Pearson bounds include the exact zero-win and zero-loss endpoints.
        double lower = wins == 0 ? 0 : BetaInverse(.025, wins, losses + 1);
        double upper = losses == 0 ? 1 : BetaInverse(.975, wins + 1, losses);
        return new Result
        {
            ["n_informative"] = informative,
            ["left_win_fraction"] = (double)wins / informative,
            ["left_win_fraction_ci95"] = new[] { lower, upper },
            ["p_two_sided"] = probability,
            ["z_two_sided"] = NormalMagnitude(probability)
        };
    }

    public static double[] Holm(double[] p)
    {
        if (p.Any(v => !double.IsFinite(v) || v < 0 || v > 1))
            throw new ArgumentException("Holm adjustment requires finite probabilities.");
        var adjusted = new double[p.Length];
        double running = 0;
        int[] order = Enumerable.Range(0, p.Length).OrderBy(i => p[i]).ToArray();
        // Holm's step-down adjustment is monotone in sorted p-value order.
        for (int rank = 0; rank < p.Length; rank++)
        {
            int originalIndex = order[rank];
            running = Math.Max(running, (p.Length - rank) * p[originalIndex]);
            adjusted[originalIndex] = Math.Min(1, running);
        }
        return adjusted;
    }
    static Result MedianInterval(IEnumerable<double> source)
    {
        double[] values = source.Where(v => !double.IsNaN(v)).Order().ToArray();
        int n = values.Length;
        if (n == 0)
            return new Result { ["n"] = 0, ["median"] = null, ["ci95"] = null, ["k_1based"] = null, ["coverage_continuous_null"] = null };
        int k = 0;
        // Exact sign-test inversion gives the median interval's order-statistic ranks.
        for (int candidate = 1; candidate <= (n + 1) / 2; candidate++)
        {
            if (2 * BinomialCdf(candidate - 1, n) <= .05)
                k = candidate;
        }
        return new Result { ["n"] = n, ["median"] = Median(values), ["ci95"] = k > 0 ? new[] { values[k - 1], values[n - k] } : new[] { double.NegativeInfinity, double.PositiveInfinity }, ["k_1based"] = k, ["coverage_continuous_null"] = k > 0 ? 1 - 2 * BinomialCdf(k - 1, n) : 1.0 };
    }
    static double BinomialCdf(int k, int n)
    {
        if (k < 0)
            return 0;
        if (k >= n)
            return 1;
        // Preserve the exactly representable all-wins/all-losses probability.
        return k == 0 ? Math.Pow(.5, n) : Binomial.CDF(.5, n, k);
    }

    public static double BetaInverse(double probability, double a, double b) => Beta.InvCDF(a, b, probability);

    public static double RegularizedBeta(double x, double a, double b)
    {
        if (x <= 0)
            return 0;
        if (x >= 1)
            return 1;
        return SpecialFunctions.BetaRegularized(a, b, x);
    }

    public static double NormalMagnitude(double p)
    {
        if (p == 0)
            return double.PositiveInfinity;
        if (p == 1)
            return 0;
        // Erfc inversion avoids subtracting a tiny tail probability from one.
        return Math.Sqrt(2) * SpecialFunctions.ErfcInv(p);
    }
    static object? JsonSafe(object? value)
    {
        if (value is double d && !double.IsFinite(d))
            return double.IsNaN(d) ? "undefined" : d > 0 ? "positive_infinity" : "negative_infinity";
        if (value is Result dictionary)
            return dictionary.ToDictionary(k => k.Key, k => JsonSafe(k.Value));
        if (value is System.Collections.IEnumerable list && value is not string && value is not System.Collections.IDictionary)
            return list.Cast<object?>().Select(JsonSafe).ToArray();
        return value;
    }
    static Result SafeRow(Result row) => row.ToDictionary(k => k.Key, k => JsonSafe(k.Value));
    static void WriteJson(string path, object value) => File.WriteAllText(path, JsonSerializer.Serialize(JsonSafe(value), JsonOptions) + Environment.NewLine);
}
