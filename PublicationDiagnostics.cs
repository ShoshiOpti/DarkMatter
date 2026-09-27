using System.Text.Json;
using Row = System.Collections.Generic.Dictionary<string, string>;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

/// <summary>Recomputes frozen-target diagnostics; it does not solve or refit a halo.</summary>
public static class PublicationDiagnostics
{
    const double G = 4.300917270e-6, C = 299792.458;
    static double N(Row row, string key) => Csv.Number(row, key);
    static bool B(Row row, string key) => Csv.Flag(row, key);
    static double Value(Result row, string key) => Convert.ToDouble(row[key]);
    static bool Flag(Result row, string key) => (bool)row[key]!;

    public static Result Run(string dataRoot, string outputRoot)
    {
        string source = Path.Combine(dataRoot, "publication", "diagnostics");
        string output = Path.Combine(outputRoot, "publication_diagnostics");
        var inputs = Csv.Read(Path.Combine(source, "point_requirements.csv"));
        var galaxies = Csv.Read(Path.Combine(source, "galaxy_requirements.csv")).ToDictionary(r => r["galaxy"]);
        var expectedPoints = Csv.Read(Path.Combine(source, "expected", "point_closure.csv"))
            .ToDictionary(r => (r["galaxy"], r["index"]));
        var expectedGalaxies = Csv.Read(Path.Combine(source, "expected", "galaxy_closure.csv"))
            .ToDictionary(r => r["galaxy"]);
        var expectedSummary = Data.Json(Path.Combine(source, "expected", "closure_summary.json"));
        double mass = Data.Json(Path.Combine(source, "inverse_summary.json")).GetProperty("m_eV").GetDouble();
        double kappa = 1.054571817e-34 * Math.Pow(299792458, 2) /
            (1.602176634e-19 * 3.085677581491367e19 * 1000) / mass;
        double alpha = -3.0 / 16 * 4 * Math.PI * G * Math.Pow(kappa / C, 2);
        double r = .5 * (Math.Tanh(.5) / .5 + 1 / Math.Pow(Math.Cosh(.5), 2));
        double f = (1 + 6 * r) / (1 + 11 * r), dc = (1 - r) * (1 + 11 * r);
        double ap = r * (3 - 2 * r), bg = 12 * r * r * (3 * r - 2) / (1 + 11 * r);
        int checks = 0;
        double maxScaledError = 0;
        double maximumPhaseRoundingAllowance = 0;
        Dictionary<string, double> phaseRoundingAllowances = [];
        void Close(double actual, double expected, string label)
        {
            checks++;
            if (actual == expected || double.IsNaN(actual) && double.IsNaN(expected)) return;
            double tolerance = (label.Contains("mass_eV", StringComparison.OrdinalIgnoreCase) ? 1e-35 : 1e-10) + 3e-10 * Math.Abs(expected);
            string field = label[(label.LastIndexOf(':') + 1)..];
            // The archive's pandas CSV reader rounded a few h values within one
            // ulp of 1 to exactly 1. Propagate a four-unit-scale-ulp interval
            // through only the phase formulas, whose sqrt(1-h) is not Lipschitz
            // at that endpoint. Other algebra and all discrete decisions retain
            // the tighter ordinary tolerance; this is not a scientific error bar.
            tolerance += phaseRoundingAllowances.GetValueOrDefault(field);
            double error = Math.Abs(actual - expected);
            maxScaledError = Math.Max(maxScaledError, error / tolerance);
            if (!double.IsFinite(error) || error > tolerance)
                throw new InvalidDataException($"Publication diagnostic mismatch: {label}: {actual:R} != {expected:R}");
        }
        void Require(bool value, string label)
        {
            checks++;
            if (!value) throw new InvalidDataException("Publication diagnostic failed: " + label);
        }
        var points = new List<Result>();
        foreach (var row in inputs)
        {
            double radius = N(row, "R_kpc"), v2 = Math.Pow(N(row, "target_velocity_physical_kms"), 2);
            double rho = N(row, "target_density_msun_kpc3"), d1 = N(row, "target_density_prime");
            double d2 = N(row, "target_density_second"), d3 = N(row, "target_density_third");
            double h1 = d1 / rho, h2 = d2 / rho - h1 * h1;
            double h3 = d3 / rho - 3 * d2 * d1 / (rho * rho) + 2 * h1 * h1 * h1;
            double lap = .5 * h2 + .25 * h1 * h1 + h1 / radius;
            double lap1 = .5 * h3 + .5 * h1 * h2 + h2 / radius - h1 / (radius * radius);
            double potential = alpha * rho - .5 * kappa * kappa * lap;
            double derivative = alpha * d1 - .5 * kappa * kappa * lap1;
            Close(potential, N(row, "scalar_potential_kms2"), "scalar potential");
            Close(derivative, N(row, "scalar_potential_derivative"), "scalar potential derivative");
            // Retained derivative-error/domain flags are inputs, not new finite-difference certificates.
            double x = -radius * derivative / v2, q = v2 + radius * derivative;
            Close(q, N(row, "required_tangential_q_kms2"), "signed support");
            double h = N(row, "projected_orthogonal_density_fraction");
            Require(h >= 0 && h <= 1, "orthogonal fraction");
            double parserInterval = 4 * Math.ScaleB(1, -52);
            var centralPhase = Phase(h, f, dc, ap, bg);
            var lowerPhase = Phase(Math.Max(0, h - parserInterval), f, dc, ap, bg);
            var upperPhase = Phase(Math.Min(1, h + parserInterval), f, dc, ap, bg);
            phaseRoundingAllowances = centralPhase.ToDictionary(pair => pair.Key,
                pair => Math.Max(Math.Abs(pair.Value - lowerPhase[pair.Key]), Math.Abs(pair.Value - upperPhase[pair.Key])));
            maximumPhaseRoundingAllowance = Math.Max(maximumPhaseRoundingAllowance, phaseRoundingAllowances.Values.Max());
            double u = Math.Sqrt(1 - f) * Math.Sqrt(1 - h) - Math.Sqrt(f) * Math.Sqrt(h);
            double v = Math.Sqrt(f) * Math.Sqrt(1 - h) + Math.Sqrt(1 - f) * Math.Sqrt(h);
            double previous = dc * Math.Pow(1 - 2 * f, 2) * h * h + 4 * ap * h * (1 - h);
            double selected = dc * Math.Pow((1 - 2 * f) * h + 2 * Math.Sqrt(f * (1 - f) * h * (1 - h)), 2);
            double reduction = (previous - selected) / (bg + previous);
            Close(u * u + v * v, 1, "norm preservation");
            Require(reduction >= -1e-12, "selected phase contact minimum");
            var result = new Result
            {
                ["galaxy"] = row["galaxy"], ["index"] = N(row, "index"), ["R_kpc"] = radius,
                ["is_outer"] = B(row, "is_outer"), ["reliable_scalar_derivative"] = B(row, "reliable_scalar_derivative"),
                ["active_envelope"] = B(galaxies[row["galaxy"]], "envelope_active"),
                ["scalar_potential_kms2"] = potential, ["scalar_potential_derivative"] = derivative,
                ["required_tangential_q_kms2"] = q, ["equilibrium_x"] = x,
                ["x_numerical_error"] = N(row, "derivative_absolute_error_kms2") / v2,
                ["pointwise_equilibrium_mass_eV"] = x > 0 ? mass * Math.Sqrt(x) : double.NaN,
                ["no_finite_mass_stationarity"] = x <= 0,
                ["fraction_independent_relative_speed_lower_kms"] = B(row, "positive_support") ? 2 * Math.Sqrt(Math.Max(q, 0)) : double.NaN,
                ["optimized_component1_over_sqrt_density"] = u, ["optimized_component2_over_sqrt_density"] = v,
                ["optimized_component2_fraction"] = v * v,
                ["quadrature_excess_contact_normalized"] = previous, ["optimized_excess_contact_normalized"] = selected,
                ["fractional_total_contact_reduction"] = reduction
            };
            var expected = expectedPoints[(row["galaxy"], row["index"])];
            foreach (string field in result.Keys.Where(k => expected.ContainsKey(k) && result[k] is double))
                Close(Value(result, field), N(expected, field), row["galaxy"] + ":" + field);
            foreach (string field in result.Keys.Where(k => expected.ContainsKey(k) && result[k] is bool))
                Require(Flag(result, field) == B(expected, field), field);
            points.Add(result);
        }
        phaseRoundingAllowances.Clear();
        Require(points.Count == 3034 && points.Count(p => Flag(p, "is_outer")) == 659, "complete inverse sample");
        var summaries = new List<Result>();
        foreach (var group in points.GroupBy(p => (string)p["galaxy"]!).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            var ordered = group.OrderBy(p => Value(p, "R_kpc")).ToArray();
            var outer = ordered.Where(p => Flag(p, "is_outer") && Flag(p, "reliable_scalar_derivative")).ToArray();
            var solution = Minimax(outer.Select(p => Value(p, "equilibrium_x")).ToArray());
            var summary = new Result
            {
                ["galaxy"] = group.Key, ["n_outer"] = outer.Length,
                ["active_envelope"] = Flag(ordered[0], "active_envelope"),
                ["x_min"] = solution.Min, ["x_max"] = solution.Max, ["c"] = solution.Scale,
                ["error"] = solution.Error, ["finite_mass"] = solution.Scale > 0,
                ["mass_eV"] = solution.Scale > 0 ? mass / Math.Sqrt(solution.Scale) : double.NaN,
                ["wrong_sign_rows"] = outer.Count(p => Flag(p, "no_finite_mass_stationarity")),
                ["observed_interval_node_crossings"] = ordered.Zip(ordered.Skip(1)).Count(pair => Value(pair.First, "optimized_component1_over_sqrt_density") * Value(pair.Second, "optimized_component1_over_sqrt_density") < 0),
                ["median_contact_reduction"] = Data.Median(outer.Select(p => Value(p, "fractional_total_contact_reduction"))),
                ["median_relative_speed_lower_kms"] = Data.Median(outer.Select(p => Value(p, "fraction_independent_relative_speed_lower_kms")))
            };
            var expected = expectedGalaxies[group.Key];
            foreach (var field in summary.Keys.Where(k => summary[k] is double or int))
                Close(Value(summary, field), N(expected, field), group.Key + ":" + field);
            summaries.Add(summary);
        }
        var outerPoints = points.Where(p => Flag(p, "is_outer") && Flag(p, "reliable_scalar_derivative")).ToArray();
        var wrong = outerPoints.Where(p => Flag(p, "no_finite_mass_stationarity")).ToArray();
        var closure = new Result
        {
            ["galaxies"] = summaries.Count, ["outer_rows"] = outerPoints.Length,
            ["active_outer_rows"] = outerPoints.Count(p => Flag(p, "active_envelope")),
            ["wrong_sign_rows"] = wrong.Length,
            ["wrong_sign_galaxies"] = wrong.Select(p => p["galaxy"]).Distinct().Count(),
            ["numerically_robust_wrong_sign_rows"] = wrong.Count(p => Value(p, "equilibrium_x") + 10 * Value(p, "x_numerical_error") < 0),
            ["finite_mass_per_galaxy_count"] = summaries.Count(p => Flag(p, "finite_mass")),
            ["per_galaxy_minimax_median"] = Data.Median(summaries.Select(p => Value(p, "error"))),
            ["active_per_galaxy_minimax_median"] = Data.Median(summaries.Where(p => Flag(p, "active_envelope")).Select(p => Value(p, "error"))),
            ["optimized_node_crossing_galaxies"] = summaries.Count(p => Value(p, "observed_interval_node_crossings") > 0),
            ["optimized_observed_interval_node_crossings"] = summaries.Sum(p => Value(p, "observed_interval_node_crossings")),
            ["median_outer_galaxy_contact_reduction"] = Data.Median(summaries.Where(p => Flag(p, "active_envelope")).Select(p => Value(p, "median_contact_reduction"))),
            ["median_galaxy_relative_speed_lower_kms"] = Data.Median(summaries.Select(p => Value(p, "median_relative_speed_lower_kms")))
        };
        foreach (var (key, value) in closure) Close(Convert.ToDouble(value), expectedSummary.GetProperty(key).GetDouble(), "summary:" + key);
        var parent = ParentCompatibility(source, output, Close, Require);
        var persistence = Persistence(source, output, Require);
        // Independent analytical boundary examples, including the infinite-mass infimum.
        Require(Minimax([1, 3]) == (1, 3, .5, .5), "positive minimax endpoints");
        Require(Minimax([-1, 3]).Error == 1 && Minimax([-1, 3]).Scale == 0, "wrong-sign minimax");
        Csv.Write(Path.Combine(output, "point_closure.csv"), points);
        Csv.Write(Path.Combine(output, "galaxy_closure.csv"), summaries);
        var report = new Result
        {
            ["passed"] = true, ["checks"] = checks, ["maximum_tolerance_fraction"] = maxScaledError,
            ["maximum_phase_csv_rounding_allowance"] = maximumPhaseRoundingAllowance,
            ["phase_csv_rounding_scope"] = "Only six phase fields: exact propagation of h +/- four unit-scale ulps through the analytic square-root formulas, for archived pandas parsing near h=1.",
            ["closure"] = closure, ["parent_compatibility"] = parent, ["persistence"] = persistence,
            ["scope"] = "Native algebra from retained density derivatives and derivative-error/domain certificates. No new core solve, numerical differentiation, fit, time evolution, or significance. Persistence audits summarize retained accepted records; full original producers are in reference/."
        };
        Data.SaveJson(Path.Combine(output, "report.json"), report);
        return report;
    }

    static Dictionary<string, double> Phase(double h, double f, double dc, double ap, double bg)
    {
        double u = Math.Sqrt(1 - f) * Math.Sqrt(1 - h) - Math.Sqrt(f) * Math.Sqrt(h);
        double v = Math.Sqrt(f) * Math.Sqrt(1 - h) + Math.Sqrt(1 - f) * Math.Sqrt(h);
        double previous = dc * Math.Pow(1 - 2 * f, 2) * h * h + 4 * ap * h * (1 - h);
        double selected = dc * Math.Pow((1 - 2 * f) * h + 2 * Math.Sqrt(f * (1 - f) * h * (1 - h)), 2);
        return new Dictionary<string, double>
        {
            ["optimized_component1_over_sqrt_density"] = u,
            ["optimized_component2_over_sqrt_density"] = v,
            ["optimized_component2_fraction"] = v * v,
            ["quadrature_excess_contact_normalized"] = previous,
            ["optimized_excess_contact_normalized"] = selected,
            ["fractional_total_contact_reduction"] = (previous - selected) / (bg + previous)
        };
    }

    static (double Min, double Max, double Scale, double Error) Minimax(double[] x)
    {
        if (x.Length == 0 || x.Any(v => !double.IsFinite(v))) throw new InvalidDataException("Invalid minimax sample.");
        double min = x.Min(), max = x.Max();
        return min <= 0 ? (min, max, 0, 1) : (min, max, 2 / (min + max), (max - min) / (max + min));
    }

    static Result ParentCompatibility(string source, string output, Action<double, double, string> close, Action<bool, string> require)
    {
        var rows = Csv.Read(Path.Combine(source, "stationary_compatibility_rows.csv"));
        var expected = Data.Json(Path.Combine(source, "expected", "compatibility_summary.json"));
        double r = .5 * (Math.Tanh(.5) / .5 + 1 / Math.Pow(Math.Cosh(.5), 2));
        double f = (1 + 6 * r) / (1 + 11 * r);
        double d = Math.Pow(Math.Cosh(.5), 2) * (1 - r) * (1 + 11 * r);
        double coefficient = 2 * Math.PI * G * d / (C * C);
        close(f, expected.GetProperty("benchmark").GetProperty("f_star").GetDouble(), "parent f_star");
        close(d, expected.GetProperty("benchmark").GetProperty("D").GetDouble(), "parent D");
        require(Math.Abs(coefficient / expected.GetProperty("benchmark").GetProperty("coefficient_2piGD_over_c2").GetDouble() - 1) < 1e-12, "parent coefficient");
        var results = new List<Result>();
        foreach (var row in rows)
        {
            double a = Math.Sqrt(N(row, "rho_core")), b = Math.Sqrt(N(row, "rho_envelope"));
            double u = Math.Sqrt(1 - f) * a - Math.Sqrt(f) * b, v = Math.Sqrt(f) * a + Math.Sqrt(1 - f) * b;
            double ab = N(row, "amplitude_laplacian_term_AB"), ba = N(row, "amplitude_laplacian_term_BA");
            double lhs = ab - ba, rhs = coefficient * ((1 - 2 * f) * b * b + 2 * Math.Sqrt(f * (1 - f)) * a * b) * u * v;
            double error = Math.Max(Math.Abs(lhs - N(row, "fd_lhs_step1")), Math.Max(Math.Abs(lhs - N(row, "fd_lhs_step2")), Math.Abs(N(row, "fd_lhs_step1") - N(row, "fd_lhs_step2"))));
            double ratio = N(row, "rho_envelope") / N(row, "rho_core");
            bool positive = ratio < (2.0 / 3) * (1 - 1e-10) || ratio > 35 * (1 + 1e-10);
            bool negative = ratio > (5.0 / 7) * (1 + 1e-10) && ratio < 24 * (1 - 1e-10);
            bool resolved = Math.Abs(lhs) > 10 * error + 1e-10 * (Math.Abs(ab) + Math.Abs(ba));
            bool obstruction = B(row, "valid_archive") && resolved && (lhs < 0 && positive || lhs > 0 && negative);
            close(lhs, N(row, "lhs"), "parent lhs"); close(rhs, N(row, "rhs_planck"), "parent rhs");
            close(error, N(row, "fd_abs_error"), "parent derivative error");
            require(obstruction == B(row, "uniform_branch_obstruction"), "uniform-branch obstruction");
            results.Add(new Result { ["galaxy"] = row["galaxy"], ["index"] = N(row, "index"),
                ["is_outer"] = B(row, "is_outer"), ["valid_archive"] = B(row, "valid_archive"),
                ["envelope_active"] = B(row, "envelope_active"), ["lhs"] = lhs, ["rhs_planck"] = rhs,
                ["fd_abs_error"] = error, ["uniform_branch_obstruction"] = obstruction });
        }
        var outer = results.Where(p => Flag(p, "is_outer") && Flag(p, "valid_archive") && Flag(p, "envelope_active")).ToArray();
        var blocked = outer.Where(p => Flag(p, "uniform_branch_obstruction")).ToArray();
        require(outer.Length == 488 && blocked.Length == 222 && blocked.Select(p => p["galaxy"]).Distinct().Count() == 56, "parent published counts");
        Csv.Write(Path.Combine(output, "parent_compatibility.csv"), results);
        return new Result { ["rows"] = rows.Count, ["trusted_active_outer_rows"] = outer.Length,
            ["branch_obstruction_rows"] = blocked.Length, ["branch_obstruction_galaxies"] = blocked.Select(p => p["galaxy"]).Distinct().Count(),
            ["scope"] = "Necessary equation recomputed from retained amplitude-Laplacian terms and finite-difference checks; no new derivative evaluation." };
    }

    static Result Persistence(string source, string output, Action<bool, string> require)
    {
        var audit = Data.Json(Path.Combine(source, "persistence_audit.json"));
        var runs = audit.GetProperty("runs").EnumerateArray().ToArray();
        require(audit.GetProperty("status").GetString() == "pass" && runs.Length == 14, "accepted persistence index");
        require(runs.Select(r => r.GetProperty("galaxy").GetString()).Distinct().Count() == 7, "seven persistence galaxies");
        require(runs.All(r => r.GetProperty("actual_horizon_Gyr").GetDouble() >= 9.999999), "persistence horizon");
        var records = runs.Select(r => new Result
        {
            ["galaxy"] = r.GetProperty("galaxy").GetString(), ["mode"] = r.GetProperty("mode").GetString(),
            ["tag"] = r.GetProperty("tag").GetString(), ["horizon_Gyr"] = r.GetProperty("actual_horizon_Gyr").GetDouble(),
            ["max_relative_norm_drift"] = r.GetProperty("max_relative_norm_drift").GetDouble(),
            ["max_relative_energy_drift"] = r.GetProperty("max_relative_energy_drift").GetDouble(),
            ["max_enclosed_mass_change_over_N"] = r.GetProperty("max_enclosed_mass_change_over_N").GetDouble()
        }).ToArray();
        Csv.Write(Path.Combine(output, "accepted_persistence.csv"), records);
        return new Result { ["accepted_runs"] = runs.Length, ["rejected_runs_retained"] = audit.GetProperty("rejected_numerical_runs").GetArrayLength(),
            ["scope"] = "Retained audit records only; no fresh evolution or stability certificate." };
    }
}
