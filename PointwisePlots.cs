using System.Globalization;
using System.Security;
using System.Text;
using System.Text.Json;
using Row = System.Collections.Generic.Dictionary<string, string>;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static class PointwisePlots
{
    const string Folder = "pointwise_comparison_2026_09_25";
    static readonly string[] Comparators = ["scalar_core", "compact_plummer", "nfw", "baryon_only"];
    static readonly Dictionary<string, string> Labels = new()
    {
        ["scalar_core"] = "Scalar core",
        ["extended"] = "Fixed envelope",
        ["compact_plummer"] = "Compact Plummer",
        ["nfw"] = "NFW",
        ["baryon_only"] = "Core-calibrated baryons"
    };
    static readonly Dictionary<string, string> Colors = new()
    {
        ["scalar_core"] = "#69717a",
        ["extended"] = "#176a8a",
        ["compact_plummer"] = "#5f8294",
        ["nfw"] = "#bf682a",
        ["baryon_only"] = "#343434"
    };
    static double N(Row row, string key) => Csv.Number(row, key);
    static string F(double value, string format = "0.###") => value.ToString(format, CultureInfo.InvariantCulture);

    public static string[] Run(string dataRoot, string outputRoot)
    {
        string source = Path.Combine(dataRoot, "numerics", Folder), result = Path.Combine(outputRoot, "pointwise_comparison", "results");
        string output = Path.Combine(outputRoot, "figures"), plotData = Path.Combine(outputRoot, "pointwise_comparison", "plots");
        Directory.CreateDirectory(output);
        Directory.CreateDirectory(plotData);
        string[] files = [Path.Combine(source, "data", "pointwise_outer.csv"), Path.Combine(result, "primary_cluster_results.csv"),
            Path.Combine(result, "score_sensitivity.csv"), Path.Combine(result, "pointwise_paired_differences.csv"), Path.Combine(result, "residual_distribution_summary.csv")];
        var hashes = files.ToDictionary(p => p, Data.FileSha256);
        var outer = Csv.Read(files[0]);
        var primary = Csv.Read(files[1]).ToDictionary(r => r["comparator"]);
        var sensitivity = Csv.Read(files[2]);
        var paired = Csv.Read(files[3]);
        var summary = Csv.Read(files[4]).ToDictionary(r => r["model"]);
        foreach (string model in Comparators)
            Require(N(primary[model], "n_galaxies") == 131 && N(primary[model], "n_points") == 659, "Primary cluster population");
        var checks = new Result { ["no_refit"] = true, ["statistics_loaded_from_regenerated_csv"] = true, ["no_new_p_values"] = true };
        string first = Path.Combine(output, "pointwise_cluster_comparison.svg"), radial = Path.Combine(output, "pointwise_radial_differences.svg"), covariance = Path.Combine(output, "pointwise_covariance_sensitivity.svg");
        Cluster(outer, primary, sensitivity, summary, first, plotData, checks);
        Radial(paired, radial, checks);
        Covariance(primary, sensitivity, covariance, plotData, checks);
        var comparisons = new Result();
        foreach (string name in new[] { "plotted_ecdf_values.csv", "plotted_forest_values.csv", "plotted_covariance_values.csv" })
            comparisons[name] = CompareCsv(Path.Combine(source, "plots", name), Path.Combine(plotData, name));
        using var reference = JsonDocument.Parse(File.ReadAllText(Path.Combine(source, "plots", "plot_checks.json")));
        using var actual = JsonDocument.Parse(JsonSerializer.Serialize(checks));
        int compared = 0;
        foreach (string section in new[] { "ecdf", "forest", "radial", "covariance" })
            CompareJson(reference.RootElement.GetProperty(section), actual.RootElement.GetProperty(section), section, ref compared);
        checks["reference_plot_tables"] = comparisons;
        checks["reference_plot_check_fields"] = compared;
        checks["inputs_unchanged"] = files.All(p => hashes[p] == Data.FileSha256(p));
        Require((bool)checks["inputs_unchanged"]!, "Plotting must not change inputs");
        checks["input_manifest"] = hashes.Select(p => new { path = p.Key, sha256 = p.Value }).ToArray();
        checks["passed"] = true;
        Data.SaveJson(Path.Combine(plotData, "plot_checks.json"), checks);
        return [first, radial, covariance];
    }

    static void Cluster(List<Row> outer, Dictionary<string, Row> primary, List<Row> sensitivity, Dictionary<string, Row> summary,
        string path, string plotData, Result checks)
    {
        var svg = new Svg(1700, 1030, "Current outer-point comparison", "659 original measurements in 131 galaxies · frozen predictions · 25 September 2026");
        var ecdfs = new Dictionary<string, Result>();
        var forestChecks = new Dictionary<string, double[]>();
        var ecdfRows = new List<Result>();
        var forestRows = new List<Result>();
        var sidecar = new List<Result>();
        var ecdfData = new Dictionary<string, (double[] X, double[] Y)>();
        double maximum = 0;
        foreach (string model in new[] { "scalar_core", "extended", "nfw" })
        {
            var values = BuildWeightedEcdf(outer, model, summary[model], ecdfRows);
            maximum = Math.Max(maximum, values.X[^1]);
            ecdfData[model] = (values.X, values.Y);
            ecdfs[model] = values.Checks;
        }
        var ecdf = new Axis(svg, 110, 125, 590, 430, 0, maximum * 1.08, 0, 1.02,
            "(a) Every original outer observation", "Absolute standardized residual |z|", "Equal-galaxy cumulative fraction",
            [0, 1, 3, 10, 30], [0, .2, .4, .6, .8, 1], x => SymLog(x, 3));
        int series = 0;
        foreach (var (model, values) in ecdfData)
        {
            svg.Step(values.X.Select(ecdf.X).ToArray(), values.Y.Select(ecdf.Y).ToArray(), Colors[model], model == "nfw" ? "7,5" : "");
            svg.Line(165 + series * 195, 630, 193 + series * 195, 630, Colors[model], model == "nfw" ? "7,5" : "", 2);
            svg.Text(200 + series * 195, 635, Labels[model], 15);
            for (int i = 0; i < values.X.Length; i++)
                sidecar.Add(Series(0, model, i, values.X[i], values.Y[i]));
            series++;
        }
        svg.Text(126, 150, "131 galaxies / 659 points", 15);
        svg.Text(110, 665, "Full tail retained; x is linear to 3, then logarithmic.", 15, color: "#555");
        var forest = new Axis(svg, 985, 125, 615, 250, -12, 32, -.55, 2.55,
            "(b) Mean paired loss differences", "Envelope minus comparator", "", [-10, 0, 10, 20, 30], [0, 1, 2],
            yLabels: new()
            {
                [0] = "NFW",
                [1] = "Compact Plummer",
                [2] = "Scalar core"
            });
        forest.ZeroX();
        var baryon = new Axis(svg, 985, 490, 615, 95, -640, 35, -.6, .6,
            "Core-calibrated baryon ablation: separate horizontal scale", "Mean squared-standardized loss difference", "", [-600, -400, -200, 0], [0],
            yLabels: new()
            {
                [0] = "Baryons"
            });
        baryon.ZeroX();
        foreach (string model in Comparators)
        {
            Row row = primary[model];
            double effect = N(row, "effect"), low = N(row, "bootstrap_t_ci_low"), high = N(row, "bootstrap_t_ci_high");
            Require(low <= effect && effect <= high, "Forest interval ordering");
            Axis ax = model == "baryon_only" ? baryon : forest;
            double y = model == "scalar_core" ? 2 : model == "compact_plummer" ? 1 : 0;
            ax.HorizontalInterval(effect, y, low, high, Colors[model]);
            forestChecks[model] = [low, effect, high];
            forestRows.Add(new()
            {
                ["comparator"] = model,
                ["effect"] = effect,
                ["bootstrap_t_ci_low"] = low,
                ["bootstrap_t_ci_high"] = high
            });
            var entry = Series(model == "baryon_only" ? 2 : 1, model, 0, effect, y);
            entry["x_lower"] = low;
            entry["x_upper"] = high;
            sidecar.Add(entry);
        }
        svg.Text(810, 665, "Negative favors envelope. Bars: marginal 95% galaxy-bootstrap-t intervals.", 15, color: "#555");
        svg.Rect(85, 700, 1530, 276, "#f3f6f8");
        Row mae = sensitivity.Single(r => r["score"] == "absolute_kms" && r["comparator"] == "scalar_core");
        double coreMae = N(summary["scalar_core"], "mean_abs_kms_equal_galaxy"), envelopeMae = N(summary["extended"], "mean_abs_kms_equal_galaxy");
        Close(envelopeMae - coreMae, N(mae, "effect"), "Mean absolute speed-error contrast");
        svg.Text(110, 736, $"Equal-galaxy outer MAE: scalar core {F(coreMae, "0.00")} → envelope {F(envelopeMae, "0.00")} km/s.", 19);
        svg.Text(110, 764, $"Reduction {F(-N(mae, "effect"), "0.00")} km/s; marginal 95% bootstrap-t interval [{F(-N(mae, "bootstrap_t_ci_high"), "0.00")}, {F(-N(mae, "bootstrap_t_ci_low"), "0.00")}].", 17);
        var mc = new Dictionary<string, Result>();
        for (int i = 0; i < Comparators.Length; i++)
        {
            string model = Comparators[i];
            Row r = primary[model];
            double x = 110 + i * 375;
            bool floor = Csv.Flag(r, "monte_carlo_floor");
            if (floor)
            {
                Require(N(r, "wild_exceedances") == 0, "Monte Carlo floor has zero exceedances");
                Close(N(r, "p_wild"), 1 / (N(r, "wild_draws") + 1), "Plus-one Monte Carlo floor");
            }
            svg.Text(x, 807, Labels[model], 16);
            svg.Text(x, 833, $"Holm p = {F(N(r, "p_holm4"), "0.#####")}; |Z| = {F(N(r, "Z_holm4"), "0.00")}", 16);
            svg.Text(x, 858, floor ? "Monte Carlo resolution reached" : N(r, "p_holm4") < .05 ? "Holm-adjusted test below 0.05" : "No rejection at the 5% threshold", 14, color: floor ? "#815c24" : "#555");
            svg.Text(x, 880, $"{F(N(r, "wild_exceedances"), "N0")} / {F(N(r, "wild_draws"), "N0")} wild-test exceedances", 13, color: "#555");
            mc[model] = new()
            {
                ["monte_carlo_floor"] = floor,
                ["wild_exceedances"] = N(r, "wild_exceedances"),
                ["wild_draws"] = N(r, "wild_draws"),
                ["p_holm4"] = N(r, "p_holm4"),
                ["Z_holm4"] = N(r, "Z_holm4"),
                ["raw_MC_95low"] = N(r, "raw_MC_95low"),
                ["raw_MC_95high"] = N(r, "raw_MC_95high")
            };
        }
        svg.Text(110, 918, "Tests: null-imposed wild bootstrap, four-test Holm adjustment. Intervals are marginal, not simultaneous or test inversions.", 15);
        svg.Text(110, 945, "At zero exceedances, p uses the plus-one convention; Z is its Gaussian equivalent, not a measured significance or bound.", 15);
        svg.Text(850, 1009, "Conditional on frozen predictions; no refitting or new significance tests in these figures.", 14, "middle", "#555");
        checks["ecdf"] = ecdfs;
        checks["forest"] = forestChecks;
        checks["monte_carlo_reporting"] = mc;
        checks["absolute_speed_error"] = new
        {
            core_mae = coreMae,
            envelope_mae = envelopeMae,
            reduction = -N(mae, "effect"),
            ci_low = -N(mae, "bootstrap_t_ci_high"),
            ci_high = -N(mae, "bootstrap_t_ci_low")
        };
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), sidecar);
        Csv.Write(Path.Combine(plotData, "plotted_ecdf_values.csv"), ecdfRows);
        Csv.Write(Path.Combine(plotData, "plotted_forest_values.csv"), forestRows);
    }

    static (double[] X, double[] Y, Result Checks) BuildWeightedEcdf(
        List<Row> outer, string model, Row summary, List<Result> exportedRows)
    {
        const int galaxyCount = 131;
        const int pointCount = 659;
        Row[] points = outer.Where(r => r["model"] == model)
            .OrderBy(r => Math.Abs(N(r, "standardized_residual")))
            .ThenBy(r => r["galaxy"], StringComparer.Ordinal)
            .ThenBy(r => N(r, "point_index"))
            .ToArray();
        Require(points.Length == pointCount && points.Select(r => r["galaxy"]).Distinct().Count() == galaxyCount
            && points.All(r => Csv.Flag(r, "valid_prediction")), "Every original finite outer observation");
        var counts = points.GroupBy(r => r["galaxy"]).ToDictionary(group => group.Key, group => group.Count());
        var x = new double[pointCount + 1];
        var y = new double[pointCount + 1];
        var weights = new double[pointCount];
        double mean = 0;
        for (int i = 0; i < points.Length; i++)
        {
            x[i + 1] = Math.Abs(N(points[i], "standardized_residual"));
            // Each galaxy contributes the same total weight, regardless of its radial point count.
            weights[i] = 1.0 / (galaxyCount * counts[points[i]["galaxy"]]);
            y[i + 1] = y[i] + weights[i];
            mean += x[i + 1] * weights[i];
            exportedRows.Add(new()
            {
                ["galaxy"] = points[i]["galaxy"],
                ["point_index"] = N(points[i], "point_index"),
                ["model"] = model,
                ["abs_z"] = x[i + 1],
                ["weight"] = weights[i],
                ["cumulative_weight"] = y[i + 1]
            });
        }
        Require(Math.Abs(y[^1] - 1) < 2e-14, "Unit ECDF weight");
        var perGalaxyWeights = points.Select((row, i) => (Galaxy: row["galaxy"], Weight: weights[i])).GroupBy(row => row.Galaxy);
        foreach (var group in perGalaxyWeights)
            Require(Math.Abs(group.Sum(row => row.Weight) - 1.0 / galaxyCount) < 1e-15, "Equal total weight per galaxy");
        Close(mean, N(summary, "mean_abs_standardized_equal_galaxy"), "Weighted absolute residual");
        var checks = new Result
        {
            ["points"] = pointCount,
            ["galaxies"] = galaxyCount,
            ["total_weight"] = y[^1],
            ["maximum_abs_z"] = x[^1],
            ["weighted_mean_abs_z"] = mean
        };
        return (x, y, checks);
    }

    static void Radial(List<Row> paired, string path, Result checks)
    {
        var svg = new Svg(1500, 700, "Original outer-point loss contrasts", "All 659 measurements in each panel; observations within a galaxy are correlated");
        var audit = new Dictionary<string, Result>();
        var sidecar = new List<Result>();
        int panel = 0;
        foreach (string model in new[] { "scalar_core", "nfw" })
        {
            Row[] rows = paired.Where(r => r["comparator"] == model).ToArray();
            Require(rows.Length == 659 && rows.Select(r => r["galaxy"]).Distinct().Count() == 131 && rows.All(r => Csv.Flag(r, "comparison_valid")), "Radial panel population");
            double[] x = rows.Column("radius_fraction_of_last"), y = rows.Column("delta_squared");
            Require(x.All(v => v > 0 && v <= 1 + 1e-13) && y.All(double.IsFinite), "Finite original radial contrast");
            double lo = SymLog(y.Min(), 1), hi = SymLog(y.Max(), 1), pad = .07 * (hi - lo);
            double ymin = InvSymLog(lo - pad, 1), ymax = InvSymLog(hi + pad, 1);
            Require(y.All(v => v > ymin && v < ymax), "Every radial point lies inside displayed limits");
            var ax = new Axis(svg, 115 + panel * 740, 130, 560, 380, Math.Min(.25, x.Min() - .02), 1.025, ymin, ymax,
                $"({(panel == 0 ? "a" : "b")}) Envelope minus {Labels[model]}", "Radius / last measured radius", "Point loss difference zE² − zref²",
                [.25, .4, .6, .8, 1], SignedTicks(ymin, ymax), ty: v => SymLog(v, 1));
            for (int i = 0; i < x.Length; i++)
            {
                svg.Point(ax.X(x[i]), ax.Y(y[i]), Colors[model], 2.5, .5);
                sidecar.Add(Series(panel, model, i, x[i], y[i], rows[i]["galaxy"], N(rows[i], "point_index")));
            }
            ax.ZeroY();
            svg.Text(ax.Left + 14, 155, "659 points / 131 galaxies", 15);
            audit[model] = new()
            {
                ["points"] = 659,
                ["galaxies"] = 131,
                ["x_min"] = x.Min(),
                ["x_max"] = x.Max(),
                ["delta_min"] = y.Min(),
                ["delta_max"] = y.Max(),
                ["all_points_inside_axes"] = true
            };
            panel++;
        }
        svg.Text(110, 625, "Below zero favors the envelope. Vertical axes are linear from −1 to 1 and logarithmic outside.", 17);
        svg.Text(110, 662, "Descriptive, unbinned contrasts; no radial relation is fitted and no measurements are deleted.", 17);
        checks["radial"] = audit;
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), sidecar);
    }

    static void Covariance(Dictionary<string, Row> primary, List<Row> sensitivity, string path, string plotData, Result checks)
    {
        var svg = new Svg(1520, 1140, "Within-galaxy covariance sensitivity", "Assumed common correlation, not estimated; all 131 galaxies and 659 measurements retained");
        var audit = new Dictionary<string, Result>();
        var plotted = new List<Result>();
        var sidecar = new List<Result>();
        for (int panel = 0; panel < Comparators.Length; panel++)
        {
            string model = Comparators[panel];
            Row[] rows = [primary[model], .. sensitivity.Where(r => r["score"] == "correlated_quadratic" && r["comparator"] == model).OrderBy(r => N(r, "rho"))];
            Require(rows.Length == 4, "Four covariance cases");
            double[] x = [0, .25, .5, .75], y = rows.Column("effect"), low = rows.Column("bootstrap_t_ci_low"), high = rows.Column("bootstrap_t_ci_high");
            for (int i = 1; i < 4; i++)
                Close(N(rows[i], "rho"), x[i], "Declared covariance rho");
            Require(Enumerable.Range(0, 4).All(i => low[i] <= y[i] && y[i] <= high[i]), "Covariance interval ordering");
            double lo = Math.Min(low.Min(), 0), hi = Math.Max(high.Max(), 0), pad = .08 * (hi - lo), ymin = lo - pad, ymax = hi + pad;
            var ax = new Axis(svg, 125 + panel % 2 * 750, 135 + panel / 2 * 445, 550, 300, -.07, .82, ymin, ymax,
                $"({(char)('a' + panel)}) Envelope minus {Labels[model]}", "Assumed within-galaxy correlation ρ", "Mean correlated quadratic-loss difference", x, NiceTicks(ymin, ymax));
            ax.ZeroY();
            for (int i = 0; i < 4; i++)
            {
                ax.VerticalInterval(x[i], y[i], low[i], high[i], Colors[model]);
                plotted.Add(new()
                {
                    ["comparator"] = model,
                    ["rho"] = x[i],
                    ["effect"] = y[i],
                    ["bootstrap_t_ci_low"] = low[i],
                    ["bootstrap_t_ci_high"] = high[i]
                });
                var row = Series(panel, model, i, x[i], y[i]);
                row["y_lower"] = low[i];
                row["y_upper"] = high[i];
                sidecar.Add(row);
            }
            audit[model] = new()
            {
                ["rho"] = x,
                ["effects"] = y,
                ["low"] = low,
                ["high"] = high,
                ["points"] = 659,
                ["galaxies"] = 131
            };
        }
        svg.Text(115, 1040, "Bars: marginal 95% whole-galaxy bootstrap-t intervals. Panels use separate vertical scales.", 17);
        svg.Text(115, 1080, "The common correlation is assumed. These are covariance sensitivities, not additional significance tests.", 17);
        checks["covariance"] = audit;
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), sidecar);
        Csv.Write(Path.Combine(plotData, "plotted_covariance_values.csv"), plotted);
    }

    static Result Series(int panel, string model, int index, double x, double y, string galaxy = "", double? pointIndex = null) =>
        new()
        {
            ["panel"] = panel,
            ["label"] = Labels[model],
            ["model"] = model,
            ["index"] = index,
            ["x"] = x,
            ["y"] = y,
            ["x_lower"] = null,
            ["x_upper"] = null,
            ["y_lower"] = null,
            ["y_upper"] = null,
            ["galaxy"] = galaxy,
            ["point_index"] = pointIndex
        };
    static double SymLog(double x, double threshold) => Math.Abs(x) <= threshold ? x / threshold * .7 / .9 : Math.Sign(x) * (.7 / .9 + Math.Log10(Math.Abs(x) / threshold));
    static double InvSymLog(double x, double threshold) => Math.Abs(x) <= .7 / .9 ? x * threshold * .9 / .7 : Math.Sign(x) * threshold * Math.Pow(10, Math.Abs(x) - .7 / .9);
    static double[] SignedTicks(double min, double max) => new[] { -10000d, -1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000 }.Where(v => v >= min && v <= max).ToArray();
    static double[] NiceTicks(double min, double max)
    {
        double raw = (max - min) / 5, power = Math.Pow(10, Math.Floor(Math.Log10(raw)));
        double step = new[] { 1d, 2, 2.5, 5, 10 }.First(v => v * power >= raw) * power;
        return Enumerable.Range(0, 12).Select(i => Math.Ceiling(min / step) * step + i * step).Where(v => v <= max).ToArray();
    }
    static void Require(bool condition, string label)
    {
        if (!condition)
            throw new InvalidDataException("Pointwise plot check: " + label);
    }
    static void Close(double a, double b, string label) => Require(double.IsFinite(a) && double.IsFinite(b) && Math.Abs(a - b) <= 2e-10 + 2e-11 * Math.Abs(b), label);
    static Result CompareCsv(string reference, string generated)
    {
        var expectedRows = Csv.Read(reference);
        var actualRows = Csv.Read(generated);
        Require(expectedRows.Count == actualRows.Count, "Reference plotted row count");
        int fields = 0;
        for (int i = 0; i < expectedRows.Count; i++)
        {
            foreach (var (key, expected) in expectedRows[i])
            {
                string actual = actualRows[i][key];
                if (double.TryParse(expected, NumberStyles.Float, CultureInfo.InvariantCulture, out double expectedValue)
                    && double.TryParse(actual, NumberStyles.Float, CultureInfo.InvariantCulture, out double actualValue))
                    Close(expectedValue, actualValue, reference + ":" + key);
                else
                    Require(expected == actual, "Reference plotted categorical value " + key);
                fields++;
            }
        }
        return new()
        {
            ["passed"] = true,
            ["rows"] = expectedRows.Count,
            ["fields"] = fields
        };
    }

    static void CompareJson(JsonElement expected, JsonElement actual, string label, ref int count)
    {
        switch (expected.ValueKind)
        {
            case JsonValueKind.Object:
                foreach (var item in expected.EnumerateObject())
                    CompareJson(item.Value, actual.GetProperty(item.Name), label + "." + item.Name, ref count);
                break;
            case JsonValueKind.Array:
                Require(expected.GetArrayLength() == actual.GetArrayLength(), label + " length");
                for (int i = 0; i < expected.GetArrayLength(); i++)
                    CompareJson(expected[i], actual[i], label + "[]", ref count);
                break;
            case JsonValueKind.Number:
                Close(expected.GetDouble(), actual.GetDouble(), label);
                count++;
                break;
            default:
                Require(expected.ToString() == actual.ToString(), label);
                count++;
                break;
        }
    }

    sealed class Axis
    {
        readonly Svg svg; readonly double top, width, height, xmin, xmax, ymin, ymax; readonly Func<double, double> tx, ty;
        public double Left
        {
            get;
        }
        public Axis(Svg svg, double left, double top, double width, double height, double xmin, double xmax, double ymin, double ymax,
            string title, string xlabel, string ylabel, double[] xticks, double[] yticks, Func<double, double>? tx = null, Func<double, double>? ty = null, Dictionary<double, string>? yLabels = null)
        {
            this.svg = svg;
            Left = left;
            this.top = top;
            this.width = width;
            this.height = height;
            this.tx = tx ?? (v => v);
            this.ty = ty ?? (v => v);
            this.xmin = this.tx(xmin);
            this.xmax = this.tx(xmax);
            this.ymin = this.ty(ymin);
            this.ymax = this.ty(ymax);
            svg.Text(left, top - 28, title, 19);
            foreach (double tick in xticks.Where(v => v >= xmin && v <= xmax))
            {
                svg.Line(X(tick), top, X(tick), top + height, "#e6e6e6");
                svg.Text(X(tick), top + height + 26, F(tick), 14, "middle");
            }
            foreach (double tick in yticks.Where(v => v >= ymin && v <= ymax))
            {
                svg.Line(left, Y(tick), left + width, Y(tick), "#e6e6e6");
                svg.Text(left - 14, Y(tick) + 5, yLabels?.GetValueOrDefault(tick) ?? F(tick), 15, "end");
            }
            svg.Line(left, top, left, top + height, "#444");
            svg.Line(left, top + height, left + width, top + height, "#444");
            svg.Text(left + width / 2, top + height + 61, xlabel, 17, "middle");
            if (ylabel.Length > 0)
                svg.Rotated(left - 76, top + height / 2, ylabel, 16);
        }
        public double X(double v) => Left + (tx(v) - xmin) / (xmax - xmin) * width;
        public double Y(double v) => top + (ymax - ty(v)) / (ymax - ymin) * height;
        public void ZeroX() => svg.Line(X(0), top, X(0), top + height, "#777");
        public void ZeroY() => svg.Line(Left, Y(0), Left + width, Y(0), "#555", width: 1.4);
        public void HorizontalInterval(double x, double y, double low, double high, string color)
        {
            svg.Line(X(low), Y(y), X(high), Y(y), color, width: 2);
            svg.Line(X(low), Y(y) - 5, X(low), Y(y) + 5, color, width: 2);
            svg.Line(X(high), Y(y) - 5, X(high), Y(y) + 5, color, width: 2);
            svg.Point(X(x), Y(y), color, 5);
        }
        public void VerticalInterval(double x, double y, double low, double high, string color)
        {
            svg.Line(X(x), Y(low), X(x), Y(high), color, width: 2);
            svg.Line(X(x) - 5, Y(low), X(x) + 5, Y(low), color, width: 2);
            svg.Line(X(x) - 5, Y(high), X(x) + 5, Y(high), color, width: 2);
            svg.Point(X(x), Y(y), color, 5);
        }
    }
    sealed class Svg
    {
        readonly StringBuilder content;
        static string E(string text) => SecurityElement.Escape(text) ?? "";
        public Svg(int width, int height, string title, string subtitle)
        {
            content = new($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><title>{E(title)}</title><rect width=\"100%\" height=\"100%\" fill=\"white\"/><style>text{{font-family:Arial,sans-serif}}</style>");
            Text(width / 2d, 33, title, 24, "middle");
            Text(width / 2d, 65, subtitle, 16, "middle", "#555");
        }
        public void Text(double x, double y, string text, int size, string anchor = "start", string color = "#222") => content.Append($"<text x=\"{F(x)}\" y=\"{F(y)}\" font-size=\"{size}\" text-anchor=\"{anchor}\" fill=\"{color}\">{E(text)}</text>");
        public void Rotated(double x, double y, string text, int size) => content.Append($"<text transform=\"translate({F(x)},{F(y)}) rotate(-90)\" font-size=\"{size}\" text-anchor=\"middle\" fill=\"#222\">{E(text)}</text>");
        public void Line(double x1, double y1, double x2, double y2, string color, string dash = "", double width = 1) => content.Append($"<line x1=\"{F(x1)}\" y1=\"{F(y1)}\" x2=\"{F(x2)}\" y2=\"{F(y2)}\" stroke=\"{color}\" stroke-width=\"{F(width)}\" stroke-dasharray=\"{dash}\"/>");
        public void Rect(double x, double y, double w, double h, string color) => content.Append($"<rect x=\"{F(x)}\" y=\"{F(y)}\" width=\"{F(w)}\" height=\"{F(h)}\" fill=\"{color}\"/>");
        public void Point(double x, double y, string color, double r, double opacity = 1) => content.Append($"<circle cx=\"{F(x)}\" cy=\"{F(y)}\" r=\"{F(r)}\" fill=\"{color}\" opacity=\"{F(opacity)}\"/>");
        public void Step(double[] x, double[] y, string color, string dash)
        {
            var path = new StringBuilder($"M{F(x[0])},{F(y[0])}");
            for (int i = 1; i < x.Length; i++)
                path.Append($"L{F(x[i])},{F(y[i - 1])}L{F(x[i])},{F(y[i])}");
            content.Append($"<path d=\"{path}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.3\" stroke-dasharray=\"{dash}\"/>");
        }
        public void Save(string path) => File.WriteAllText(path, content.ToString() + "</svg>");
    }
}


