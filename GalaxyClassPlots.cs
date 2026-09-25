using System.Globalization;
using System.Security;
using System.Text;

namespace DarkUniverse;

public static class GalaxyClassPlots
{
    static readonly (string Classification, int Id, string Name, int Y)[] Classes =
        new[] { "S0", "Sa", "Sab", "Sb", "Sbc", "Sc", "Scd", "Sd", "Sdm", "Sm", "Im", "BCD" }
        .Select((n, i) => ("detailed", i, n, i)).Concat(new[] { ("broad", 0, "S0", 13), ("broad", 1, "Spirals", 14), ("broad", 2, "Im/BCD", 15) }).ToArray();
    static readonly (string Id, string Label, string Color, string Marker, bool Filled, double Offset)[] Models =
    [ ("A_core", "A core", "#24658c", "circle", false, -.16), ("A_nfw", "A NFW", "#be5e28", "triangle", false, -.08),
      ("B_core", "B core", "#24658c", "circle", true, 0), ("B_nfw", "B/C NFW", "#be5e28", "triangle", true, .08),
      ("C_core", "C core", "#147e83", "square", true, .16) ];
    static double N(Dictionary<string, string> row, string key) => Csv.Number(row, key);
    static Dictionary<string, string> One(List<Dictionary<string, string>> table, string classification, int id, params (string Key, string Value)[] where) =>
        table.Single(r => r["classification"] == classification && N(r, "class_id") == id && where.All(w => r[w.Key] == w.Value));

    public static string[] Run(string root, string outputRoot)
    {
        string input = Path.Combine(root, "numerics", "galaxy_classes"), output = Path.Combine(outputRoot, "figures");
        Directory.CreateDirectory(output);
        var comparisons = Csv.Read(Path.Combine(input, "class_comparisons.csv"));
        var scores = Csv.Read(Path.Combine(input, "class_model_scores.csv"));
        var residuals = Csv.Read(Path.Combine(input, "residual_regions.csv"));
        if (comparisons.Count != 60)
            throw new InvalidDataException("Galaxy classes require the declared 60-test family.");
        foreach (var b in scores.Where(r => r["model_id"] == "B_nfw"))
        {
            var c = One(scores, b["classification"], (int)N(b, "class_id"), ("model_id", "C_nfw"), ("fit_mode", b["fit_mode"]), ("region", b["region"]));
            foreach (string key in new[] { "n_galaxies", "n_failures", "median_loss", "mean_loss" })
                if (!N(b, key).Equals(N(c, key)))
                    throw new InvalidDataException("C NFW differs from B NFW: " + key);
        }
        DrawWins(comparisons, Path.Combine(output, "galaxy_class_wins.svg"));
        DrawLosses(scores, Path.Combine(output, "galaxy_class_losses.svg"));
        double limit = DrawResiduals(residuals, comparisons, Path.Combine(output, "galaxy_class_residuals.svg"));
        Data.SaveJson(Path.Combine(outputRoot, "galaxy_classes", "plot_audit.json"), new
        {
            class_count = 15,
            win_family_size = 60,
            win_panels = new[] { "B_core_vs_B_nfw", "C_core_vs_C_nfw", "C_core_vs_B_core" },
            intervals = "Marginal exact 95% Clopper-Pearson; not simultaneous",
            residual_color_limit = limit,
            significance_marker = "Filled if saved 60-family Holm p < 0.05; hollow otherwise",
            scope = "Reconstructed saved-table rendering; all twelve detailed and three broad classes retained"
        });
        return new[] { "galaxy_class_wins", "galaxy_class_losses", "galaxy_class_residuals" }.Select(name => Path.Combine(output, name + ".svg")).ToArray();
    }

    static void DrawWins(List<Dictionary<string, string>> table, string path)
    {
        var pairs = new[] { ("B_core_vs_B_nfw", "B core / NFW", "#24658c"), ("C_core_vs_C_nfw", "C core / NFW", "#147e83"), ("C_core_vs_B_core", "C core / B core", "#75528d") };
        var svg = new Svg(1310, 900, "Galaxy class prediction wins");
        var plotted = new List<Dictionary<string, object?>>();
        const double top = 105, dy = 42, left = 205, width = 325;
        svg.Text(24, 65, "Class (eligible n)", 15);
        for (int panel = 0; panel < pairs.Length; panel++)
        {
            var (pair, title, color) = pairs[panel];
            double ox = left + panel * 355;
            double X(double p) => ox + (p + .06) / 1.54 * width;
            svg.Text(ox + width / 2, 66, title, 18, "middle");
            svg.Rect(ox, top + 12.5 * dy, width, 3.1 * dy, "#edf0f3");
            svg.Line(ox, top + 12.5 * dy, ox + width, top + 12.5 * dy, "#999");
            foreach (double x in new[] { 0d, .5, 1d })
            {
                svg.Line(X(x), top - .6 * dy, X(x), top + 15.6 * dy, x == .5 ? "#aaa" : "#eee", x == .5 ? "4,4" : "");
                svg.Text(X(x), top + 16.1 * dy, x.ToString("0.#", CultureInfo.InvariantCulture), 14, "middle");
            }
            svg.Text(X(1.12), 89, "wins / n", 13);
            foreach (var cls in Classes)
            {
                var row = One(table, cls.Classification, cls.Id, ("pair_id", pair));
                double y = top + cls.Y * dy;
                int n = (int)N(row, "n_informative"), wins = (int)N(row, "left_wins");
                double p = N(row, "left_win_fraction"), low = N(row, "win_ci_low"), high = N(row, "win_ci_high");
                bool significant = double.IsFinite(N(row, "p_holm")) && N(row, "p_holm") < .05;
                if (panel == 0)
                    svg.Text(left - 14, y + 5, $"{cls.Name}  ({N(row, "n_total"):0})", 15, "end");
                if (n > 0)
                {
                    if (wins + N(row, "right_wins") != n || Math.Abs(p - (double)wins / n) > 1e-12 || !(0 <= low && low <= p && p <= high && high <= 1))
                        throw new InvalidDataException("Invalid class win interval.");
                    svg.Line(X(low), y, X(high), y, color);
                    svg.Line(X(low), y - 3, X(low), y + 3, color);
                    svg.Line(X(high), y - 3, X(high), y + 3, color);
                    svg.Marker(X(p), y, "circle", color, significant, 4.2);
                }
                else
                    svg.Text(X(.5), y + 5, "—", 16, "middle", "#888");
                svg.Text(X(1.12), y + 5, $"{wins}/{n}", 14);
                plotted.Add(new()
                {
                    ["panel"] = panel,
                    ["pair_id"] = pair,
                    ["classification"] = cls.Classification,
                    ["class_id"] = cls.Id,
                    ["class_name"] = cls.Name,
                    ["x"] = p,
                    ["y"] = cls.Y,
                    ["x_lower"] = low,
                    ["x_upper"] = high,
                    ["n_total"] = N(row, "n_total"),
                    ["n_informative"] = n,
                    ["left_wins"] = wins,
                    ["p_holm"] = N(row, "p_holm"),
                    ["filled"] = significant
                });
            }
        }
        svg.Text(740, 818, "Fraction won by the first-named core", 16, "middle");
        svg.Marker(478, 847, "circle", "#444", true, 4);
        svg.Text(491, 852, "Holm p < 0.05", 14);
        svg.Marker(710, 847, "circle", "#444", false, 4);
        svg.Text(723, 852, "Holm p ≥ 0.05", 14);
        svg.Text(740, 885, "Marginal exact 95% intervals · 60-test Holm family · retrospective", 14, "middle");
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), plotted);
    }

    static void DrawLosses(List<Dictionary<string, string>> table, string path)
    {
        var modes = new[] { ("full", "all", "All-radius fitted agreement"), ("inner", "outer", "Withheld outer prediction") };
        var finite = Classes.SelectMany(cls => Models.SelectMany(m => modes.Select(mode => N(One(table, cls.Classification, cls.Id,
            ("model_id", m.Id), ("fit_mode", mode.Item1), ("region", mode.Item2)), "median_loss")))).Where(v => double.IsFinite(v) && v > 0).ToArray();
        double lo = Math.Floor(Math.Log10(finite.Min()) - .06), hi = Math.Ceiling(Math.Log10(finite.Max()) + .06);
        var svg = new Svg(1330, 930, "Galaxy class descriptive losses");
        var plotted = new List<Dictionary<string, object?>>();
        const double top = 105, dy = 42, left = 230, width = 465;
        svg.Text(16, 64, "Class (full / outer n)", 15);
        for (int panel = 0; panel < modes.Length; panel++)
        {
            var (mode, region, title) = modes[panel];
            double ox = left + panel * 535;
            double X(double value) => ox + (Math.Log10(value) - lo) / (hi - lo) * width;
            svg.Text(ox + width / 2, 65, title, 18, "middle");
            svg.Rect(ox, top + 12.5 * dy, width, 3.1 * dy, "#edf0f3");
            svg.Line(ox, top + 12.5 * dy, ox + width, top + 12.5 * dy, "#999");
            for (int p = (int)lo; p <= hi; p++)
            {
                double x = X(Math.Pow(10, p));
                svg.Line(x, top - .6 * dy, x, top + 15.6 * dy, "#e8e8e8");
                svg.Text(x, top + 16.2 * dy, $"10^{p}", 13, "middle");
            }
            foreach (var cls in Classes)
            {
                if (panel == 0)
                {
                    var full = One(table, cls.Classification, cls.Id, ("model_id", "B_core"), ("fit_mode", "full"), ("region", "all"));
                    var outer = One(table, cls.Classification, cls.Id, ("model_id", "B_core"), ("fit_mode", "inner"), ("region", "outer"));
                    svg.Text(left - 16, top + cls.Y * dy + 5, $"{cls.Name}  ({N(full, "n_galaxies"):0}/{N(outer, "n_galaxies"):0})", 15, "end");
                }
                bool empty = true;
                foreach (var model in Models)
                {
                    var row = One(table, cls.Classification, cls.Id, ("model_id", model.Id), ("fit_mode", mode), ("region", region));
                    double value = N(row, "median_loss"), y = top + (cls.Y + model.Offset) * dy;
                    if (N(row, "n_galaxies") > 0)
                    {
                        empty = false;
                        if (double.IsFinite(value) && value > 0)
                            svg.Marker(X(value), y, model.Marker, model.Color, model.Filled, 4.2);
                        else
                            svg.Text(value == 0 ? ox : ox + width, y + 4, value == 0 ? "0" : "∞", 12, value == 0 ? "start" : "end", model.Color);
                    }
                    plotted.Add(new()
                    {
                        ["panel"] = panel,
                        ["classification"] = cls.Classification,
                        ["class_id"] = cls.Id,
                        ["class_name"] = cls.Name,
                        ["model_id"] = model.Id,
                        ["fit_mode"] = mode,
                        ["region"] = region,
                        ["x"] = value,
                        ["y"] = cls.Y + model.Offset,
                        ["n_galaxies"] = N(row, "n_galaxies"),
                        ["n_failures"] = N(row, "n_failures"),
                        ["median_loss"] = value
                    });
                }
                if (empty)
                    svg.Text(ox + width / 2, top + cls.Y * dy + 5, "—", 16, "middle", "#888");
            }
            svg.Text(ox + width / 2, 819, "Median galaxy loss χ²/N", 16, "middle");
        }
        for (int i = 0; i < Models.Length; i++)
        {
            var m = Models[i];
            double x = 300 + i * 180;
            svg.Marker(x, 858, m.Marker, m.Color, m.Filled, 4.5);
            svg.Text(x + 14, 863, m.Label, 14);
        }
        svg.Text(760, 906, "Galaxy-equal medians; memberships differ · NFW C = B · logarithmic scale", 14, "middle");
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), plotted);
    }

    static double DrawResiduals(List<Dictionary<string, string>> table, List<Dictionary<string, string>> comparisons, string path)
    {
        var models = new[] { ("B_core", "B core"), ("C_core", "C core"), ("B_nfw", "B/C NFW") };
        var regions = new[] { ("full", "inner_third", "Inner"), ("full", "middle_third", "Middle"), ("full", "outer_third", "Outer"), ("inner", "withheld_outer", "Held-out") };
        var plotted = new List<Dictionary<string, object?>>();
        double maximum = Classes.SelectMany(cls => models.SelectMany(m => regions.Select(r => One(table, cls.Classification, cls.Id,
            ("model_id", m.Item1), ("fit_mode", r.Item1), ("region", r.Item2))))).Where(r => N(r, "n_galaxies") > 0)
            .Select(r => N(r, "median_signed_standardized")).Where(double.IsFinite).Select(Math.Abs).Max();
        double limit = Math.Max(.5, Math.Ceiling(maximum * 2) / 2);
        var svg = new Svg(1350, 985, "Galaxy class residual patterns");
        const double top = 91, left = 222, cell = 82, dy = 42;
        svg.Text(18, 65, "Class (full / outer n)", 15);
        for (int panel = 0; panel < models.Length; panel++)
        {
            var (model, title) = models[panel];
            double ox = left + panel * 357;
            svg.Text(ox + cell * 2, 66, title, 18, "middle");
            foreach (var cls in Classes)
            {
                if (panel == 0)
                {
                    var full = One(table, cls.Classification, cls.Id, ("model_id", "B_core"), ("fit_mode", "full"), ("region", "inner_third"));
                    var outer = One(comparisons, cls.Classification, cls.Id, ("pair_id", "B_core_vs_B_nfw"));
                    svg.Text(left - 14, top + (cls.Y + .5) * dy + 5, $"{cls.Name} ({N(full, "n_galaxies"):0}/{N(outer, "n_total"):0})", 15, "end");
                }
                for (int col = 0; col < regions.Length; col++)
                {
                    var (mode, region, label) = regions[col];
                    var row = One(table, cls.Classification, cls.Id, ("model_id", model), ("fit_mode", mode), ("region", region));
                    double value = N(row, "n_galaxies") > 0 ? N(row, "median_signed_standardized") : double.NaN;
                    svg.Rect(ox + col * cell, top + cls.Y * dy, cell, dy, double.IsFinite(value) ? Diverging(value, limit) : "#f4f5f7", "white");
                    svg.Text(ox + (col + .5) * cell, top + (cls.Y + .5) * dy + 5,
                        double.IsFinite(value) ? (double.IsNegative(value) ? "-" : "+") + Math.Abs(value).ToString("0.0", CultureInfo.InvariantCulture) : "—", 15, "middle",
                        !double.IsFinite(value) ? "#999" : Math.Abs(value) > .59 * limit ? "white" : "#262626");
                    plotted.Add(new()
                    {
                        ["panel"] = panel,
                        ["classification"] = cls.Classification,
                        ["class_id"] = cls.Id,
                        ["class_name"] = cls.Name,
                        ["model_id"] = model,
                        ["fit_mode"] = mode,
                        ["region"] = region,
                        ["x"] = col,
                        ["y"] = cls.Y,
                        ["n_galaxies"] = N(row, "n_galaxies"),
                        ["median_signed_standardized"] = value,
                        ["color_limit"] = limit
                    });
                }
            }
            svg.Line(ox + 3 * cell, top, ox + 3 * cell, top + 16 * dy, "#444");
            svg.Line(ox, top + 13 * dy, ox + 4 * cell, top + 13 * dy, "#888");
            for (int col = 0; col < 4; col++)
                svg.Text(ox + (col + .5) * cell, 790, regions[col].Item3, 14, "middle");
        }
        const double barX = 410, barY = 854, barW = 650;
        for (int i = 0; i < 256; i++)
            svg.Rect(barX + barW * i / 256, barY, barW / 256 + 1, 21, Diverging((2.0 * i / 255 - 1) * limit, limit));
        foreach (double value in new[] { -limit, -limit / 2, 0, limit / 2, limit })
            svg.Text(barX + (value / limit + 1) * barW / 2, 896, value.ToString("0.#", CultureInfo.InvariantCulture), 13, "middle");
        svg.Text(735, 925, "Signed standardized residual: model − observation", 15, "middle");
        svg.Text(735, 966, "Three full-fit radial thirds; last column is withheld prediction · galaxy-equal medians", 14, "middle");
        svg.Save(path);
        Csv.Write(Path.ChangeExtension(path, "series.csv"), plotted);
        return limit;
    }

    static string Diverging(double value, double limit)
    {
        (int R, int G, int B)[] colors = [(5, 48, 97), (33, 102, 172), (67, 147, 195), (146, 197, 222), (209, 229, 240), (247, 247, 247), (253, 219, 199), (244, 165, 130), (214, 96, 77), (178, 24, 43), (103, 0, 31)];
        double t = Math.Min(255, (int)(Math.Clamp((value / limit + 1) / 2, 0, 1) * 256)) / 255.0 * 10;
        int i = Math.Min(9, (int)t);
        double f = t - i;
        int C(int a, int b) => (int)Math.Round(a + (b - a) * f);
        return $"#{C(colors[i].R, colors[i + 1].R):x2}{C(colors[i].G, colors[i + 1].G):x2}{C(colors[i].B, colors[i + 1].B):x2}";
    }

    sealed class Svg
    {
        readonly StringBuilder text;
        static string F(double x) => x.ToString("0.###", CultureInfo.InvariantCulture);
        static string E(string s) => SecurityElement.Escape(s) ?? "";
        public Svg(int width, int height, string title) => text = new($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><title>{E(title)}</title><rect width=\"100%\" height=\"100%\" fill=\"white\"/><style>text{{font-family:Arial,sans-serif}}</style>");
        public void Text(double x, double y, string s, int size, string anchor = "start", string color = "#262626") => text.Append($"<text x=\"{F(x)}\" y=\"{F(y)}\" fill=\"{color}\" font-size=\"{size}\" text-anchor=\"{anchor}\">{E(s)}</text>");
        public void Line(double x1, double y1, double x2, double y2, string color, string dash = "") => text.Append($"<line x1=\"{F(x1)}\" y1=\"{F(y1)}\" x2=\"{F(x2)}\" y2=\"{F(y2)}\" stroke=\"{color}\" stroke-width=\"1\" stroke-dasharray=\"{dash}\"/>");
        public void Rect(double x, double y, double w, double h, string color, string stroke = "none") => text.Append($"<rect x=\"{F(x)}\" y=\"{F(y)}\" width=\"{F(w)}\" height=\"{F(h)}\" fill=\"{color}\" stroke=\"{stroke}\" stroke-width=\".6\"/>");
        public void Marker(double x, double y, string shape, string color, bool filled, double r)
        {
            string style = $"fill=\"{(filled ? color : "white")}\" stroke=\"{color}\" stroke-width=\"1.2\"";
            text.Append(shape == "circle" ? $"<circle cx=\"{F(x)}\" cy=\"{F(y)}\" r=\"{F(r)}\" {style}/>" : shape == "square"
                ? $"<rect x=\"{F(x - r)}\" y=\"{F(y - r)}\" width=\"{F(2 * r)}\" height=\"{F(2 * r)}\" {style}/>"
                : $"<path d=\"M{F(x)},{F(y - r - 1)}L{F(x - r - 1)},{F(y + r)}L{F(x + r + 1)},{F(y + r)}Z\" {style}/>");
        }
        public void Save(string path) => File.WriteAllText(path, text.ToString() + "</svg>");
    }
}


