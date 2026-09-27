using System.Globalization;
using System.Net;
using System.Text;

namespace DarkUniverse;

/// <summary>Renders effects and a readable comparison table directly from this run's results.</summary>
public static class CurrentResultsPlots
{
    public static Dictionary<string, object?> Run(string outputRoot)
    {
        string folder = Path.Combine(outputRoot, "current_results");
        Directory.CreateDirectory(folder);
        var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Current statistical results</title><style>body{font:16px system-ui;max-width:1100px;margin:35px auto;padding:0 20px}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}img{max-width:100%;height:auto}p{line-height:1.5}</style><h1>Statistics calculated in this run</h1><p>Negative effects favor the named primary model. Effects are equal-galaxy mean differences in squared standardized loss. Bars are marginal 95% pairs-bootstrap-t intervals; they do not invert the separate wild tests and are not simultaneous intervals. All analyses condition on frozen predictions.</p>");
        var records = new List<object>();
        foreach (string family in new[] { "envelope", "pressure", "charged" })
        {
            string relative = family == "envelope" ? "pointwise_comparison/results/primary_cluster_results.csv" : $"publication_statistics/{family}/continuous_comparisons.csv";
            string source = Path.Combine(outputRoot, relative);
            string hash = Data.FileSha256(source);
            var rows = Csv.Read(source);
            var document = new PlotDocument(family + ": current mean-loss comparisons", 1, 1)
            {
                Scope = "Calculated from this run's results. Marginal 95% intervals; no new fit or physical significance."
            };
            var panel = new PlotPanel("Primary minus comparator", "Mean galaxy loss difference", "")
            {
                YMin = -.6, YMax = rows.Count - .4, YTicks = [],
                Note = family == "envelope" ? "All comparisons use the same declared galaxy sample." :
                    "Profiled comparison uses feasible galaxies only and has no primary p-value."
            };
            double lower = Math.Min(0, rows.Min(row => Csv.Number(row, "bootstrap_t_ci_low")));
            double upper = Math.Max(0, rows.Max(row => Csv.Number(row, "bootstrap_t_ci_high")));
            double rawStep = (upper - lower) / 5;
            if (rawStep > 0 && double.IsFinite(rawStep))
            {
                double power = Math.Pow(10, Math.Floor(Math.Log10(rawStep)));
                double fraction = rawStep / power;
                double step = power * (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10);
                panel.XTicks = [];
                for (double tick = Math.Ceiling(lower / step) * step; tick <= upper; tick += step)
                    panel.XTicks[tick] = tick.ToString("G4", CultureInfo.InvariantCulture);
            }
            panel.Add("", [0, 0], [-.6, rows.Count - .4], "#aaaaaa");
            html.Append($"<h2>{family}</h2><p><a href=\"../{relative}\">Full calculated CSV</a></p><img src=\"{family}.svg\" alt=\"{family} mean-loss intervals\"><table><tr><th>Comparator</th><th>Galaxies</th><th>Effect [95% CI]</th><th>Holm p</th><th>Scope</th></tr>");
            for (int i = 0; i < rows.Count; i++)
            {
                var row = rows[i];
                bool feasibleOnly = row.TryGetValue("full_sample_continuous_defined", out string? defined) && !bool.Parse(defined);
                double effect = Csv.Number(row, "effect"), low = Csv.Number(row, "bootstrap_t_ci_low"), high = Csv.Number(row, "bootstrap_t_ci_high");
                if (!double.IsFinite(effect) || !double.IsFinite(low) || !double.IsFinite(high) || low > high)
                    throw new InvalidDataException("Invalid current-result plotting interval.");
                string name = row["comparator"];
                panel.YTicks[i] = Label(name);
                var mark = panel.Add(feasibleOnly ? "Feasible subset" : "Full sample", [effect], [(double)i], feasibleOnly ? "#ba742a" : "#24658c", true);
                mark.XLower = [low]; mark.XUpper = [high]; mark.Hollow = feasibleOnly;
                string p = feasibleOnly ? "Undefined" : Format(Csv.Number(row, family == "envelope" ? "p_holm4" : "p_wild_holm6"));
                string scope = feasibleOnly ? "Feasible subset; no primary test" : Csv.Flag(row, "monte_carlo_floor") ? "Monte Carlo resolution reached" : "Full sample";
                html.Append($"<tr><td>{WebUtility.HtmlEncode(name)}</td><td>{row["n_galaxies"]}</td><td>{Format(effect)} [{Format(low)}, {Format(high)}]</td><td>{p}</td><td>{scope}</td></tr>");
            }
            document.Panels.Add(panel);
            document.Save(Path.Combine(folder, family + ".svg"));
            html.Append("</table>");
            if (Data.FileSha256(source) != hash) throw new InvalidDataException("Statistical results changed during rendering.");
            records.Add(new { family, source = relative, source_sha256 = hash, svg = family + ".svg", comparisons = rows.Count });
        }
        html.Append("<p>Direction, effect size, fit failures and numerical resolution should be read together. A small nominal probability does not establish the cubic mechanism.</p></html>");
        File.WriteAllText(Path.Combine(folder, "index.html"), html.ToString());
        var report = new Dictionary<string, object?> { ["status"] = "pass", ["statistics_loaded_from_current_run"] = true,
            ["figures"] = records, ["scope"] = "Fresh rendering of recomputed frozen-prediction statistics; no new inference or fitted states." };
        Data.SaveJson(Path.Combine(folder, "report.json"), report);
        return report;
    }
    static string Format(double x) => x.ToString("G5", CultureInfo.InvariantCulture);
    static string Label(string model) => model switch
    {
        "inherited_baryons" or "baryon_only" => "Baryons", "profiled_baryons" => "Profiled", "derived_quartic" => "Quartic",
        "original_core" or "scalar_core" => "Core", "free_control" => "Free", "compact_plummer" => "Compact", "envelope" => "Envelope", "nfw" => "NFW", _ => model
    };
}
