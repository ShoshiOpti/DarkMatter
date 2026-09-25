using System.Text;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

internal static class PointwiseReport
{
    public static string Write(string outputRoot)
    {
        string folder = Path.Combine(outputRoot, "pointwise_comparison");
        string results = Path.Combine(folder, "results");
        var primary = Csv.Read(Path.Combine(results, "primary_cluster_results.csv"));
        var scores = Csv.Read(Path.Combine(results, "score_sensitivity.csv"));
        var distributions = Csv.Read(Path.Combine(results, "residual_distribution_summary.csv"));
        var pooled = Csv.Read(Path.Combine(results, "pooled_point_sensitivity.csv")).Single(r => r["comparator"] == "nfw");
        var absolute = scores.Single(r => r["score"] == "absolute_kms" && r["comparator"] == "scalar_core");
        var nfw = primary.Single(r => r["comparator"] == "nfw");
        double Number(Row r, string key) => Csv.Number(r, key);
        double MeanAbsolute(string model) => Number(distributions.Single(r => r["model"] == model), "mean_abs_kms_equal_galaxy");
        double core = MeanAbsolute("scalar_core");
        double envelope = MeanAbsolute("extended");
        double reduction = -Number(absolute, "effect");
        // The saved contrast is envelope minus core; report the positive reduction.
        double low = -Number(absolute, "bootstrap_t_ci_high");
        double high = -Number(absolute, "bootstrap_t_ci_low");
        double nfwMae = MeanAbsolute("nfw");
        double nfwP = Number(nfw, "p_wild");
        double nfwAdjustedP = Number(nfw, "p_holm4");
        double nfwZ = Number(nfw, "Z_holm4");

        int floors = primary.Count(r => Csv.Flag(r, "monte_carlo_floor"));
        foreach (var row in primary)
            if (Csv.Flag(row, "monte_carlo_floor") != (Number(row, "wild_exceedances") == 0))
                throw new InvalidDataException("Monte Carlo floor flag is inconsistent with exceedance count.");
        var summary = new Dictionary<string, object?>
        {
            ["analysis"] = "Pointwise continuous-loss comparison, 25 September 2026",
            ["outer_measurements"] = 659,
            ["galaxies"] = 131,
            ["resampling_unit"] = "whole galaxy",
            ["replicates_per_stream"] = 99999,
            ["core_mae_kms"] = core,
            ["envelope_mae_kms"] = envelope,
            ["nfw_mae_kms"] = nfwMae,
            ["envelope_mae_reduction_vs_core_kms"] = reduction,
            ["reduction_marginal_bootstrap_t_ci95_kms"] = new[] { low, high },
            ["nfw_primary_wild_p"] = nfwP,
            ["nfw_primary_holm4_p"] = nfwAdjustedP,
            ["nfw_primary_two_sided_gaussian_equivalent_Z"] = nfwZ,
            ["monte_carlo_resolution_reached_comparisons"] = primary.Where(r => Csv.Flag(r, "monte_carlo_floor")).Select(r => r["comparator"]).ToArray(),
            ["resolution_note"] = "Zero exceedances in 99,999 replicates. Plus-one p and its Gaussian-equivalent Z are finite-resolution display conventions, not measured extreme tails or bounds on true significance.",
            ["interval_note"] = "Marginal pairs-bootstrap-t intervals do not invert the null-imposed wild-bootstrap test and are not simultaneous Holm confidence intervals.",
            ["historical_sign_tests"] = "Retained separately in output/significance; they test win frequency rather than mean error magnitude."
        };
        Data.SaveJson(Path.Combine(folder, "summary.json"), summary);
        var text = new StringBuilder("# Current pointwise statistical results — 25 September 2026\n\n");
        text.AppendLine("All 659 original outer measurements across 131 galaxies contribute. " +
            "Differences are averaged within galaxy, then galaxies receive equal weight. " +
            "Each bootstrap draw resamples or perturbs whole galaxies, preserving the grouping of their radial measurements. " +
            "Saved fits and catalogue observations are unchanged.\n");
        text.AppendLine($"Envelope versus core: mean absolute velocity error falls from **{core:F5} to {envelope:F5} km/s**, " +
            $"a reduction of **{reduction:F5} km/s**, " +
            $"with a marginal 95% whole-galaxy pairs-bootstrap-t interval of **[{low:F5}, {high:F5}] km/s**.\n");
        text.AppendLine($"NFW has smaller observed mean absolute error ({nfwMae:F5} km/s). " +
            $"The primary mean squared standardized-loss comparison gives **p = {nfwP:F5}, Holm4 p = {nfwAdjustedP:F5}, Z = {nfwZ:F5}**.\n");
        text.AppendLine("| Envelope minus comparator | Mean squared-loss difference | Marginal 95% bootstrap-t interval | Holm4 p | Two-sided Z | Monte Carlo status |");
        text.AppendLine("|---|---:|---|---:|---:|---|");
        foreach (var row in primary)
            text.AppendLine(FormatComparison(row));
        text.AppendLine($"\n{floors} comparisons have zero exceedances in 99,999 wild-bootstrap draws. " +
            "Their raw plus-one p is 0.00001; Holm4 p is 0.00004. " +
            "Their Z values are only the Gaussian equivalents of that reporting convention. " +
            "These are not precise extreme significances or bounds on an unknown tail probability. " +
            "Raw tail-count uncertainty is retained in primary_cluster_results.csv.\n");
        text.AppendLine($"The NFW mean-loss interval [{Number(nfw, "bootstrap_t_ci_low"):F3}, {Number(nfw, "bootstrap_t_ci_high"):F3}] and its p-value use different procedures; " +
            "the interval is not an inversion of the primary test. " +
            $"The pooled-point-weighted percentile interval is [{Number(pooled, "percentile_ci_low"):F3}, {Number(pooled, "percentile_ci_high"):F3}]. " +
            "Assumed radial-correlation sensitivities are recorded for rho = 0.25, 0.5, 0.75.\n");
        text.AppendLine("UGC01281 remains a failed independently profiled baryon fit. " +
            "That full-sample comparison has no finite effect inference or p-value; " +
            "its separate feasible-only sensitivity uses 130 galaxies / 654 measurements.\n");
        text.AppendLine("The analysis is retrospective and conditional on saved predictions. " +
            "Whole-galaxy resampling assumes independent galaxies and does not propagate a new fitting procedure. " +
            "Older sign-test results remain separate historical comparisons.");
        File.WriteAllText(Path.Combine(folder, "SUMMARY.md"), text.ToString());
        Console.WriteLine($"Current pointwise analysis: 659 measurements / 131 galaxies; core MAE {core:F2} -> envelope {envelope:F2} km/s; reduction {reduction:F2} [{low:F2}, {high:F2}].");
        Console.WriteLine($"NFW mean-loss p = {nfwP:F5}, Z = {nfwZ:F2}; {floors} primary comparisons flagged at the Monte Carlo resolution limit.");
        return "<h2>Current pointwise results — 25 September 2026</h2><p>659 outer measurements across 131 galaxies; whole-galaxy resampling.</p>"
            + $"<p>Core → envelope mean absolute error: <strong>{core:F2} → {envelope:F2} km/s</strong>. Reduction <strong>{reduction:F2} km/s</strong>, marginal 95% interval <strong>[{low:F2}, {high:F2}]</strong>. NFW primary p = <strong>{nfwP:F5}</strong>, Z = <strong>{nfwZ:F2}</strong>.</p>"
            + $"<p><strong>{floors} primary comparisons reach the Monte Carlo resolution limit.</strong> Their p/Z values are finite-resolution reporting conventions. Intervals and tests use different procedures; weighting and assumed radial correlations affect the NFW comparison.</p>"
            + "<p><a href=\"pointwise_comparison/SUMMARY.md\">Current results and interpretation</a> · <a href=\"pointwise_comparison/results/primary_cluster_results.csv\">Primary table</a> · <a href=\"pointwise_comparison/results/score_sensitivity.csv\">Score and covariance sensitivities</a></p>";
    }

    static string FormatComparison(Row row)
    {
        double Value(string column) => Csv.Number(row, column);
        string label = row["comparator"] switch
        {
            "scalar_core" => "Original core",
            "compact_plummer" => "Compact control",
            "nfw" => "NFW",
            "baryon_only" => "Core-calibrated baryons",
            _ => row["comparator"]
        };
        string status = Csv.Flag(row, "monte_carlo_floor")
            ? "Resolution reached; p/Z display convention"
            : $"{Value("wild_exceedances"):F0} exceedances";

        return $"| {label} | {Value("effect"):F6} | " +
            $"[{Value("bootstrap_t_ci_low"):F6}, {Value("bootstrap_t_ci_high"):F6}] | " +
            $"{Value("p_holm4"):G6} | {Value("Z_holm4"):F6} | {status} |";
    }
}
