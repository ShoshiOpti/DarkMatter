using System.Text;
using System.Net;

namespace DarkUniverse;

/// <summary>Separates a derived field-equation check from a fitted-profile bookkeeping audit.</summary>
public static class ReconstructionBridge
{
    public static Dictionary<string, object?> Run(string dataRoot, string outputRoot)
    {
        var equations = FieldReconstructionChecks.Run(Path.Combine(outputRoot, "field_reconstruction"));
        var profiles = ReconstructionAudit.Run(dataRoot, outputRoot);
        string folder = Path.Combine(outputRoot, "reconstruction");
        Directory.CreateDirectory(folder);
        var report = new Dictionary<string, object?>
        {
            ["status"] = "pass", ["equations"] = equations, ["profiles"] = profiles,
            ["scope"] = "Retained-order original two-field GPP checks and frozen-target identities. No new galaxy fit, halo time integration, occupation law or observational significance."
        };
        Data.SaveJson(Path.Combine(folder, "report.json"), report);
        var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Reconstructed dynamics and current plots</title><style>body{font:16px system-ui;max-width:1150px;margin:35px auto;padding:0 20px;line-height:1.6}img{width:100%;height:auto}section{border-top:1px solid #ddd;margin-top:30px}</style><h1>Reconstructed dynamics and current plots</h1><p>A normalized internal decomposition preserves total density and the leading Poisson source. Its internal stress can change subsequent density evolution. Projection onto one component instead removes density. These are distinct operations.</p><p>The field fixture below is dimensionless and periodic; it checks instantaneous equations, not a galaxy simulation. Observational panels are algebraic audits of already fitted profiles; they do not estimate internal phases or a new halo population.</p><p><a href=\"report.json\">Combined report</a></p>");
        foreach (string child in new[] { "field_reconstruction", "reconstruction_audit" })
        {
            string directory = Path.Combine(outputRoot, child);
            html.Append($"<p><a href=\"../{child}/report.json\">{WebUtility.HtmlEncode(child)} report</a></p>");
            foreach (string file in Directory.EnumerateFiles(directory, "*.svg").Order())
            {
                string relative = "../" + child + "/" + Path.GetFileName(file);
                html.Append($"<section><h2>{WebUtility.HtmlEncode(Path.GetFileNameWithoutExtension(file).Replace('_', ' '))}</h2><a href=\"{relative}\">SVG</a> &middot; <a href=\"{Path.ChangeExtension(relative, "series.csv")}\">Plotted coordinates</a><img src=\"{relative}\" alt=\"Scientific reconstruction diagnostic\"></section>");
            }
        }
        html.Append("<p>The ordinary publication run separately recomputes all three statistical families and validates the 25 frozen figures against the existing result contract. Read statistical_impact.json in that run for the unchanged-result audit. No new p-value is assigned to an algebraic identity.</p></html>");
        File.WriteAllText(Path.Combine(folder, "index.html"), html.ToString());
        return report;
    }

    public static Dictionary<string, object?> WriteStatisticalImpact(string outputRoot)
    {
        string compatibilityPath = Path.Combine(outputRoot, "publication", "result_compatibility.json");
        var compatibility = Data.Json(compatibilityPath);
        if (compatibility.GetProperty("status").GetString() != "pass")
            throw new InvalidDataException("Cannot report unchanged inference without passing the figure/result contract.");
        var report = new Dictionary<string, object?>
        {
            ["status"] = "pass", ["new_observational_fit"] = false,
            ["change_in_predictions_from_coordinate_reconstruction"] = "None for the same total-density state at leading Poisson order; not a claim of equivalence for independently evolved states.",
            ["frozen_result_contract"] = compatibility,
            ["statistics_validation_sha256"] = Data.FileSha256(Path.Combine(outputRoot, "publication_statistics", "validation.json")),
            ["model_summaries"] = new[] { "pressure", "charged" }.ToDictionary(f => f,
                f => Csv.Read(Path.Combine(outputRoot, "publication_statistics", f, "model_summary.csv"))),
            ["inference"] = "Three original families recomputed using their declared 99999-draw protocols. No new hypothesis test for change of field coordinates. Existing failures, retrospective conditioning and multiple-comparison families retained.",
            ["interpretation_change"] = "Full empirical minus original-core force also includes the fitted core-amplitude change. The discarded orthogonal density of the inverse lift corresponds to its envelope and requires comparison with the rescaled projected core. Restoring internal variables is not a derived envelope occupation law."
        };
        Data.SaveJson(Path.Combine(outputRoot, "reconstruction", "statistical_impact.json"), report);
        return report;
    }
}
