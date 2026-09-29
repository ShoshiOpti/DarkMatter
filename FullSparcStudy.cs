using System.Globalization;
using System.Net;
using System.Text;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public sealed record FullSparcPrediction(string Galaxy, string Phase, string Model, int Index, double Radius,
    double Observed, double Sigma, double Predicted, bool Outer, bool Quality, bool NegativeBaryonForce)
{
    public double Q => Math.Pow((Predicted - Observed) / Sigma, 2);
    public Result Row() => new()
    {
        ["galaxy"] = Galaxy, ["phase"] = Phase, ["model"] = Model, ["index"] = Index, ["r_kpc"] = Radius,
        ["observed"] = Observed, ["sigma"] = Sigma, ["predicted"] = Predicted, ["outer"] = Outer,
        ["quality_sample"] = Quality, ["negative_catalogue_baryon_force"] = NegativeBaryonForce
    };
}

/// <summary>29 September all-SPARC native refitting, inference and plotting, separate from historical publication replay.</summary>
public static class FullSparcStudy
{
    public const string Scope = "Native all-175 conditional halo-template refit; frozen theory density inputs; no shared-microphysics or baryon-forced field solution.";
    public static Result Run(string dataRoot, string outputRoot, bool publishLanding = true)
    {
        string input = Path.Combine(dataRoot, "sparc_full"), output = Path.Combine(outputRoot, "sparc_full");
        Directory.CreateDirectory(output);
        var data = LoadGalaxies(input); var profiles = LoadProfiles(input);
        var profileCheck = CheckProfiles(input, profiles);
        var predictions = new List<FullSparcPrediction>(); var states = new List<Result>();
        int index = 0;
        foreach (var d in data)
        {
            foreach (string phase in new[] { "inner", "all" })
            {
                var models = new Dictionary<string, double[]> { ["baryons"] = d.BaryonSquared.Select(b => Math.Sqrt(Math.Max(b, 0))).ToArray() };
                foreach (var (model, profile) in profiles)
                {
                    var fit = FullSparcFitter.Fit(d, profile, phase == "inner"); models[model] = fit.Predicted;
                    states.Add(new()
                    {
                        ["galaxy"] = d.Name, ["phase"] = phase, ["model"] = model,
                        ["L_kpc"] = fit.Length, ["M_total_Msun"] = fit.TotalMass, ["amplitude_GM_over_L"] = fit.Amplitude,
                        ["fit_chi2"] = fit.Chi2, ["length_bound"] = fit.LengthBound, ["weak_field_bound"] = fit.PotentialBound,
                        ["mass_outside_last_observation_fraction"] = fit.ExteriorMassFraction,
                        ["central_potential_over_c2"] = fit.CentralPotential
                    });
                }
                foreach (var (model, velocity) in models)
                    for (int j = 0; j < velocity.Length; j++) predictions.Add(new(d.Name, phase, model, j, d.R[j], d.Y[j],
                        d.Error[j], velocity[j], j >= d.InnerCount, d.Quality, d.BaryonSquared[j] < 0));
            }
            if (++index % 25 == 0) Console.WriteLine($"Native SPARC fits: {index}/175 galaxies");
        }
        Csv.Write(Path.Combine(output, "predictions.csv"), predictions.Select(p => p.Row()));
        Csv.Write(Path.Combine(output, "fit_states.csv"), states);
        var fitCheck = CheckPredictions(input, predictions);
        var score = Scores(data, predictions);
        Csv.Write(Path.Combine(output, "galaxy_scores.csv"), score);
        var lookup = score.ToDictionary(s => ((string)s["galaxy"]!, (string)s["phase"]!, (string)s["model"]!));
        var summaries = new Dictionary<string, object>(); var inference = new Dictionary<string, object>();
        var nuisance = new Dictionary<string, object>();
        var random = Data.Json(Path.Combine(input, "random_contract.json"));
        foreach (bool quality in new[] { false, true })
        {
            string label = quality ? "quality131" : "all175";
            var sample = data.Where(d => !quality || d.Quality).ToArray();
            var counts = FullSparcStatistics.BootstrapCounts(random, sample.Select(d => d.Group).Distinct().Count());
            Console.WriteLine($"Verified full PCG64 multinomial stream; calculating {label} inference.");
            var summary = new List<Result>();
            foreach (string phase in new[] { "inner", "all" })
                foreach (string model in new[] { "baryons", "scalar", "new_analytic", "old_exact", "mean_response" })
                {
                    var s = sample.Select(d => lookup[(d.Name, phase, model)]).ToArray();
                    summary.Add(new() { ["phase"] = phase, ["model"] = model, ["galaxies"] = sample.Length,
                        ["rows"] = s.Sum(r => (int)r["rows"]!), ["mean_loss"] = s.Average(r => (double)r["loss"]!),
                        ["mae"] = s.Average(r => (double)r["mae"]!), ["total_chi2"] = s.Sum(r => (double)r["chi2"]!) });
                }
            summaries[label] = summary;
            var tests = new List<Result>();
            foreach (var (a, b) in new[] { ("new_analytic", "baryons"), ("mean_response", "baryons"), ("new_analytic", "scalar"), ("mean_response", "scalar") })
            {
                var difference = sample.Select(d => (double)lookup[(d.Name, "inner", a)]["loss"]! - (double)lookup[(d.Name, "inner", b)]["loss"]!).ToArray();
                var result = FullSparcStatistics.Inference(difference, sample.Select(d => d.Group).ToArray(), counts);
                result["model"] = a; result["baseline"] = b; tests.Add(result);
            }
            FullSparcStatistics.AdjustHolm(tests); inference[label] = tests;
            nuisance[label] = NuisanceReplay(input, sample, counts);
        }
        var results = new Result { ["scope"] = Scope, ["summaries"] = summaries, ["inference"] = inference };
        Data.SaveJson(Path.Combine(output, "RESULTS.json"), results);
        Data.SaveJson(Path.Combine(output, "nuisance_sensitivity.json"), new { scope = "Native scoring/resampling of archived matched-nuisance predictions; the nuisance optimizer is not rerun.", samples = nuisance });
        var bounds = FullFitBounds(predictions);
        Data.SaveJson(Path.Combine(output, "full_fit_significance.json"), bounds);
        CheckStatistics(input, inference, nuisance, bounds);
        Render(output, data, predictions, score, states, bounds, inference, nuisance);
        File.Copy(Path.Combine(input, "FULL_FIT_SIGNIFICANCE.md"), Path.Combine(output, "FULL_FIT_SIGNIFICANCE.md"));
        var report = new Result
        {
            ["status"] = "pass", ["scope"] = Scope, ["galaxies"] = data.Length, ["rows"] = data.Sum(d => d.R.Length),
            ["native_fit_count"] = states.Count, ["profiles"] = profileCheck, ["prediction_regression"] = fitCheck,
            ["random_stream"] = "Complete native 99999 x 148 and 99999 x 111 count streams match NumPy SHA-256 contracts.",
            ["statistics_regression"] = "Primary bootstrap/sign counts and bounds checked against independent Python references.",
            ["field_profiles"] = "Versioned numerical theory inputs; original solver and mesh diagnostics supplied for provenance, not invoked by .NET.",
            ["nuisance_profiles"] = "Archived nuisance-fitted predictions rescored and resampled natively; no fresh nuisance optimization.",
            ["figures"] = 32, ["atlas_panels"] = data.Length, ["outer_mass"] = "Integrated to infinity; never truncated at observed Rmax.",
            ["full_fit_bounds"] = bounds
        };
        Data.SaveJson(Path.Combine(output, "verification.json"), report);
        if (publishLanding)
            File.WriteAllText(Path.Combine(outputRoot, "index.html"), "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"0;url=sparc_full/index.html\"><title>Full SPARC</title><a href=\"sparc_full/index.html\">Open full SPARC results</a></html>");
        return report;
    }

    public static FullSparcGalaxy[] LoadGalaxies(string input)
    {
        var meta = new Dictionary<string, (bool quality, string group)>();
        foreach (string line in File.ReadLines(Path.Combine(input, "inputs", "SPARC_Lelli2016c.mrt")))
        {
            string[] a = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (a.Length != 19 || !double.TryParse(a[5], CultureInfo.InvariantCulture, out double inclination) || !int.TryParse(a[17], out int q)) continue;
            meta.Add(a[0], (q < 3 && inclination >= 30, a[4] == "4" ? "distance:UrsaMajor" : a[0]));
        }
        var raw = meta.Keys.ToDictionary(n => n, _ => new List<double[]>());
        foreach (string line in File.ReadLines(Path.Combine(input, "inputs", "MassModels_Lelli2016c.mrt")))
        {
            string[] a = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (a.Length == 10 && raw.TryGetValue(a[0], out var list)) list.Add(a.Skip(1).Select(x => double.Parse(x, CultureInfo.InvariantCulture)).ToArray());
        }
        var result = raw.OrderBy(p => p.Key, StringComparer.Ordinal).Select(pair =>
        {
            var a = pair.Value.OrderBy(r => r[1]).ToArray(); var m = meta[pair.Key];
            return new FullSparcGalaxy(pair.Key, a.Select(r => r[1]).ToArray(), a.Select(r => r[2]).ToArray(), a.Select(r => r[3]).ToArray(),
                a.Select(r => r[4] * Math.Abs(r[4]) + .5 * r[5] * Math.Abs(r[5]) + .7 * r[6] * Math.Abs(r[6])).ToArray(), m.quality && a.Length >= 8, m.group);
        }).ToArray();
        if (result.Length != 175 || result.Sum(d => d.R.Length) != 3391 || result.Count(d => d.Quality) != 131 ||
            result.Any(d => d.Error.Any(e => !double.IsFinite(e) || e <= 0) || d.R.Any(r => r <= 0)))
            throw new InvalidDataException("Unexpected full SPARC sample.");
        return result;
    }

    public static Dictionary<string, IFullSparcProfile> LoadProfiles(string input)
    {
        var rows = Csv.Read(Path.Combine(input, "field_profiles.csv")); var x = rows.Column("x");
        var diag = Data.Json(Path.Combine(input, "profile_diagnostics.json"));
        var core = new FullSparcDensity(x, rows.Column("core_radial_density"), Math.Sqrt(2 * Math.Abs(diag.GetProperty("fine").GetProperty("mu").GetDouble())));
        var mean = new FullSparcDensity(x, rows.Column("response_radial_density"), diag.GetProperty("response_tail_decay_rate").GetDouble());
        var p = Data.Json(Path.Combine(input, "inputs", "compression.json")).GetProperty("parameters").GetProperty("envelope_mode_0");
        double sc = p.GetProperty("core_dilation").GetDouble(), fraction = p.GetProperty("envelope_fraction").GetDouble(),
            a = p.GetProperty("a_over_L").GetDouble(), b = p.GetProperty("b_over_L").GetDouble();
        return new() { ["scalar"] = core, ["new_analytic"] = new FullSparcEnvelope(core, sc, fraction, a, b),
            ["old_exact"] = new FullSparcEnvelope(core, sc, fraction, a, b, true), ["mean_response"] = mean };
    }

    static Result CheckProfiles(string input, Dictionary<string, IFullSparcProfile> profiles)
    {
        var expected = Data.Json(Path.Combine(input, "profile_contract.json")); var x = expected.GetProperty("radii").EnumerateArray().Select(v => v.GetDouble()).ToArray();
        double maxMass = 0, maxDepth = 0;
        foreach (var (name, p) in profiles)
        {
            var gold = expected.GetProperty("profiles").GetProperty(name);
            var y = gold.GetProperty("mass").EnumerateArray().Select(v => v.GetDouble()).ToArray();
            maxDepth = Math.Max(maxDepth, Math.Abs(p.Depth - gold.GetProperty("depth").GetDouble()));
            for (int i = 0; i < x.Length; i++) maxMass = Math.Max(maxMass, Math.Abs(p.Mass(x[i]) - y[i]));
        }
        if (maxMass > 2e-12 || maxDepth > 2e-12) throw new InvalidDataException($"Native profile contract failed: mass {maxMass}, depth {maxDepth}.");
        return new() { ["max_CDF_error"] = maxMass, ["max_depth_error"] = maxDepth, ["tolerance"] = 2e-12 };
    }

    static Result CheckPredictions(string input, List<FullSparcPrediction> predictions)
    {
        var reference = Csv.Read(Path.Combine(input, "reference", "predictions.csv")).ToDictionary(r => (r["galaxy"], r["phase"], r["model"], int.Parse(r["index"])));
        double error = 0; string worst = "";
        foreach (var p in predictions)
        {
            var r = reference[(p.Galaxy, p.Phase, p.Model, p.Index)];
            if (p.Observed != Csv.Number(r, "observed") || p.Sigma != Csv.Number(r, "sigma") || p.Radius != Csv.Number(r, "r_kpc")) throw new InvalidDataException("Observation identity mismatch.");
            double difference = Math.Abs(p.Predicted - Csv.Number(r, "predicted"));
            if (difference > error) { error = difference; worst = $"{p.Galaxy}/{p.Phase}/{p.Model}/{p.Index}"; }
        }
        if (error > .001) throw new InvalidDataException($"Native prediction differs by {error:G8} km/s at {worst} (limit .001).");
        return new() { ["rows_compared"] = predictions.Count, ["maximum_velocity_error_kms"] = error, ["worst"] = worst, ["tolerance_kms"] = .001 };
    }

    static List<Result> Scores(FullSparcGalaxy[] data, List<FullSparcPrediction> predictions)
    {
        var group = data.ToDictionary(d => d.Name, d => d.Group);
        return predictions.Where(p => p.Phase == "all" || p.Outer).GroupBy(p => (p.Galaxy, p.Phase, p.Model)).Select(g => new Result
        {
            ["galaxy"] = g.Key.Galaxy, ["phase"] = g.Key.Phase, ["model"] = g.Key.Model, ["quality_sample"] = g.First().Quality,
            ["group"] = group[g.Key.Galaxy], ["rows"] = g.Count(), ["loss"] = g.Average(p => p.Q),
            ["mae"] = g.Average(p => Math.Abs(p.Predicted - p.Observed)), ["chi2"] = g.Sum(p => p.Q)
        }).ToList();
    }

    static Result NuisanceReplay(string input, FullSparcGalaxy[] sample, byte[] counts)
    {
        var saved = Csv.Read(Path.Combine(input, "reference", "full_nuisance_predictions.csv")).ToDictionary(r => (r["galaxy"], r["model"], int.Parse(r["index"])));
        double Loss(FullSparcGalaxy d, string model) => Enumerable.Range(d.InnerCount, d.R.Length - d.InnerCount)
            .Average(i => Math.Pow((Csv.Number(saved[(d.Name, model, i)], "predicted") - d.Y[i]) / d.Error[i], 2));
        var baryons = sample.Select(d => Loss(d, "baryons")).ToArray(); var halo = sample.Select(d => Loss(d, "new_analytic")).ToArray();
        var result = FullSparcStatistics.Inference(halo.Zip(baryons, (a, b) => a - b).ToArray(), sample.Select(d => d.Group).ToArray(), counts);
        result["baryon_mean_loss"] = baryons.Average(); result["new_mean_loss"] = halo.Average();
        result["fit_source"] = "Archived matched-nuisance fits; observations scored natively, not a new nuisance fit.";
        return result;
    }

    static Result FullFitBounds(List<FullSparcPrediction> predictions)
    {
        var all = predictions.Where(p => p.Phase == "all").ToArray();
        Result Calculate(bool physical, bool quality)
        {
            var selected = all.Where(p => (!physical || !p.NegativeBaryonForce) && (!quality || p.Quality)).ToArray();
            var b = selected.Where(p => p.Model == "baryons").ToArray();
            var result = FullSparcStatistics.Bounds(b.Length, b.Sum(p => p.Q), selected.Where(p => p.Model == "new_analytic").Sum(p => p.Q));
            result["galaxies"] = b.Select(p => p.Galaxy).Distinct().Count(); return result;
        }
        return new() { ["primary_physical_rows"] = Calculate(true, false), ["all3391_zero_speed_convention"] = Calculate(false, false),
            ["quality_physical_rows"] = Calculate(true, true),
            ["null"] = "Fixed grey mean, correct marginal sigma; no uncertainty in distance, inclination, stellar mass or baryonic templates.",
            ["calibration"] = "T<=Q0; Gaussian union bound min(1,2N SF(sqrt(T/N))); variance-only Markov bound min(1,N/T). Arbitrary correlations permitted.",
            ["excluded_from_physical_statistic_only"] = "Two central UGC01281 negative-force rows. Original green fits unchanged; all 175 galaxies and exterior mass retained.",
            ["Z_convention"] = "Two-sided equivalent; conservative p upper bounds and Z lower bounds, not exact LRT tails or action-detection significance." };
    }

    static void CheckStatistics(string input, Dictionary<string, object> inference, Dictionary<string, object> nuisance, Result bounds)
    {
        var reference = Data.Json(Path.Combine(input, "reference", "RESULTS.json"));
        var nr = Data.Json(Path.Combine(input, "reference", "nuisance_sensitivity.json"));
        foreach (string sample in new[] { "all175", "quality131" })
        {
            var actual = (List<Result>)inference[sample]; var gold = reference.GetProperty("inference").GetProperty(sample).EnumerateArray().ToArray();
            for (int i = 0; i < actual.Count; i++)
            {
                foreach (string k in new[] { "bootstrap_exceedances", "group_wins", "group_losses", "group_ties" })
                    if (Convert.ToInt32(actual[i][k]) != gold[i].GetProperty(k).GetInt32()) throw new InvalidDataException($"Statistical count changed: {sample}/{i}/{k}.");
                if (Math.Abs((double)actual[i]["mean_difference"]! - gold[i].GetProperty("mean_difference").GetDouble()) > 1e-4)
                    throw new InvalidDataException("Mean-loss contract mismatch.");
                foreach (string k in new[] { "Holm_p", "Holm_Z", "Holm_group_sign_p", "Holm_group_sign_Z" })
                    if (Math.Abs((double)actual[i][k]! / gold[i].GetProperty(k).GetDouble() - 1) > 2e-10)
                        throw new InvalidDataException("p/Z contract mismatch: " + k);
            }
            var n = (Result)nuisance[sample];
            if ((double)n["two_sided_p"]! != nr.GetProperty(sample).GetProperty("two_sided_p").GetDouble()) throw new InvalidDataException("Nuisance p mismatch.");
        }
        var bf = Data.Json(Path.Combine(input, "reference", "full_fit_significance.json"));
        foreach (string sample in new[] { "primary_physical_rows", "all3391_zero_speed_convention", "quality_physical_rows" })
        {
            var actual = (Result)bounds[sample]!; var gold = bf.GetProperty(sample);
            foreach (string k in new[] { "gaussian_marginal_conservative_p", "variance_only_conservative_p", "gaussian_marginal_equivalent_two_sided_Z", "chi2_green" })
                if (Math.Abs((double)actual[k]! / gold.GetProperty(k).GetDouble() - 1) > 2e-6) throw new InvalidDataException("Full-fit bound mismatch: " + k);
        }
    }

    static void Render(string output, FullSparcGalaxy[] data, List<FullSparcPrediction> predictions, List<Result> scores,
        List<Result> states, Result bounds, Dictionary<string, object> inference, Dictionary<string, object> nuisance)
    {
        string figures = Path.Combine(output, "figures"); Directory.CreateDirectory(figures);
        var curves = predictions.GroupBy(p => (p.Galaxy, p.Phase, p.Model)).ToDictionary(g => g.Key, g => g.OrderBy(p => p.Index).Select(p => p.Predicted).ToArray());
        PlotPanel Panel(FullSparcGalaxy d)
        {
            var panel = new PlotPanel(d.Name + (d.Quality ? "" : " (outside quality subset)"), "Radius [kpc]", "Speed [km/s]")
                { XMin = 0, XMax = d.R[^1] * 1.025, YMin = Math.Min(0, d.Y.Zip(d.Error, (y, e) => y - e).Min()), Note = "Shaded/open points: outer evaluation; green uses all points." };
            panel.Shades.Add((d.R[d.InnerCount - 1], d.R[^1] * 1.025, "#eef0f2"));
            foreach (var (model, phase, label, color, dash) in new[] {
                ("baryons","inner","Baryons only","#737a82","5 3"), ("scalar","inner","Scalar (inner fit)","#326b91","2 3"),
                ("new_analytic","inner","Analytic envelope (inner fit)","#af3370",""), ("mean_response","inner","Mean response (inner fit)","#bd782a","6 3 1 3"),
                ("new_analytic","all","Analytic envelope (all fitted)","#288069","5 3") })
                panel.Add(label, d.R, curves[(d.Name, phase, model)], color).Dash = dash;
            foreach (bool outer in new[] { false, true })
            {
                int[] ix = Enumerable.Range(0, d.R.Length).Where(i => (i >= d.InnerCount) == outer).ToArray();
                panel.ErrorBars("", ix.Select(i => d.R[i]).ToArray(), ix.Select(i => d.Y[i]).ToArray(), ix.Select(i => d.Error[i]).ToArray()).Hollow = outer;
            }
            return panel;
        }
        void Sheet(string name, string title, IEnumerable<FullSparcGalaxy> galaxies)
        {
            var doc = new PlotDocument(title, 2, 3) { Scope = "Fixed catalogue calibration. Full exterior mass retained. Green: fitted data; pink: outer prediction." };
            doc.Panels.AddRange(galaxies.Select(Panel)); doc.Save(Path.Combine(figures, name + ".svg"));
        }
        string[] six = ["DDO154", "UGC11557", "F568-3", "NGC4559", "UGC06614", "NGC2955"];
        Sheet("figure3_updated", "Figure 3 revisited: complete-mass profiles", six.Select(n => data.Single(d => d.Name == n)));
        for (int i = 0; i < data.Length; i += 6) Sheet($"atlas_{i / 6 + 1:00}", $"Full SPARC | {i + 1}-{Math.Min(i + 6, data.Length)} of 175", data.Skip(i).Take(6));
        var scoreLookup = scores.ToDictionary(s => ((string)s["galaxy"]!, (string)s["phase"]!, (string)s["model"]!));
        var summary = new PlotDocument("Full-profile outer prediction and unobserved mass", 2, 1) { Scope = "Every galaxy retained; unobserved exterior mass is a model extrapolation." };
        var scatter = new PlotPanel("Outer prediction losses", "Baryon mean standardized loss", "Envelope mean standardized loss") { XLog = true, YLog = true };
        scatter.Add("Equal loss", [1e-5, 1e5], [1e-5, 1e5], "#999999").Dash = "3 3";
        foreach (bool quality in new[] { true, false })
        {
            var selected = data.Where(d => d.Quality == quality).ToArray();
            scatter.Add(quality ? "Original 131" : "Additional 44", selected.Select(d => (double)scoreLookup[(d.Name, "inner", "baryons")]["loss"]!).ToArray(),
                selected.Select(d => (double)scoreLookup[(d.Name, "inner", "new_analytic")]["loss"]!).ToArray(), quality ? "#af3370" : "#555555", true).Hollow = !quality;
        }
        summary.Panels.Add(scatter);
        var exterior = states.Where(s => (string)s["phase"]! == "inner" && (string)s["model"]! == "new_analytic").Select(s => (double)s["mass_outside_last_observation_fraction"]!).ToArray();
        var hist = new PlotPanel("Exterior mass retained", "Fraction outside last measured radius", "Galaxies") { XMin = 0, XMax = 1, YMin = 0 };
        hist.Add("All 175", Enumerable.Range(0, 20).Select(i => (i + .5) / 20).ToArray(), Enumerable.Range(0, 20).Select(i => (double)exterior.Count(v => Math.Min(19, (int)(v * 20)) == i)).ToArray(), "#288069").Bars = true;
        summary.Panels.Add(hist); summary.Save(Path.Combine(figures, "full_profile_summary.svg"));
        var b = (Result)bounds["primary_physical_rows"]!;
        string F(object? v) => Convert.ToDouble(v).ToString("G7", CultureInfo.InvariantCulture);
        var html = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Full SPARC reproduction</title><style>body{font:16px system-ui;max-width:1180px;margin:30px auto;padding:0 20px;line-height:1.5}img{width:100%;height:auto}td,th{padding:8px;text-align:left;border-bottom:1px solid #ddd}table{border-collapse:collapse}details{margin:15px 0}@media print{details{display:block}details>img{display:block}.atlas{break-before:page}}</style><h1>Full SPARC: native .NET results</h1>");
        html.Append("<p>All 175 galaxies / 3,391 observations. Halo mass and scale freshly fitted in .NET from versioned theory profiles. Full exterior mass retained. This is a conditional template test, not a universal-parameter solution.</p><p><a href=\"RESULTS.json\">Outer statistics</a> · <a href=\"full_fit_significance.json\">Green-fit significance</a> · <a href=\"predictions.csv\">Every prediction</a> · <a href=\"fit_states.csv\">Fitted states</a> · <a href=\"verification.json\">Verification</a> · <a href=\"nuisance_sensitivity.json\">Nuisance sensitivity</a> · <a href=\"FULL_FIT_SIGNIFICANCE.md\">Derivation</a></p>");
        html.Append($"<h2>Green fitted envelope versus fixed grey baryons</h2><p>Delta chi-square = {F(b["delta_chi2"])}. These are conservative upper bounds on p and lower bounds on equivalent two-sided Z, not exact likelihood-ratio tails.</p><table><tr><th>Assumptions</th><th>p upper bound</th><th>Z lower bound</th></tr><tr><td>Gaussian marginal errors</td><td>{F(b["gaussian_marginal_conservative_p"])}</td><td>{F(b["gaussian_marginal_equivalent_two_sided_Z"])}</td></tr><tr><td>Only correct marginal variances</td><td>{F(b["variance_only_conservative_p"])}</td><td>{F(b["variance_only_equivalent_two_sided_Z"])}</td></tr></table>");
        html.Append("<p>Both bounds permit arbitrary correlations but condition on the grey means and quoted error scales being correct. They omit baryonic/calibration uncertainty and do not establish the action. Two central UGC01281 negative-force rows are undefined for the physical null; the primary bound uses 3,389 rows, with all 175 galaxies retained. The all-row convention is separately reported. Green residual Q/N remains about 13.70.</p>");
        html.Append("<h2>Outer evaluation: distinct statistics</h2><table><tr><th>Sample</th><th>Comparison</th><th>Holm bootstrap p</th><th>Z</th><th>Holm sign p</th><th>Z</th></tr>");
        foreach (var (sample, value) in inference) foreach (var test in (List<Result>)value)
            html.Append($"<tr><td>{sample}</td><td>{test["model"]} / {test["baseline"]}</td><td>{F(test["Holm_p"])}</td><td>{F(test["Holm_Z"])}</td><td>{F(test["Holm_group_sign_p"])}</td><td>{F(test["Holm_group_sign_Z"])}</td></tr>");
        html.Append("</table><p>99999 group-pairs bootstrap draws; zero exceedances mark a Monte Carlo floor. Sign tests measure frequency of improvement, not detection of two-field physics. These values belong to inner-fit/outer-evaluation curves, not the green all-data curve.</p><h2>Matched nuisance sensitivity</h2>");
        foreach (var (sample, value) in nuisance) { var n = (Result)value; html.Append($"<p>{sample}: p={F(n["two_sided_p"])}, Z={F(n["equivalent_Z"])} (exploratory, unadjusted). Archived nuisance fits rescored and resampled natively.</p>"); }
        html.Append("<h2>Figure 3</h2><img src=\"figures/figure3_updated.svg\" alt=\"Six original Figure 3 galaxies\"><h2>Full-sample summary</h2><img src=\"figures/full_profile_summary.svg\" alt=\"Outer loss and exterior mass\"><h2>All 175 galaxies</h2><p>Each sheet is a vector SVG with a companion scientific-coordinate CSV. Open or print any sheet at full resolution.</p>");
        for (int i = 1; i <= 30; i++) html.Append($"<details class=\"atlas\"><summary>Atlas sheet {i:00} — <a href=\"figures/atlas_{i:00}.svg\">SVG</a> · <a href=\"figures/atlas_{i:00}.series.csv\">coordinates</a></summary><img loading=\"lazy\" src=\"figures/atlas_{i:00}.svg\" alt=\"Atlas sheet {i}\"></details>");
        html.Append("</html>"); File.WriteAllText(Path.Combine(output, "index.html"), html.ToString());
    }
}
