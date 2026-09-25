using Row = System.Collections.Generic.Dictionary<string, string>;
namespace DarkUniverse;

public static class ObservationPlots
{
    static readonly string[] Chosen = ["DDO154", "UGC11557", "F568-3", "NGC4559", "UGC06614", "NGC2955"];
    static readonly string[] Colors = ["#24658c", "#147e83", "#be5e28", "#777777", "#75528d"];
    static readonly string[] ModelNames = ["Fixed core", "Variable core", "NFW"];
    static string root = "", output = "";
    static readonly List<string> paths = [];
    static List<Row> Read(string relative) => Csv.Read(Path.Combine(root, "numerics", relative));
    static string T(Row r, string k) => Csv.Text(r, k);
    static double N(Row r, string k) => Csv.Number(r, k);
    static bool B(Row r, string k) => Csv.Flag(r, k);
    static List<Row> Select(IEnumerable<Row> rows, params (string, string)[] filters) => rows.Where(r => filters.All(f => T(r, f.Item1) == f.Item2)).ToList();
    static PlotDocument Doc(string title, int columns, params PlotPanel[] panels)
    {
        var d = new PlotDocument(title, columns, (panels.Length + columns - 1) / columns);
        d.Panels.AddRange(panels);
        return d;
    }
    static void Save(PlotDocument doc, string name)
    {
        var p = Path.Combine(output, name + ".svg");
        doc.Save(p);
        paths.Add(p);
    }
    static PlotPanel Panel(string title, string x, string y, bool xlog = false, bool ylog = false) => new(title, x, y) { XLog = xlog, YLog = ylog };
    static double Loss(Row r)
    {
        double v = N(r, "holdout_chi2_per_row");
        return B(r, "holdout_admissible") && double.IsFinite(v) && v >= 0 ? v : double.PositiveInfinity;
    }
    static void Identity(PlotPanel p, double lo, double hi) => p.Add("", [lo, hi], [lo, hi], "#999999").Dash = "3 3";
    static void Scatter(PlotPanel p, string label, IEnumerable<Row> rows, string x, string y, string color)
    {
        var a = rows.ToArray();
        p.Add(label, a.Column(x), a.Column(y), color, true);
    }
    static void Observed(PlotPanel p, List<Row> rows, string radius, string observed, string error, Func<Row, bool> train)
    {
        foreach (bool used in new[] { true, false })
        {
            var a = rows.Where(r => train(r) == used).ToList();
            if (a.Count == 0)
                continue;
            p.ErrorBars(used ? "Used in fit" : "Withheld outer", a.Column(radius), a.Column(observed), a.Column(error)).Hollow = !used;
        }
        var tr = rows.Where(train).ToList();
        var test = rows.Where(r => !train(r)).ToList();
        if (tr.Count > 0 && test.Count > 0)
            p.Shades.Add(((N(tr[^1], radius) + N(test[0], radius)) / 2, N(rows[^1], radius) * 1.025, "#f3eee6"));
    }
    public static string[] Run(string dataRoot, string outputRoot)
    {
        root = dataRoot;
        output = Path.Combine(outputRoot, "figures");
        paths.Clear();
        Fixed("results", false);
        Fixed("wide_bounds", true);
        Tiers();
        Guarded();
        Nuisance();
        Residuals();
        Paired();
        Investigation();
        Morphology(false);
        Morphology(true);
        Growth();
        Baryons();
        Family("results");
        Family("wide_results");
        Pilot();
        return paths.ToArray();
    }
    static void Fixed(string folder, bool legacy)
    {
        var fits = Read($"observations/{folder}/galaxy_fit_results.csv").Where(r => N(r, "ml_disk") == .5).ToList();
        var points = Read($"observations/{folder}/observed_points_and_predictions.csv");
        var curves = Read($"observations/{folder}/model_curves.csv");
        var names = Data.Json(Path.Combine(root, $"numerics/observations/{folder}/fit_summary.json")).GetProperty("displayed").EnumerateArray().Select(v => v.GetString()!).ToArray();
        var doc = new PlotDocument("Fixed catalogue inputs; additive mass profiles", 2, 3) { Scope = "Historical fixed-input procedure; descriptive fit agreement" };
        foreach (var name in names)
        {
            var p = Panel(name, "Radius [kpc]", "Circular speed [km/s]");
            p.YMin = 0;
            var a = Select(points, ("galaxy", name));
            var c = Select(curves, ("galaxy", name));
            foreach (var (key, label, col) in new[] { ("Vbar_kms", "Baryons", Colors[3]), ("Vnfw_kms", "NFW", Colors[2]), ("Vcore_kms", "Fixed isolated core", Colors[0]) })
                p.Add(label, c.Column("r_kpc"), c.Column(key), col);
            p.ErrorBars("SPARC", a.Column("r_kpc"), a.Column("Vobs_kms"), a.Column("e_Vobs_kms"));
            doc.Panels.Add(p);
        }
        string prefix = legacy ? "legacy/wide_bounds/" : "";
        Save(doc, prefix + "sparc_rotation_comparison");
        var left = Panel("All 153 selected galaxies", "Descriptive chi²/(N−2)", "Fraction of galaxies", true);
        left.Ecdf("Fixed core", fits.Column("core_reduced_chi2"), Colors[0]);
        left.Ecdf("NFW", fits.Column("nfw_reduced_chi2"), Colors[2]);
        var right = Panel("Same two scale parameters", "NFW chi²/(N−2)", "Core chi²/(N−2)", true, true);
        Scatter(right, "Galaxies", fits, "nfw_reduced_chi2", "core_reduced_chi2", Colors[0]);
        var vals = fits.Column("core_reduced_chi2").Concat(fits.Column("nfw_reduced_chi2")).ToArray();
        Identity(right, vals.Min() * .7, vals.Max() * 1.4);
        Save(Doc("Fixed-input templates; historical diagnostic", 2, left, right), prefix + "sparc_sample_diagnostics");
    }
    static void Tiers()
    {
        var rows = Read("fit_tiers/tier_predictions.csv");
        var dense = Read("fit_tiers/dense_curves.csv");
        foreach (var (tier, shortName) in new[] { ("fixed_inputs", "baseline"), ("expanded_family", "expanded") })
            foreach (string mode in new[] { "inner", "full" })
            {
                var doc = new PlotDocument($"{shortName}: {(mode == "inner" ? "withheld outer prediction" : "all-radius fitted agreement")}", 2, 3);
                foreach (string name in Chosen)
                {
                    var refRows = Select(rows, ("galaxy", name), ("tier", "fixed_inputs"), ("model", "core"), ("fit_mode", "inner")).OrderBy(r => N(r, "radius_kpc")).ToList();
                    var p = Panel(name, "Catalogue radius [kpc]", "Catalogue speed [km/s]");
                    p.XMin = 0;
                    p.YMin = 0;
                    var candidates = rows.Where(r => T(r, "galaxy") == name).Column("prediction_kms").Concat(dense.Where(r => T(r, "galaxy") == name).Column("prediction_kms")).Concat(refRows.Select(r => N(r, "observed_kms") + N(r, "error_kms"))).Concat(refRows.Select(r => Math.Sqrt(Math.Max(0, N(r, "catalogue_baryon_v2"))))).Where(double.IsFinite);
                    p.YMax = candidates.Max() * 1.08;
                    p.XMax = N(refRows[^1], "radius_kpc") * 1.025;
                    if (tier == "expanded_family")
                    {
                        var b = Select(dense, ("galaxy", name), ("tier", "matched_nuisance"), ("model", "core"), ("fit_mode", mode));
                        p.Add("Core: nuisance only", b.Column("radius_kpc"), b.Column("prediction_kms"), "#528ab0").Dash = "5 2 1 2";
                    }
                    foreach (string model in new[] { "core", "nfw" })
                    {
                        var c = Select(dense, ("galaxy", name), ("tier", tier), ("model", model), ("fit_mode", mode));
                        var s = p.Add(model == "core" ? "Core" : "NFW", c.Column("radius_kpc"), c.Column("prediction_kms"), model == "core" ? Colors[tier == "fixed_inputs" ? 0 : 1] : Colors[2]);
                        if (model == "nfw")
                            s.Dash = "6 3";
                        if (model == "core")
                            p.Add("Catalogue baryons", c.Column("radius_kpc"), c.Select(r => Math.Sqrt(N(r, "catalogue_baryon_v2"))).ToArray(), Colors[3]).Dash = "2 3";
                    }
                    Observed(p, refRows, "radius_kpc", "observed_kms", "error_kms", r => mode == "full" || !B(r, "is_outer"));
                    doc.Panels.Add(p);
                }
                Save(doc, $"tier_{shortName}_{(mode == "inner" ? "outer" : "full")}");
            }
        var panels = new List<PlotPanel>();
        foreach (string mode in new[] { "full", "inner" })
        {
            var p = Panel(mode == "full" ? "All radii fitted: 153 galaxies" : "Outer prediction: 131 galaxies", "Data loss chi²/N", "Fraction of galaxies", true);
            p.YMin = 0;
            p.YMax = 1.015;
            foreach (var (tier, model, label, col) in new[] { ("fixed_inputs", "core", "Core: fixed inputs", Colors[0]), ("matched_nuisance", "core", "Core: fitted inputs", "#528ab0"), ("expanded_family", "core", "Core: inputs + shape", Colors[1]), ("fixed_inputs", "nfw", "NFW: fixed inputs", Colors[2]), ("expanded_family", "nfw", "NFW: fitted inputs", "#a65825") })
            {
                var a = Select(rows, ("tier", tier), ("model", model), ("fit_mode", mode)).Where(r => mode == "full" || B(r, "is_outer"));
                var loss = a.GroupBy(r => T(r, "galaxy")).Select(g => g.Average(r => Math.Pow((N(r, "prediction_kms") - N(r, "observed_kms")) / N(r, "error_kms"), 2))).ToArray();
                if (loss.Length != (mode == "full" ? 153 : 131) || loss.Any(v => !double.IsFinite(v)))
                    throw new InvalidDataException("Tier population mismatch");
                p.Ecdf(label, loss, col, true);
            }
            panels.Add(p);
        }
        Save(Doc("Fixed and fitted observational inputs", 2, panels.ToArray()), "tier_population_comparison");
    }
    static void Guarded()
    {
        var baseRows = Read("regime_baseline/results/guard_0.001_merged_fits.csv");
        var free = Read("population_family/guard_0.001/free_family_fits.csv");
        List<Row>[] groups = [Select(baseRows, ("model", "core")), free, Select(baseRows, ("model", "nfw"))];
        var p = new[] { Panel("Full-data fits", "Data chi² / point", "Fraction of 153 galaxies", true), Panel("Outermost radii withheld", "Test chi² / point", "Fraction of 131 galaxies", true), Panel("Variable-core eta boundaries", "Boundary", "Galaxies"), Panel("Potential-depth guard", "Unguarded |Phi(0)|/c²", "Guarded |Phi(0)|/c²", true, true) };
        for (int m = 0; m < 3; m++)
            for (int h = 0; h < 2; h++)
            {
                var a = groups[m].Where(r => B(r, "holdout") == (h == 1)).ToList();
                p[h].Ecdf(ModelNames[m], a.Select(r => h == 1 ? Loss(r) : N(r, "chi2_data_per_row")), Colors[m]);
            }
        p[2].XTicks = new()
        {
            [0] = "Near zero",
            [1] = "Interior",
            [2] = "Upper bound"
        };
        p[2].YMin = 0;
        p[2].YMax = 135;
        foreach (bool h in new[] { false, true })
        {
            var a = free.Where(r => B(r, "holdout") == h).ToList();
            p[2].Add(h ? "Inner training (131)" : "Full data (153)", [h ? .17 : -.17, 1 + (h ? .17 : -.17), 2 + (h ? .17 : -.17)], new[] { "lower", "none", "upper" }.Select(k => (double)a.Count(r => T(r, "eta_boundary") == k)).ToArray(), h ? "#80b9bb" : Colors[1]).Bars = true;
        }
        Scatter(p[3], "Variable core", free.Where(r => !B(r, "holdout")), "original_potential_depth_ratio", "potential_depth_ratio", Colors[1]);
        Identity(p[3], 1e-10, 100);
        p[3].XMin = 1e-10;
        p[3].XMax = 100;
        p[3].YMin = 1e-10;
        p[3].YMax = .005;
        p[3].Add("Guard 0.001", [1e-10, 100], [.001, .001], Colors[2]).Dash = "5 3";
        Save(Doc("Isolated profiles; variable core adds a local shape parameter", 2, p), "guarded_family_comparison");
    }
    static void Nuisance()
    {
        var rows = Read("observations/uncertainty/penalized_fits.csv");
        string[] configs = ["fixed", "stellar_0.10", "geometry", "all_0.05", "all_0.10", "all_0.20"];
        var a = Panel("Full-data fits: 153 galaxies", "Nuisance configuration", "Median data chi²/N");
        a.XTicks = configs.Select((s, i) => (s, i)).ToDictionary(v => (double)v.i, v => v.s);
        var b = Panel("Inner-only fits: outer predictions", "Outer chi²/Ntest", "Fraction of 131 galaxies", true);
        b.YMin = 0;
        b.YMax = 1.02;
        foreach (var (model, col, offset) in new[] { ("core", Colors[0], -.1), ("nfw", Colors[2], .1) })
        {
            a.Add(model, Enumerable.Range(0, 6).Select(i => i + offset).ToArray(), configs.Select(c => Data.Median(rows.Where(r => T(r, "model") == model && T(r, "configuration") == c && !B(r, "holdout")).Column("chi2_data_per_row"))).ToArray(), col);
            foreach (string cfg in new[] { "fixed", "all_0.10" })
                b.Ecdf(model + ": " + cfg, rows.Where(r => T(r, "model") == model && T(r, "configuration") == cfg && B(r, "holdout")).Select(Loss), col);
        }
        Save(Doc("Historical matched-prior sensitivity", 2, a, b), "nuisance_and_holdout");
    }
    static void Residuals()
    {
        var radial = Read("observations/diagnostics/radial_summary.csv");
        var rows = Select(Read("observations/diagnostics/galaxy_diagnostics.csv"), ("run", "ml_0.5"));
        var p = new[] { Panel("Median and IQR; R/Rmax bins", "Radial bin", "Galaxy-mean standardized residual"), Panel("Matched outer slopes", "Observed log slope", "Fitted log slope"), Panel("Coverage association", "Rmax/Rdisk", "log10(core chi²/NFW chi²)", true), Panel("NFW bound solutions", "Observed log slope", "Rmax/Rpeak", false, true) };
        p[0].XTicks = new()
        {
            [0] = "Inner",
            [1] = "Middle",
            [2] = "Outer"
        };
        foreach (var (model, col, off) in new[] { ("core", Colors[0], -.08), ("nfw", Colors[2], .08) })
        {
            var a = new[] { "inner", "middle", "outer" }.Select(bin => radial.Single(r => T(r, "run") == "ml_0.5" && T(r, "binning") == "r_over_rmax" && T(r, "model") == model && T(r, "radial_bin") == bin)).ToList();
            var s = p[0].Add(model, [off, 1 + off, 2 + off], a.Column("mean_z_median"), col, true);
            s.Lower = a.Column("mean_z_q25");
            s.Upper = a.Column("mean_z_q75");
            Scatter(p[1], model, rows, "observed_outer_slope", model + "_outer_slope", col);
        }
        var slopes = rows.Column("observed_outer_slope").Concat(rows.Column("core_outer_slope")).Concat(rows.Column("nfw_outer_slope")).Where(double.IsFinite).ToArray();
        Identity(p[1], slopes.Min() - .05, slopes.Max() + .05);
        foreach (bool flag in new[] { false, true })
        {
            Scatter(p[2], flag ? "Baryonic fraction >= 0.5" : "Baryonic fraction < 0.5", rows.Where(r => (N(r, "baryon_fraction_median") >= .5) == flag), "r_max_over_Rdisk", "log10_core_over_nfw_chi2", Colors[flag ? 2 : 0]);
            Scatter(p[3], flag ? "On bound" : "Interior", rows.Where(r => B(r, "nfw_bound_hit") == flag), "observed_outer_slope", "nfw_r_max_over_peak", Colors[flag ? 0 : 2]);
        }
        Save(Doc("Historical fixed-input residual structure", 2, p), "residual_physics_diagnostics");
    }
    static void Paired()
    {
        var stats = Read("significance/primary_summary.csv");
        var paired = Read("significance/paired_galaxy_losses.csv");
        var a = Panel("Conditional win counts: exact 95% intervals", "First model's win fraction", "Pair");
        a.XMin = 0;
        a.XMax = 1;
        a.YTicks = new()
        {
            [2] = "NFW/fixed",
            [1] = "NFW/variable",
            [0] = "Variable/fixed"
        };
        for (int i = 0; i < 3; i++)
        {
            var r = stats[i];
            double wins = N(r, i < 2 ? "right_wins" : "left_wins"), lo = N(r, "win_ci95_lower"), hi = N(r, "win_ci95_upper");
            if (i < 2)
                (lo, hi) = (1 - hi, 1 - lo);
            var s = a.Add($"{wins}/131; Holm Z={N(r, "z_holm_two_sided"):F2}", [wins / 131], [2 - i], Colors[i < 2 ? 2 : 1], true);
            s.XLower = [lo];
            s.XUpper = [hi];
        }
        var delta = paired.Where(r => T(r, "analysis") == "guard_0.001" && T(r, "split") == "holdout" && T(r, "left_model") == "variable_core" && T(r, "right_model") == "fixed_core").Column("delta_left_minus_right").Order().ToArray();
        var b = Panel("Ranked paired loss change", "Galaxy rank", "Variable loss − fixed loss");
        b.YSymLog = true;
        foreach (bool negative in new[] { true, false })
        {
            var ix = Enumerable.Range(0, delta.Length).Where(i => (delta[i] < 0) == negative).ToArray();
            b.Add(negative ? "Variable better" : "Fixed better", ix.Select(i => (double)i + 1).ToArray(), ix.Select(i => delta[i]).ToArray(), Colors[negative ? 1 : 2], true);
        }
        b.Note = $"Median {Data.Median(delta):F3}; galaxy-equal mean {delta.Average():F3}";
        Save(Doc("Retrospective galaxy-level prediction diagnostics", 2, a, b), "paired_prediction_diagnostics");
    }
    static void Investigation()
    {
        var baseline = Read("regime_baseline/results/guard_0.001_predictions.csv");
        var free = Read("population_family/guard_0.001/free_family_predictions.csv");
        List<Row>[] models = [Select(baseline, ("model", "core")), free, Select(baseline, ("model", "nfw"))];
        var raw = Sparc.Load(root);
        var summary = new List<(string Galaxy, int Model, bool Inner, double Q, double Mae)>();
        var names = models[0].Where(r => B(r, "holdout")).Select(r => T(r, "galaxy")).Distinct().Order().ToArray();
        foreach (string name in names)
            for (int m = 0; m < 3; m++)
                foreach (bool inner in new[] { false, true })
                {
                    var r = models[m].Where(r => T(r, "galaxy") == name && B(r, "holdout") == inner).OrderBy(r => N(r, "r_catalogue_kpc")).ToList();
                    int count = Sparc.HoldoutCount(r.Count);
                    var d = r.TakeLast(count).Select(r => (N(r, "Vpred_catalogue_kms") - N(r, "Vobs_kms")) / N(r, "Vobs_kms")).ToArray();
                    summary.Add((name, m, inner, d.Average(), d.Average(Math.Abs)));
                }
        PlotPanel Galaxy(string name, bool inner)
        {
            var p = Panel(name + (inner ? " | outer withheld" : " | all radii fitted"), "Catalogue radius [kpc]", "Catalogue speed [km/s]");
            p.XMin = 0;
            p.YMin = 0;
            var r = models[0].Where(r => T(r, "galaxy") == name && B(r, "holdout")).OrderBy(r => N(r, "r_catalogue_kpc")).ToList();
            p.YMax = models.SelectMany(v => v.Where(r => T(r, "galaxy") == name)).Column("Vpred_catalogue_kms").Concat(r.Select(r => N(r, "Vobs_kms") + N(r, "e_Vobs_kms"))).Max() * 1.075;
            p.XMax = N(r[^1], "r_catalogue_kpc") * 1.025;
            for (int m = 0; m < 3; m++)
            {
                var a = models[m].Where(v => T(v, "galaxy") == name && B(v, "holdout") == inner).OrderBy(v => N(v, "r_catalogue_kpc")).ToList();
                p.Add(ModelNames[m], a.Column("r_catalogue_kpc"), a.Column("Vpred_catalogue_kms"), Colors[m]).Dash = m == 0 ? "" : m == 1 ? "5 2 1 2" : "6 3";
            }
            p.Add("Catalogue baryons", raw.Points[name].Select(v => v.RadiusKpc).ToArray(), raw.Points[name].Select(v => Math.Sqrt(v.BaryonsSquared())).ToArray(), Colors[3]).Dash = "2 3";
            Observed(p, r, "r_catalogue_kpc", "Vobs_kms", "e_Vobs_kms", v => !inner || B(v, "training_row"));
            return p;
        }
        Save(Doc("Outer predictions from inner-only fits", 2, Chosen.Select(n => Galaxy(n, true)).ToArray()), "sparc_guarded_predictions");
        for (int plate = 0; plate < 2; plate++)
            Save(Doc("Fit and prediction: identical galaxy axes", 2, Chosen.Skip(plate * 3).Take(3).SelectMany(n => new[] { Galaxy(n, false), Galaxy(n, true) }).ToArray()), "sparc_fit_prediction_pair_" + (plate == 0 ? "a" : "b"));
        var panels = new[] { Panel("Fixed core: identical outer rows", "Inner-trained residual [%]", "Full-fit residual [%]"), Panel("Same 659 outer points", "Mean absolute fractional residual [%]", "Fraction of 131 galaxies"), Panel("Mass inside last radius", "Fitted field mass enclosed [%]", "Fraction of 131 galaxies"), Panel("Same-mass spherical ceiling", "Required / total fitted field mass", "Fraction of 131 galaxies", true) };
        foreach (var (lo, hi, label, col) in new[] { (1, 9, "Spiral (110)", Colors[1]), (10, 11, "Irregular/BCD (19)", Colors[2]), (0, 0, "S0 (2)", "#333333") })
        {
            var s = names.Where(n => raw.Catalogue[n].Type >= lo && raw.Catalogue[n].Type <= hi).ToArray();
            panels[0].Add(label, s.Select(n => 100 * summary.Single(v => v.Galaxy == n && v.Model == 0 && v.Inner).Q).ToArray(), s.Select(n => 100 * summary.Single(v => v.Galaxy == n && v.Model == 0 && !v.Inner).Q).ToArray(), col, true);
        }
        for (int m = 0; m < 3; m++)
            foreach (bool inner in new[] { true, false })
                panels[1].Ecdf(ModelNames[m] + (inner ? " inner" : " full"), summary.Where(v => v.Model == m && v.Inner == inner).Select(v => v.Mae * 100), Colors[m], true);
        var tails = Read("morphology/tail_results/tail_galaxy_diagnostics.csv");
        for (int m = 0; m < 2; m++)
        {
            var a = Select(tails, ("model", m == 0 ? "fixed_core" : "variable_core"));
            panels[2].Ecdf(ModelNames[m], a.Select(r => 100 * N(r, "mass_fraction_last")), Colors[m], true);
            panels[3].Ecdf(ModelNames[m], a.Select(r => N(r, "mass_fraction_last") * (1 + N(r, "last_required_additional_enclosed_over_fitted"))), Colors[m], true);
        }
        panels[3].XLog = false;
        panels[3].XSymLog = true;
        panels[3].XSymLogThreshold = 2;
        panels[3].XSymLogScale = 1.2;
        panels[3].XTicks = new()
        {
            [0] = "0",
            [1] = "1",
            [2] = "2",
            [10] = "10",
            [100] = "100",
            [1000] = "1000"
        };
        var doc = Doc("Conditional mass-budget diagnostics", 2, panels);
        doc.Scope = "Saved templates; spherical mass ceiling is not a new equilibrium";
        Save(doc, "outer_prediction_mass_budget");
    }
    static void Morphology(bool legacy)
    {
        var rows = Read("morphology/diagnostics/galaxy_model_diagnostics.csv").Where(r => N(r, "epsilon") == .001 && T(r, "analysis") == "outer_holdout").ToList();
        var p = new[] { Panel("Saved outer-block predictions", "Mean signed velocity residual [%]", "Fraction of galaxies"), Panel("Coverage and morphology", "Rmax/Rdisk", "Fixed-core residual [%]", true), Panel("Morphology covariate overlap", "Median observed speed [km/s]", "Rmax/Rdisk", true, true), Panel("Outer slope: median and IQR", "Model", "Predicted slope − observed slope") };
        string[] model = ["fixed_core", "variable_core", "nfw"];
        string[] morph = ["spiral_T1_9", "irregular_T10_11", "S0_T0"];
        for (int m = 0; m < 3; m++)
            for (int g = 0; g < 2; g++)
            {
                var a = Select(rows, ("model", model[m]), ("morphology", morph[g]));
                p[0].Ecdf(ModelNames[m] + (g == 0 ? " spiral" : " irregular"), a.Select(r => N(r, "fractional_residual_mean") * 100), Colors[m]);
                double med = Data.Median(a.Column("slope_residual"));
                var s = p[3].Add(m == 0 ? (g == 0 ? "Spiral" : "Irregular") : "", [m + (g - .5) * .18], [med], Colors[m], true);
                s.Lower = [Data.Quantile(a.Column("slope_residual"), .25)];
                s.Upper = [Data.Quantile(a.Column("slope_residual"), .75)];
            }
        for (int g = 0; g < 3; g++)
        {
            var a = Select(rows, ("model", "fixed_core"), ("morphology", morph[g]));
            p[1].Add(morph[g], a.Column("coverage"), a.Select(r => 100 * N(r, "fractional_residual_mean")).ToArray(), Colors[g == 0 ? 1 : g == 1 ? 2 : 3], true);
            if (g < 2)
                Scatter(p[2], morph[g], a, "median_speed", "coverage", Colors[g == 0 ? 1 : 2]);
        }
        p[3].XTicks = new()
        {
            [0] = "Fixed core",
            [1] = "Variable core",
            [2] = "NFW"
        };
        Save(Doc("Conditional morphology diagnostics", 2, p), (legacy ? "legacy/morphology/" : "") + "morphology_diagnostics");
    }
    static void Growth()
    {
        var rows = Read("outer_mechanism/growth/growth_galaxy_diagnostics.csv").Where(r => N(r, "guard") == .001 && N(r, "anchor_offset") == 0 && T(r, "model") != "variable_core_wide").ToList();
        var a = Select(rows, ("model", "fixed_core"));
        var b = Select(rows, ("model", "variable_core"));
        var p = Panel("Growth beyond training edge", "Enclosed mass growth G(Rmax)", "Fraction of galaxies", true);
        var q = Panel("Growth-only exterior diagnostic", "Outer loss", "Fraction of galaxies", true);
        p.Ecdf("Fixed core", a.Column("core_mass_growth"), Colors[0]);
        p.Ecdf("Variable core", b.Column("core_mass_growth"), Colors[1]);
        p.Ecdf("NFW", a.Column("nfw_mass_growth"), Colors[2]);
        q.Ecdf("Fixed: saved", a.Column("core_original_loss"), Colors[0]);
        q.Ecdf("Variable: saved", b.Column("core_original_loss"), Colors[1]);
        q.Ecdf("NFW: saved", a.Column("nfw_original_loss"), Colors[2]);
        q.Ecdf("Fixed: growth diagnostic", a.Column("core_with_nfw_growth_loss"), "#729ab4");
        q.Ecdf("Variable: growth diagnostic", b.Column("core_with_nfw_growth_loss"), "#79b6b8");
        q.XLog = false;
        q.XSymLog = true;
        q.XSymLogThreshold = 1;
        q.XSymLogScale = .7;
        q.XMin = 0;
        var doc = Doc("Conditional outer mass-growth diagnostic", 2, p, q);
        doc.Scope = "Algebraic counterfactual using saved profiles; not new equilibria";
        Save(doc, "outer_growth_mechanism");
    }
    static void Baryons()
    {
        var baseline = Read("baryons_only/baseline_predictions.csv");
        var halo = Read("fit_tiers/tier_predictions.csv");
        var saved = Read("baryons_only/galaxy_scores.csv");
        var models = new[] { ("bar_fixed", "fixed", ""), ("bar_profiled", "profiled", ""), ("A_core", "fixed_inputs", "core"), ("A_nfw", "fixed_inputs", "nfw"), ("B_core", "matched_nuisance", "core"), ("B_nfw", "matched_nuisance", "nfw"), ("C_core", "expanded_family", "core") };
        var losses = new Dictionary<(string, string), double[]>();
        int checks = 0;
        double maxError = 0;
        foreach (var (id, tier, model) in models)
            foreach (string mode in new[] { "full", "inner" })
            {
                var rows = model == "" ? Select(baseline, ("baseline", tier), ("fit_mode", mode)) : Select(halo, ("tier", tier), ("model", model), ("fit_mode", mode));
                var results = new List<double>();
                foreach (var group in rows.GroupBy(r => T(r, "galaxy")))
                {
                    var selected = group.Where(r => mode == "full" || B(r, "is_outer")).ToList();
                    bool failed = selected.Any(r => !double.IsFinite(N(r, "prediction_kms"))) || (tier == "profiled" && group.Any(r => B(r, "pipeline_failure")));
                    double loss = failed ? double.PositiveInfinity : selected.Average(r => Math.Pow((N(r, "prediction_kms") - N(r, "observed_kms")) / N(r, "error_kms"), 2));
                    var old = saved.Single(r => T(r, "model_id") == id && T(r, "fit_mode") == mode && T(r, "region") == (mode == "full" ? "all" : "outer") && T(r, "galaxy") == group.Key);
                    double expected = N(old, "loss");
                    if (double.IsFinite(loss))
                    {
                        double err = Math.Abs(loss - expected);
                        maxError = Math.Max(maxError, err);
                        if (err > 2e-10 + 3e-12 * Math.Abs(expected))
                            throw new InvalidDataException("Baryons loss mismatch");
                    }
                    else if (loss != expected)
                        throw new InvalidDataException("Baryons failure mismatch");
                    checks++;
                    results.Add(loss);
                }
                losses[(id, mode)] = results.ToArray();
            }
        var doc = new PlotDocument("Baryons-only reference and added mass profiles", 2, 2);
        foreach (bool matched in new[] { false, true })
            foreach (string mode in new[] { "full", "inner" })
            {
                var p = Panel((matched ? "Matched priors" : "Fixed inputs") + ": " + mode, "Per-galaxy loss chi²/N", "Fraction of all eligible galaxies", true);
                p.YMin = 0;
                p.YMax = 1.025;
                foreach (string id in matched ? new[] { "bar_profiled", "B_core", "C_core", "B_nfw" } : new[] { "bar_fixed", "A_core", "A_nfw" })
                {
                    var vals = losses[(id, mode)];
                    int failed = vals.Count(v => !double.IsFinite(v));
                    p.Ecdf(id + (failed > 0 ? $" ({failed} at infinity)" : ""), vals, id.StartsWith("bar") ? Colors[3] : id.Contains("nfw") ? Colors[2] : id == "C_core" ? Colors[1] : Colors[0], true);
                }
                doc.Panels.Add(p);
            }
        doc.Scope = "Undefined predictions retain their galaxies in ECDF denominators";
        Save(doc, "baryons_only_population");
        var gal = new PlotDocument("Independently fitted baryons-only versus added mass", 2, 3);
        var dense = Read("fit_tiers/dense_curves.csv");
        var galleryAxes = new Dictionary<string, double[]>();
        foreach (string name in Chosen)
        {
            var p = Panel(name, "Catalogue radius [kpc]", "Catalogue speed [km/s]");
            p.YMin = 0;
            p.XMin = 0;
            var r = Select(halo, ("tier", "fixed_inputs"), ("model", "core"), ("fit_mode", "inner"), ("galaxy", name));
            var candidates = halo.Where(v => T(v, "galaxy") == name).Column("prediction_kms").Concat(dense.Where(v => T(v, "galaxy") == name).Column("prediction_kms")).Concat(r.Select(v => N(v, "observed_kms") + N(v, "error_kms"))).Concat(r.Select(v => Math.Sqrt(N(v, "catalogue_baryon_v2")))).Where(double.IsFinite).ToArray();
            double oldMaximum = candidates.Max() * 1.08;
            double baryonMaximum = Select(baseline, ("baseline", "profiled"), ("fit_mode", "inner"), ("galaxy", name)).Column("prediction_kms").Where(double.IsFinite).DefaultIfEmpty(0).Max() * 1.08;
            p.YMax = Math.Max(oldMaximum, baryonMaximum);
            double lowest = r.Min(v => N(v, "observed_kms") - N(v, "error_kms"));
            p.YMin = lowest < 0 ? Math.Min(0, lowest - .025 * p.YMax.Value) : 0;
            p.XMax = N(r[^1], "radius_kpc") * 1.025;
            galleryAxes[name] = [p.XMin.Value, p.XMax.Value, p.YMin.Value, p.YMax.Value];
            foreach (var (id, tier, model) in models.Where(v => new[] { "bar_profiled", "B_core", "C_core", "B_nfw" }.Contains(v.Item1)))
            {
                var a = model == "" ? Select(baseline, ("baseline", tier), ("fit_mode", "inner"), ("galaxy", name)) : Select(halo, ("tier", tier), ("model", model), ("fit_mode", "inner"), ("galaxy", name));
                bool failure = model == "" && a.Any(r => B(r, "pipeline_failure"));
                p.Add(id, a.Column("radius_kpc"), failure ? a.Select(_ => double.NaN).ToArray() : a.Column("prediction_kms"), id.StartsWith("bar") ? Colors[3] : model == "nfw" ? Colors[2] : id == "C_core" ? Colors[1] : Colors[0]);
            }
            Observed(p, r, "radius_kpc", "observed_kms", "error_kms", r => !B(r, "is_outer"));
            gal.Panels.Add(p);
        }
        Save(gal, "baryons_only_outer_examples");
        Data.SaveJson(Path.Combine(output, "baryons_score_verification.json"), new
        {
            checks,
            maxError,
            passed = true,
            gallery_axes = galleryAxes
        });
    }
    static void Family(string folder)
    {
        var baseline = Select(Read("observations/uncertainty/penalized_fits.csv"), ("configuration", "all_0.10"));
        var free = Read($"population_family/{folder}/free_family_fits.csv");
        List<Row>[] groups = [Select(baseline, ("model", "core")), free, Select(baseline, ("model", "nfw"))];
        var p = new[] { Panel("Same observational priors", "Full data chi²/N", "Fraction of galaxies", true), Panel("Outer radii excluded from fitting", "Outer chi²/Ntest", "Fraction of galaxies", true), Panel("Physical scale product", "Eta (zero at 1e-10)", "L × amplitude [kpc km/s]", true, true), Panel("Paired galaxies", "NFW outer loss", "Family outer loss", true, true) };
        for (int m = 0; m < 3; m++)
            for (int h = 0; h < 2; h++)
            {
                var a = groups[m].Where(r => B(r, "holdout") == (h == 1)).ToList();
                p[h].Ecdf(ModelNames[m], a.Select(r => h == 1 ? Loss(r) : N(r, "chi2_data_per_row")), Colors[m]);
            }
        var full = free.Where(r => !B(r, "holdout")).ToList();
        p[2].Add("Free family", full.Select(r => Math.Max(1e-10, N(r, "eta"))).ToArray(), full.Select(r => N(r, "physical_length_kpc") * N(r, "physical_amplitude_kms")).ToArray(), Colors[1], true);
        var test = free.Where(r => B(r, "holdout")).ToList();
        var nfw = groups[2].Where(r => B(r, "holdout")).ToDictionary(r => T(r, "galaxy"));
        p[3].Add("Free family", test.Select(r => Loss(nfw[T(r, "galaxy")])).ToArray(), test.Select(Loss).ToArray(), Colors[1], true);
        Identity(p[3], .001, 1e5);
        Save(Doc("Historical unrestricted family diagnostic", 2, p), $"legacy/population_family/{folder}/family_population_comparison");
    }
    static void Pilot()
    {
        var rows = Read("coupled_pilot/predictions.csv");
        var doc = new PlotDocument("Conditional spherical-source two-galaxy pilot", 2, 2);
        var sharedAxes = new Dictionary<string, double[]>();
        foreach (string name in new[] { "F583-1", "UGC00731" })
        {
            var reference = Select(rows, ("galaxy", name), ("split", "inner"));
            double start = N(reference[0], "r_kpc"), end = N(reference[^1], "r_kpc"), pad = .03 * (end - start);
            var values = reference.Column("coupled_kms").Concat(reference.Column("nfw_kms")).Concat(reference.Select(r => Math.Sqrt(N(r, "baryon_v2")))).Concat(reference.Select(r => N(r, "observed_kms") - N(r, "error_kms"))).Concat(reference.Select(r => N(r, "observed_kms") + N(r, "error_kms"))).ToArray();
            sharedAxes[name] = [Math.Max(0, start - pad), end + pad, 0, values.Max() + .05 * (values.Max() - values.Min())];
        }
        foreach (string name in new[] { "F583-1", "UGC00731" })
            foreach (string mode in new[] { "inner", "full" })
            {
                var r = Select(rows, ("galaxy", name), ("split", mode));
                var p = Panel(name + " | " + mode, "Radius [kpc]", "Circular speed [km/s]");
                var limits = sharedAxes[name];
                p.XMin = limits[0];
                p.XMax = limits[1];
                p.YMin = limits[2];
                p.YMax = limits[3];
                p.Add("Coupled field + baryons", r.Column("r_kpc"), r.Column("coupled_kms"), Colors[0]);
                p.Add("NFW + baryons", r.Column("r_kpc"), r.Column("nfw_kms"), Colors[2]);
                p.Add("Baryons", r.Column("r_kpc"), r.Select(r => Math.Sqrt(N(r, "baryon_v2"))).ToArray(), Colors[3]);
                Observed(p, r, "r_kpc", "observed_kms", "error_kms", r => B(r, "training"));
                doc.Panels.Add(p);
            }
        doc.Scope = "Lines join saved observed-radius predictions; conditional pilot, not a population result";
        Save(doc, "coupled_galaxy_pilot");
        Data.SaveJson(Path.Combine(output, "pilot_shared_axes.json"), sharedAxes);
    }
}
