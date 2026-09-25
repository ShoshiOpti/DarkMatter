using System.Globalization;
using MathNet.Numerics.Interpolation;
using System.Text.Json;

namespace DarkUniverse;

internal static class ControlledPlots
{
    const string Blue = "#2467A6", Orange = "#D07828", Green = "#278573", Red = "#A34766", Gray = "#626977";
    static readonly string[] Colors = [Gray, Blue, Orange, Green, Red, "#8055a1"];
    static string source = "", output = "";
    static readonly List<string> files = [];
    static readonly Dictionary<string, double> checks = [];
    static double N(Dictionary<string, string> row, string key) => Csv.Number(row, key);
    static double J(JsonElement row, string key) => row.GetProperty(key).GetDouble();
    static double[] A(IEnumerable<Dictionary<string, string>> rows, string key) => rows.Select(r => N(r, key)).ToArray();
    static double[] A(IEnumerable<JsonElement> rows, string key) => rows.Select(r => J(r, key)).ToArray();
    static List<Dictionary<string, string>> Read(string file) => Csv.Read(Path.Combine(source, "numerics", file));
    static JsonElement Json(string file) => JsonDocument.Parse(File.ReadAllText(Path.Combine(source, "numerics", file))).RootElement.Clone();
    static double[] Grid(double start, double stop, int count) => Enumerable.Range(0, count).Select(i => start + (stop - start) * i / (count - 1)).ToArray();
    static double[] Map(double[] x, Func<double, double> f) => x.Select(f).ToArray();
    static double[] Combine(double[] a, double[] b, Func<double, double, double> f) => a.Zip(b, f).ToArray();
    static void Check(string name, double error, double tolerance)
    {
        checks[name] = error;
        if (!double.IsFinite(error) || error > tolerance)
            throw new InvalidDataException($"{name}: {error:G17} > {tolerance:G17}");
    }
    static void Compare(string name, double[] actual, double[] expected, double tolerance = 1e-10)
    {
        if (actual.Length != expected.Length)
            throw new InvalidDataException(name + ": different lengths");
        Check(name, actual.Zip(expected, (a, b) => Math.Abs(a - b)).DefaultIfEmpty(0).Max(), tolerance);
    }
    static PlotPanel Panel(PlotDocument doc, string title, string x, string y, bool xlog = false, bool ylog = false)
    {
        var p = new PlotPanel(title, x, y) { XLog = xlog, YLog = ylog };
        doc.Panels.Add(p);
        return p;
    }
    static void Save(PlotDocument doc, string name)
    {
        doc.Scope = "Controlled calculation or exact formula; not an observed galaxy fit";
        var path = Path.Combine(output, name + ".svg");
        doc.Save(path);
        files.Add(path);
    }
    static void Line(PlotPanel p, string label, double x0, double x1, double y, string color = Gray)
        => p.Add(label, [x0, x1], [y, y], color).Dash = "4,4";

    public static string[] Run(string sourceRoot, string outputDir)
    {
        source = sourceRoot;
        output = outputDir;
        files.Clear();
        checks.Clear();
        Directory.CreateDirectory(output);
        Stress();
        Spectrum();
        Dynamics();
        Carrier();
        Geometry();
        Thermal();
        TwoField();
        Synthetic();
        Register();
        File.WriteAllText(Path.Combine(output, "controlled_validation.json"), JsonSerializer.Serialize(new
        {
            scope = "Controlled saved solutions and exact identities; no galaxy fitting, new field evolution, or significance inference.",
            plots = files.Select(Path.GetFileName),
            maximum_absolute_errors = checks
        }, new JsonSerializerOptions { WriteIndented = true }));
        return files.ToArray();
    }

    static void Stress()
    {
        var r = Read("stress_comparison/results/stress_comparison.csv");
        var x = A(r, "r_over_L");
        var d = new PlotDocument("Controlled equal-density excitation: equal masses, zero detuning", 2, 2);
        var a = Panel(d, "(a) Same initial gravitational curve", "r/L", "vc/vs");
        a.Add("Mixed reference, t = 0", x, A(r, "initial_vc_over_vs"), "#222222");
        var sample = r.Where((_, i) => i >= 12 && (i - 12) % 45 == 0).ToList();
        a.Add("Internal excitation, t = 0", A(sample, "r_over_L"), A(sample, "initial_vc_over_vs"), Orange, true);
        a = Panel(d, "(b) Initial radial accelerations", "r/L", "ar L/vs² (outward positive)");
        a.Add("Gravity on field and tracers", x, A(r, "gravity_acceleration_L_over_vs2"), Blue);
        a.Add("Extra support on field only", x, A(r, "internal_stress_acceleration_L_over_vs2"), Orange);
        a = Panel(d, "(c) Evolving-potential force diagnostic", "r/L", "Δ(vc²/vs²)/(t vs/L)²");
        a.Add("Exact t → 0 coefficient", x, A(r, "initial_delta_vc2_over_t2"), "#222222");
        foreach (var (key, time, color) in new[] { ("delta_vc2_t015", .15, Blue), ("delta_vc2_t030", .3, Orange) })
        {
            var y = Map(A(r, key), v => v / (time * time));
            a.Add($"Evolution, t = {time}", x, y, color);
            Compare("stress_" + key, y, A(r, key.Replace("delta_vc2_", "delta_vc2_over_t2_")));
        }
        a = Panel(d, "(d) Fixed-profile confinement requirement", "r/L", "ΔMconf(<r)/M0");
        a.Add("Required confinement mass", x, A(r, "required_confinement_mass_over_core_mass"), Red);
        int start = -1;
        for (int i = 0; i <= r.Count; i++)
        {
            bool negative = i < r.Count && N(r[i], "required_confinement_density") < 0;
            if (negative && start < 0)
                start = i;
            if (!negative && start >= 0)
            {
                a.Shades.Add((x[start], x[i - 1], "#f6eaee"));
                start = -1;
            }
        }
        a.Note = "Shading: negative source density required";
        foreach (var p in d.Panels)
        {
            p.XMin = 0;
            p.XMax = 8;
        }
        Save(d, "gr_internal_stress_comparison");
    }

    static void Spectrum()
    {
        var s = Json("core/results/core_diagnostics.json").GetProperty("spectrum");
        var runs = s.GetProperty("runs").EnumerateArray().ToArray();
        var finest = runs[^1];
        var d = new PlotDocument("Isolated mixed core: equal masses, zero detuning", 2, 1);
        var a = Panel(d, "(a) Lowest radial internal frequencies", "Radial mode index", "Dimensionless frequency");
        a.Add("Computed l = 0 modes", [1, 2, 3], finest.GetProperty("frequencies").EnumerateArray().Select(z => z.GetDouble()).ToArray(), Blue);
        Line(a, "Continuum |ω|", 1, 3, Math.Abs(J(s, "continuum_threshold")), Gray);
        Line(a, "Upper-bound expression", 1, 3, J(s, "number_weighted_upper_bound"), Orange);
        Line(a, "sqrt(μs μp)", 1, 3, Math.Sqrt(J(finest, "mu_s") * J(finest, "mu_p")), Green);
        a = Panel(d, "(b) Radial-grid sensitivity, R = 24", "Grid spacing h", "Difference from h² extrapolate", true, true);
        for (int k = 0; k < 3; k++)
        {
            var z = runs.Take(3).ToArray();
            var y = z.Select(t => t.GetProperty("frequencies")[k].GetDouble()).ToArray();
            var extrap = (4 * y[^1] - y[^2]) / 3;
            a.Add($"Mode {k + 1}", A(z, "h"), Map(y, v => Math.Abs(v - extrap)), Colors[k + 1]);
        }
        Save(d, "internal_spectrum");
    }

    static (double F, double D, double K, double Q) Coefficients()
    {
        double r = .5 * (Math.Tanh(.5) / .5 + 1 / Math.Pow(Math.Cosh(.5), 2)), c = Math.Pow(Math.Cosh(.5), 2);
        double d = c * (1 - r) * (1 + 11 * r), f = (1 + 6 * r) / (1 + 11 * r), p = 12 * c * r * r * (3 * r - 2) / (1 + 11 * r), k = .15 / p;
        return (f, d, k, 4 * k * d * f * (1 - f));
    }
    static double DeltaB(double theta)
    {
        var c = Coefficients();
        double df = (1 - 2 * c.F) * Math.Pow(Math.Sin(theta), 2) + 2 * Math.Sqrt(c.F * (1 - c.F)) * Math.Cos(theta) * Math.Sin(theta);
        return c.K * c.D * df * df;
    }
    static void Dynamics()
    {
        var rows = Json("dynamics/J960_rtol2e-11/results.json").GetProperty("rows").EnumerateArray().Where(r => J(r, "labels") == 1024).OrderBy(r => J(r, "theta")).ToArray();
        double c = J(rows[0], "C");
        var t = A(rows, "theta");
        var measured = A(rows, "KL");
        var d = new PlotDocument("Synthetic readout: equal masses, zero detuning; t = 0.3", 1, 3);
        var a = Panel(d, "(a) Independently predicted response; no finite-angle fit", "Signed coherent rotation θ (radians)", "D(pθ || p0)", false, true);
        foreach (var (lo, hi) in new[] { (-.205, -.015), (.015, .205) })
        {
            var x = Grid(lo, hi, 301);
            a.Add("Quartic law C4 θ⁴", x, Map(x, v => c * Math.Pow(v, 4)), Orange).Dash = "5,3";
            a.Add("Finite-time pressure predictor", x, Map(x, v => c * Math.Pow(DeltaB(v) / Coefficients().Q, 2)), Green);
        }
        a.Add("Evolved field + noisy readout", t, measured, "#222222", true);
        a.YMin = 3e-11;
        a.YMax = 8e-7;
        var quartic = Map(t, v => c * Math.Pow(v, 4));
        var pressure = Map(t, v => c * Math.Pow(DeltaB(v) / Coefficients().Q, 2));
        Compare("noisy_quartic_predictor", quartic, A(rows, "quartic"), 1e-16);
        Compare("noisy_pressure_predictor", pressure, A(rows, "pressure"), 1e-16);
        var qe = Combine(quartic, measured, (q, m) => 100 * (q / m - 1));
        var pe = Combine(pressure, measured, (q, m) => 100 * (q / m - 1));
        a = Panel(d, "(b) Both prediction residuals: same vertical scale", "Signed coherent rotation θ", "Relative error (%)");
        a.Add($"Quartic: max {qe.Max(Math.Abs):F2}%", t, qe, Orange);
        a.Add($"Pressure: max {pe.Max(Math.Abs):F4}%", t, pe, Green);
        a = Panel(d, "(c) Pressure residual: expanded vertical scale", "Signed coherent rotation θ", "Relative error (%)");
        a.Add("Pressure predictor", t, pe, Green);
        a.YMin = -1.15 * pe.Max(Math.Abs);
        a.YMax = 1.15 * pe.Max(Math.Abs);
        foreach (var p in d.Panels)
        {
            p.XMin = -.21;
            p.XMax = .21;
        }
        Save(d, "noisy_response");
        rows = Json("dynamics/time_scaling.json").GetProperty("rows").EnumerateArray().ToArray();
        d = new PlotDocument("Synthetic readout: equal masses, zero detuning", 2, 1);
        var left = Panel(d, "(a) Finite-angle sixth-order time law", "Readout time t", "D(pθ || p0)", true, true);
        var right = Panel(d, "(b) Pressure-normalized collapse", "Readout time t", "2D/(J* δb² t⁶)");
        foreach (var (angle, color) in new[] { (.05, Green), (.1, Orange), (.2, "#756BB1") })
        {
            var z = rows.Where(r => J(r, "theta") == angle).ToArray();
            var x = A(z, "t");
            left.Add($"θ = {angle}", x, A(z, "KL"), color);
            double half = J(z[0], "J_half");
            left.Add($"θ = {angle}, t⁶ asymptote", x, Map(x, v => half * Math.Pow(DeltaB(angle), 2) * Math.Pow(v, 6)), color).Dash = "3,4";
            foreach (int sign in new[] { -1, 1 })
            {
                z = rows.Where(r => J(r, "theta") == sign * angle).ToArray();
                right.Add($"θ = {sign * angle}", A(z, "t"), A(z, "normalized_ratio"), color).Dash = sign < 0 ? "5,3" : "";
                Compare($"time_normalization_{sign * angle}", z.Select(r => J(r, "KL") / (J(r, "J_half") * Math.Pow(DeltaB(sign * angle), 2) * Math.Pow(J(r, "t"), 6))).ToArray(), A(z, "normalized_ratio"), 1e-10);
            }
        }
        right.YMin = .916;
        right.YMax = 1.009;
        Line(right, "Unit asymptote", .025, .3, 1);
        Save(d, "time_scaling");
    }

    static double XLogX(double x) => x > 0 ? x * Math.Log(x) : 0;
    static double BinaryKl(double p, double q) => XLogX(p) + XLogX(1 - p) - p * Math.Log(q) - (1 - p) * Math.Log(1 - q);
    static void Carrier()
    {
        var r = Read("demonstrations/compression.csv");
        var env = Read("demonstrations/entropy_envelope.csv");
        var amb = Read("demonstrations/mechanical_entropy_ambiguity.csv");
        var x = A(r, "tau");
        var es = Map(x, v => .8 * Math.Exp(2 * v));
        var ea = Map(x, v => .2 * Math.Exp(3 * v));
        var total = Combine(es, ea, (a, b) => a + b);
        var z = Combine(ea, total, (a, b) => a / b);
        var vir = Map(z, v => 2 + v);
        var entropy = Map(z, v => BinaryKl(v, .4));
        Compare("carrier_Es", es, A(r, "E_s_over_E0"));
        Compare("carrier_Ea", ea, A(r, "E_a_over_E0"));
        Compare("carrier_virial", vir, A(r, "virial_over_E"));
        Compare("carrier_entropy", entropy, A(r, "D_to_uniform"));
        var d = new PlotDocument("Exact carrier mechanics and entropy identities", 2, 2);
        var a = Panel(d, "(a) Calibrated channel scaling", "Compression τ = log λ", "Energy relative to initial total", false, true);
        a.Add("Es/E0", x, es, Blue);
        a.Add("Ea/E0", x, ea, Orange);
        a.Add("Eexc/E0", x, total, Gray);
        a = Panel(d, "(b) Virial and register response", "Compression τ", "Normalized virial Vint/Eexc");
        a.Add("2 + zτ", x, vir, Blue);
        var selected = r.Where((_, i) => i % 32 == 0).ToList();
        a.Add("Numerical ∂τ log Z0", A(selected, "tau"), A(selected, "numerical_logZ_derivative"), Orange, true);
        Line(a, "Uniform carrier", x[0], x[^1], 2.4);
        a = Panel(d, "(c) Unrestricted and fixed-core constraints", "Normalized virial Vint/Eexc", "D(ρE || I5/5)");
        var ez = A(env, "z");
        a.Add("Unrestricted sharp bound", Map(ez, v => 2 + v), Map(ez, v => BinaryKl(v, .4)), Blue);
        Compare("carrier_entropy_envelope", Map(ez, v => BinaryKl(v, .4)), A(env, "sharp_lower_bound"));
        var local = env.Where(t => Math.Abs(N(t, "z") - .4) <= .2).ToList();
        a.Add("Local quadratic", A(local, "virial_over_E"), A(local, "local_quadratic"), Gray);
        a.Add("Variable-density isotropic family", A(selected, "virial_over_E"), A(selected, "D_to_uniform"), Orange, true);
        a.Add("Unattained limit", [2], [Math.Log(5.0 / 3)], Blue, true);
        a.Add("Fixed positive core, z = 1", [3], [Math.Log(5)], Red, true);
        a.XMin = 1.98;
        a.XMax = 3.04;
        a.YMin = -.04;
        a.YMax = 1.76;
        a = Panel(d, "(d) Fixed mechanics, different entropy", "Contact orientation 2θ/π", "D(ρE || I5/5)");
        var ent = Map(A(amb, "theta_radians"), v => Math.Log(5) + XLogX((1 - Math.Cos(v)) / 2) + XLogX((1 + Math.Cos(v)) / 2));
        Compare("carrier_ambiguity", ent, A(amb, "D_formula"));
        a.Add("Same n, j, R; Es = 0, Ea = Eexc", A(amb, "theta_over_pi_half"), ent, Red);
        Line(a, "log 5", 0, 1, Math.Log(5));
        Line(a, "log(5/2)", 0, 1, Math.Log(2.5));
        Save(d, "carrier_mechanics_entropy");
    }

    static void Geometry()
    {
        var info = Json("equilibrium_family/results/validation.json");
        var curves = Read("equilibrium_family/results/normalized_shapes.csv");
        var d = new PlotDocument("Isolated stationary mixed sector: zero detuning, no excess internal stress", 2, 2);
        var den = Panel(d, "(a) Density shape", "r/R50", "ρ/ρc", false, true);
        den.XMin = 0;
        den.XMax = 5;
        den.YMin = 1e-5;
        den.YMax = 1.3;
        var speed = Panel(d, "(b) Circular-speed shape", "r/R50", "vc/sqrt(GM/R50)");
        speed.XMin = 0;
        speed.XMax = 5;
        speed.YMin = 0;
        speed.YMax = 1.03;
        var sharedDen = Panel(d, "(c) Shared coefficients; different occupations", "r/R*", "ρ/ρ*", false, true);
        sharedDen.XMin = 0;
        sharedDen.XMax = 12;
        sharedDen.YMin = 1e-5;
        sharedDen.YMax = 15;
        var sharedSpeed = Panel(d, "(d) Speeds at the same shared coefficients", "r/R*", "vc/V*");
        sharedSpeed.XMin = 0;
        sharedSpeed.XMax = 12;
        int index = 0;
        foreach (var rec in info.GetProperty("runs").EnumerateArray())
        {
            double eta = J(rec, "eta"), r50 = J(rec, "R50"), mass = J(rec, "radial_mass_integral");
            string et = eta.ToString("G", CultureInfo.InvariantCulture), color = Colors[index++], label = eta == 0 ? "η = 0 limit" : $"η = {et}";
            var p = Read($"equilibrium_family/results/profile_eta_{et}.csv");
            var x = A(p, "x");
            var fp = new MonotoneCubic(x, A(p, "f"));
            var pp = new MonotoneCubic(x, A(p, "phi_prime"));
            var c = curves.Where(t => N(t, "eta") == eta).ToList();
            var xx = A(c, "r_over_R50");
            Compare("normalized_speed_eta_" + et, Map(xx, v => Math.Sqrt(Math.Max(v * r50 * pp.At(v * r50), 0) / (mass / r50))), A(c, "vc_over_sqrt_GM_R50"), 1e-7);
            den.Add(label, xx, A(c, "rho_over_rhoc"), color);
            speed.Add(label, xx, A(c, "vc_over_sqrt_GM_R50"), color);
            if (eta > 0)
            {
                var y = Grid(0, 12, 1201);
                sharedDen.Add(label, y, Map(y, v => Math.Max(eta * eta * Math.Pow(fp.At(v * Math.Sqrt(eta)), 2), 1e-14)), color);
                sharedSpeed.Add(label, y, Map(y, v => Math.Sqrt(Math.Max(eta * v * Math.Sqrt(eta) * pp.At(v * Math.Sqrt(eta)), 0))), color);
            }
        }
        var kepler = Grid(1, 5, 401);
        speed.Add("Kepler asymptote x^-1/2", kepler, Map(kepler, v => 1 / Math.Sqrt(v)), Gray).Dash = "3,3";
        Save(d, "equilibrium_shape_family");
        info = Json("baryonic_geometry/results/validation.json");
        double q = J(info.GetProperty("parameters"), "q"), number = J(info.GetProperty("parameters"), "number");
        double[] fractions = [.1, 1, 3];
        var comparisons = new Dictionary<(double, double), List<Dictionary<string, string>>>();
        foreach (double fraction in fractions)
            foreach (double flat in new[] { 0, .5, .8 })
            {
                string f = fraction.ToString("G", CultureInfo.InvariantCulture), a = flat.ToString("G", CultureInfo.InvariantCulture);
                var c = Read($"baryonic_geometry/results/comparison_m{f}_a{a}.csv");
                comparisons[(fraction, flat)] = c;
                var radii = A(c, "r");
                Compare($"geometry_baryon_force_{f}_{a}", Map(radii, r => fraction * number / (4 * Math.PI) * r * r / Math.Pow(r * r + q * q, 1.5)), A(c, "vc2_baryon"), 1e-11);
                Compare($"geometry_force_identity_{f}_{a}", Map(A(c, "vc_total_deformed"), v => v * v), Combine(A(c, "vc2_field_deformed"), A(c, "vc2_baryon"), (v, b) => v + b), 1e-11);
            }
        var reference = comparisons[(1, 0)];
        double peak = Math.Sqrt(A(reference, "vc2_field_isolated").Max());
        d = new PlotDocument("Controlled baryonic sources, q = r0: fixed field mass and coefficients; Es = Ea = 0", 2, 2);
        var panel = Panel(d, "(a) Fixed field mass, Mb/Mf = 1", "R/r0", "vc,total / vpeak,isolated");
        panel.Add("Unchanged field + baryons", Map(A(reference, "r"), v => v / q), Map(A(reference, "vc_total_unchanged_field"), v => v / peak), Gray);
        foreach (var (flat, label, color) in new[] { (0.0, "Sphere: field readjusted", Blue), (.8, "Disk: field readjusted", Orange) })
        {
            var c = comparisons[(1, flat)];
            panel.Add(label, Map(A(c, "r"), v => v / q), Map(A(c, "vc_total_deformed"), v => v / peak), color);
        }
        var response = Panel(d, "(b) Readjustment of the field density", "R/r0", "Change from unchanged-field speed (%)");
        var flattening = Panel(d, "(c) Field flattening after readjustment", "Baryon-to-field mass Mb/Mf", "sqrt(2 <z²>/<R²>)", true);
        var difference = Panel(d, "(d) Same baryonic midplane force", "R/r0", "Disk minus sphere speed (%)");
        for (int k = 0; k < fractions.Length; k++)
        {
            double f = fractions[k];
            string color = new[] { Green, Blue, Orange }[k];
            foreach (double flat in new[] { 0, .8 })
            {
                var c = comparisons[(f, flat)];
                response.Add($"Mb/Mf = {f}, {(flat == 0 ? "sphere" : "disk")}", Map(A(c, "r"), v => v / q), Combine(A(c, "vc_total_deformed"), A(c, "vc_total_unchanged_field"), (v, b) => 100 * (v / b - 1)), color).Dash = flat == 0 ? "5,3" : "";
            }
            foreach (double flat in new[] { .5, .8 })
            {
                var c = comparisons[(f, flat)];
                difference.Add($"Mb/Mf = {f}, a/q = {flat}", Map(A(c, "r"), v => v / q), Combine(A(c, "vc_total_deformed"), A(c, "vc_total_spherical"), (v, b) => 100 * (v / b - 1)), color).Dash = flat == .5 ? "5,3" : "";
            }
        }
        index = 0;
        foreach (double flat in new[] { 0, .5, .8 })
        {
            var z = info.GetProperty("fine_solutions").EnumerateArray().Where(t => J(t, "flattening") == flat && J(t, "fraction") > 0).OrderBy(t => J(t, "fraction")).ToArray();
            flattening.Add(flat == 0 ? "Sphere" : $"Disk, a/q = {flat}", A(z, "fraction"), A(z, "rms_axis_ratio"), Colors[index++]);
        }
        foreach (var p in new[] { panel, response, difference })
        {
            p.XMin = 0;
            p.XMax = 5;
        }
        Save(d, "baryonic_geometry_backreaction");
        info = Json("baryonic_geometry/scale_sensitivity/scale_sensitivity.json");
        q = J(info.GetProperty("parameters"), "isolated_R50");
        d = new PlotDocument("Prescribed baryons: identical midplane force within each pair; zero excess stress", 2, 1);
        foreach (double scale in new[] { .25, 4 })
        {
            panel = Panel(d, $"q/r0 = {scale}, a/q = 0.8", "R/r0", "Disk minus sphere speed (%)");
            panel.XMin = 0;
            panel.XMax = 5;
            for (int k = 0; k < fractions.Length; k++)
            {
                double f = fractions[k];
                var c = Read($"baryonic_geometry/scale_sensitivity/qratio_{scale.ToString("G", CultureInfo.InvariantCulture)}/comparison_m{f.ToString("G", CultureInfo.InvariantCulture)}_a0.8.csv");
                var diff = Combine(A(c, "vc_total"), A(c, "vc_total_sphere"), (v, b) => 100 * (v / b - 1));
                Compare($"scale_difference_{scale}_{f}", diff, A(c, "disk_minus_sphere_percent"));
                panel.Add($"Mb/Mf = {f}", Map(A(c, "r"), v => v / q), diff, new[] { Green, Blue, Orange }[k]);
            }
        }
        Save(d, "baryonic_scale_sensitivity");
    }

    internal sealed class MonotoneCubic(double[] x, double[] y)
    {
        readonly CubicSpline spline = CubicSpline.InterpolatePchipSorted(x, y);

        public double At(double z)
        {
            if (z < x[0] || z > x[^1])
                return double.NaN;
            int i = Array.BinarySearch(x, z);
            return i >= 0 ? y[i] : spline.Interpolate(z);
        }
    }

    static void Thermal()
    {
        var cubic = Read("cubic_nonequilibrium/results/cubic_response.csv");
        var thermal = Read("cubic_nonequilibrium/results/fixed_hamiltonian_high_temperature.csv");
        var kin = Read("cubic_nonequilibrium/results/positive_kinetic_example.csv");
        var d = new PlotDocument("Cubic response, fixed-Hamiltonian thermal control, and positive kinetic example", 2, 2);
        var p = Panel(d, "(a) Cubic compression response", "Contact fraction z", "Cexp(D,D,D)");
        p.Add("Cubic compression response", A(cubic, "z"), A(cubic, "cubic"), Blue);
        Compare("cubic_compression", Map(A(cubic, "z"), v => v * (1 - v) * (1 - 2 * v)), A(cubic, "cubic"));
        Line(p, "Upper extreme", 0, 1, 1 / (6 * Math.Sqrt(3)));
        Line(p, "Lower extreme", 0, 1, -1 / (6 * Math.Sqrt(3)));
        p.Add("Maximum entropy", [.4], [6.0 / 125], Orange, true);
        p = Panel(d, "(b) Fixed Hamiltonian, high T", "kB T / energy unit", "Free energy / energy unit", true, true);
        p.Add("Exact beyond quadratic", A(thermal, "dimensionless_kBT"), A(thermal, "beyond_quadratic_free_energy"), Blue);
        p.Add("Cubic term κ3/(6T²)", A(thermal, "dimensionless_kBT"), A(thermal, "leading_cubic"), Orange);
        var r = Map(Grid(-1, 2, 301), v => Math.Pow(10, v));
        var rows = kin.Where(t => N(t, "r_over_a") > 0).ToList();
        p = Panel(d, "(c) Positive kinetic example: variance", "Radius r/a", "Normalized velocity variance", true);
        p.YMin = -.02;
        p.YMax = 1.05;
        p.Add("Teff/Teff(0)", r, Map(r, v => 1 / Math.Sqrt(1 + v * v)), Blue);
        p.Add("Quadrature", A(rows, "r_over_a"), Map(A(rows, "variance_quadrature"), v => v / N(kin[0], "variance_quadrature")), Blue, true);
        p = Panel(d, "(d) Positive kinetic example: entropy", "Radius r/a", "Specific entropy rise ΔSB/kB", true);
        p.Add("(7/4) log(1+r²)", r, Map(r, v => 1.75 * Math.Log(1 + v * v)), Orange);
        p.Add("Quadrature", A(rows, "r_over_a"), A(rows, "entropy_change"), Orange, true);
        Compare("kinetic_entropy", Map(A(kin, "r_over_a"), v => 1.75 * Math.Log(1 + v * v)), A(kin, "entropy_change"));
        Save(d, "cubic_thermal_controls");
        var data = Read("cubic_nonequilibrium/results/texture_response.csv");
        d = new PlotDocument("Controlled nonequilibrium texture response", 2, 1);
        var a = Panel(d, "(a) Short-time sign and convergence", "Time (reference units)", "10³ Δvf²/t²", true);
        var b = Panel(d, "(b) Small, finite-time redistribution", "Time (reference units)", "10⁶ Δvf² (reference units)");
        foreach (var (sign, color, label) in new[] { (1, Blue, "Positive jet"), (-1, Orange, "Negative jet") })
            foreach (int grid in new[] { 480, 960, 1920 })
            {
                var z = data.Where(row => N(row, "sign") == sign && N(row, "J") == grid).ToList();
                var time = A(z, "t");
                a.Add($"{label}, {grid}", time, Map(A(z, "delta_vf_squared_over_t_squared"), v => 1e3 * v), color).Dash = grid == 480 ? "2,3" : grid == 960 ? "6,3" : "";
                Compare($"texture_scaling_{sign}_{grid}", z.Select(v => N(v, "delta_vf_squared") / Math.Pow(N(v, "t"), 2)).ToArray(), A(z, "delta_vf_squared_over_t_squared"));
                if (grid == 1920)
                {
                    Line(a, label + " continuum coefficient", time[0], time[^1], 1e3 * N(z[0], "continuum_coefficient"), color);
                    b.Add(label, [0, .. time], [0, .. Map(A(z, "delta_vf_squared"), v => 1e6 * v)], color);
                }
            }
        Save(d, "cubic_texture_response");
    }

    static void TwoField()
    {
        var f = Read("baryonic_two_field/results/radial_profiles.csv").Where(r => r["run"] == "J959_R24").ToList();
        var e = Read("baryonic_two_field/results/energies.csv").Where(r => r["run"] == "J959_R24").ToList();
        var t = Read("baryonic_two_field/results/tracers.csv").Where(r => r["run"] == "J959_R24").ToList();
        var proto = Json("baryonic_two_field/results/protocol.json");
        var validation = Json("baryonic_two_field/results/validation.json");
        if (validation.GetProperty("gates").EnumerateObject().Any(g => !g.Value.GetBoolean()))
            throw new InvalidDataException("Saved two-field validation gates are not passing.");
        List<Dictionary<string, string>> Profile(string name, double time) => f.Where(row => row["case"] == name && Math.Abs(N(row, "t") - time) < 1e-8 && N(row, "r") <= 8).OrderBy(row => N(row, "r")).DistinctBy(row => N(row, "r")).ToList();
        var iso = Profile("isolated_mixed", 0);
        var coupled = Profile("coupled_mixed", 0);
        var r = A(iso, "r");
        double number = J(proto, "number"), scale = J(proto, "baryon_scale");
        var vb = Map(r, v => number / (4 * Math.PI) * v * v / Math.Pow(v * v + scale * scale, 1.5));
        var d = new PlotDocument("Controlled spherical two-field experiment: prescribed baryons, fixed field mass", 3, 1);
        var p = Panel(d, "(a) Compression at fixed field mass", "r/L", "Field density / ρ0", false, true);
        p.YMin = 1e-6;
        p.YMax = 4;
        p.Add("Isolated equilibrium", r, A(iso, "rho"), Gray);
        p.Add("Baryon-coupled equilibrium", r, A(coupled, "rho"), Blue);
        p = Panel(d, "(b) Redistribution, no extra mass", "r/L", "Mf(<r)/Mf");
        p.YMin = 0;
        p.YMax = 1.04;
        p.Add("Isolated equilibrium", r, Map(A(iso, "mass_field"), v => v / number), Gray);
        p.Add("Baryon-coupled equilibrium", r, Map(A(coupled, "mass_field"), v => v / number), Blue);
        p = Panel(d, "(c) Gravity after equilibrium response", "r/L", "sqrt(r ∂r Φ)/V");
        p.Add("Isolated field", r, Map(A(iso, "vf2"), Math.Sqrt), Gray).Dash = "2,3";
        p.Add("Frozen isolated field + baryons", r, Combine(A(iso, "vf2"), vb, (v, b) => Math.Sqrt(v + b)), Gray).Dash = "6,3";
        p.Add("Coupled field + baryons", r, Combine(A(coupled, "vf2"), vb, (v, b) => Math.Sqrt(v + b)), Blue);
        Save(d, "baryonic_two_field_background");
        d = new PlotDocument("Controlled field response: instantaneous force, internal energy", 2, 2);
        p = Panel(d, "(a) Field-force change at t = 0.3", "r/L", "Δ(r ∂r Φf)/V²");
        foreach (var (name, reference, color, label) in new[] { ("baryon_quench", "isolated_mixed", Orange, "Baryonic switch, mixed field"), ("isolated_rotated", "isolated_mixed", "#8055a1", "Internal rotation, isolated"), ("coupled_rotated", "coupled_mixed", "#bd3634", "Internal rotation, coupled") })
        {
            var a = Profile(name, .3);
            var b = Profile(reference, .3);
            Compare($"two_field_aligned_radii_{name}", A(a, "r"), A(b, "r"));
            p.Add(label, A(a, "r"), Combine(A(a, "vf2"), A(b, "vf2"), (v, w) => v - w), color);
        }
        var es = Panel(d, "(b) Internal gradient energy", "t/(L/V)", "Es/E0");
        var ea = Panel(d, "(c) Contact excess energy", "t/(L/V)", "Ea/E0");
        foreach (var (name, color, label) in new[] { ("isolated_rotated", "#8055a1", "Isolated, imposed rotation"), ("coupled_rotated", "#bd3634", "Coupled, imposed rotation") })
        {
            var z = e.Where(row => row["case"] == name).ToList();
            es.Add(label, A(z, "t"), A(z, "Es"), color);
            ea.Add(label, A(z, "t"), A(z, "Ea"), color);
        }
        Save(d, "baryonic_two_field_response");
        d = new PlotDocument("Passive tracer: L = 1 kpc, V = 30 km/s; R0 = 1.01169 L; imposed θ = 0.2", 3, 1);
        double unit = J(validation.GetProperty("physical_unit_illustration"), "time_unit_Myr");
        var time = Map(A(t, "t"), v => v * unit);
        foreach (var (key, jet, factor, title, ylabel) in new[] { ("delta_radius", "radius_jet", 1000.0, "(a) Radial displacement", "ΔR [pc]"), ("delta_radial_speed", "radial_speed_jet", 30000.0, "(b) Radial velocity", "ΔvR [m/s]"), ("delta_tangential_speed", "tangential_speed_jet", 30000.0, "(c) Tangential velocity", "Δvφ [m/s]") })
        {
            p = Panel(d, title, "Time [Myr]", ylabel);
            p.Add("Passive tracer, evolved field", time, Map(A(t, key), v => factor * v), Blue);
            p.Add("Leading small-time term", time, Map(A(t, jet), v => factor * v), Gray).Dash = "6,3";
        }
        Save(d, "baryonic_two_field_tracer");
    }

    static void Synthetic()
    {
        var saved = Json("core/results/core_diagnostics.json");
        var fit = saved.GetProperty("mock_fit");
        var mock = Read("core/results/mock_circular_speed.csv");
        var curve = Read("core/results/circular_speed_fit.csv");
        var core = Read("core/results/core_profile.csv");
        var d = new PlotDocument("Synthetic core recovery: fixed shape αh = 0.3, βh = 0.5", 2, 2);
        var p = Panel(d, "(a) Circular speed: synthetic recovery", "Radius [kpc]", "Circular speed [km/s]");
        p.XMin = 0;
        p.XMax = 8;
        p.YMin = 0;
        p.YMax = 31;
        var x = A(curve, "r_kpc");
        var prediction = A(curve, "fit_v_kms");
        var se = A(curve, "fit_se_kms");
        p.Add("Generating core", x, A(curve, "true_v_kms"), Gray).Dash = "5,3";
        p.Band("Fit mean: 1 s.e.", x, Combine(prediction, se, (v, e) => v - e), Combine(prediction, se, (v, e) => v + e), Blue);
        p.Add("Two-scale fit", x, prediction, Blue);
        var density = Panel(d, "(b) Solved density profile", "Radius [kpc]", "Density / central density", false, true);
        density.XMin = 0;
        density.XMax = 8;
        density.YMin = 1e-5;
        density.YMax = 1.2;
        density.Add("n/n(0) = u²", A(core, "r_dimensionless"), Map(A(core, "u"), v => v * v), Green);
        double r99 = J(saved.GetProperty("physical"), "R99_kpc");
        density.Add("R99", [r99, r99], [1e-5, 1.2], Gray).Dash = "3,3";
        var residual = Panel(d, "(c) Synthetic residuals", "Radius [kpc]", "Data - fit [km/s]");
        residual.XMin = 0;
        residual.XMax = 8;
        residual.YMin = -3.5;
        residual.YMax = 3.5;
        foreach (var (training, color, label) in new[] { (true, Blue, "Training mock data"), (false, Orange, "Held-out mock data") })
        {
            var z = mock.Where(row => (N(row, "is_training") != 0) == training).ToList();
            var rad = A(z, "r_kpc");
            var noise = Map(rad, _ => J(fit, "sigma_kms"));
            p.ErrorBars(label, rad, A(z, "mock_v_kms"), noise, color);
            residual.ErrorBars(label, rad, Combine(A(z, "mock_v_kms"), A(z, "fit_v_kms"), (a, b) => a - b), noise, color);
        }
        double trainChi = mock.Where(row => N(row, "is_training") != 0).Sum(row => Math.Pow((N(row, "mock_v_kms") - N(row, "fit_v_kms")) / J(fit, "sigma_kms"), 2));
        double rmse = Math.Sqrt(mock.Where(row => N(row, "is_training") == 0).Average(row => Math.Pow(N(row, "mock_v_kms") - N(row, "fit_v_kms"), 2)));
        Check("synthetic_train_chisq", Math.Abs(trainChi - J(fit, "train_chisq")), 1e-10);
        Check("synthetic_holdout_rmse", Math.Abs(rmse - J(fit, "holdout_rmse_kms")), 1e-10);
        var metrics = Panel(d, "Synthetic recovery summary", "", "");
        metrics.TextOnly = true;
        metrics.Note = $"Mock truth: L = 1 kpc, vs = 30 km/s\nFit: L = {J(fit, "L_kpc"):F3} ± {J(fit, "L_se"):F3} kpc\n       vs = {J(fit, "vs_kms"):F2} ± {J(fit, "vs_se"):F2} km/s\nTraining χ² / dof = {trainChi:F2} / 22\nHeld-out RMSE = {rmse:F2} km/s\nFixed shape: αh = 0.3, βh = 0.5";
        Save(d, "synthetic_core_fit");
    }

    static void Register()
    {
        var rows = Csv.Read(Path.Combine(source, "reference/numerics/compression.csv"));
        var cycles = Csv.Read(Path.Combine(source, "reference/numerics/cycle_cost.csv"));
        var x = A(rows, "tau");
        double threshold = Math.Log(8.0 / 3);
        var z = Map(x, t => .2 * Math.Exp(3 * t) / (.8 * Math.Exp(2 * t) + .2 * Math.Exp(3 * t)));
        var gap = Map(z, v => BinaryKl(v, .2));
        var work = Map(z, v => Math.Max(0, (5 * v - 2) * threshold) / 3);
        var mismatch = Combine(gap, work, (a, b) => a - b);
        Compare("register_relative_entropy", gap, A(rows, "D_to_reference"));
        Compare("register_ergotropy", work, A(rows, "ergotropy_over_kBT"));
        Compare("register_spectral_mismatch", mismatch, A(rows, "spectral_mismatch"));
        var length = A(cycles, "BKM_distance");
        var bound = Map(length, v => 4 * Math.Sin(v / 2) * Math.Atanh(Math.Sin(v / 2)));
        var local = Map(length, v => v * v);
        double Cycle(double v, double k)
        {
            double shift = k * (Math.PI - v) / 4, a0 = Math.PI / 4 - v / 4 + shift, a1 = Math.PI / 4 + v / 4 + shift, z0 = Math.Pow(Math.Sin(a0), 2), z1 = Math.Pow(Math.Sin(a1), 2);
            return (Math.Log(z1 / (1 - z1)) - Math.Log(z0 / (1 - z0))) * (z1 - z0);
        }
        var asymmetric = Map(length, v => Cycle(v, .7));
        Compare("register_sharp_cycle_bound", bound, A(cycles, "sharp_bound"));
        Compare("register_symmetric_cycle", Map(length, v => Cycle(v, 0)), A(cycles, "symmetric_cycle_direct"));
        Compare("register_asymmetric_cycle", asymmetric, A(cycles, "asymmetric_cycle_direct"));
        var d = new PlotDocument("Independent finite-register example: deterministic exact identities", 2, 1);
        var p = Panel(d, "(a) Register work and activation", "Control τ", "Work / kBT");
        p.YMin = -.035;
        p.YMax = 1.17;
        p.Add("Free-energy gap D(ρτ || ρ0)", x, gap, Blue);
        p.Add("Extractable unitary work", x, work, Orange);
        p.Add("Spectral mismatch", x, mismatch, Green).Dash = "5,3";
        p.Add("τc = log(8/3)", [threshold, threshold], [0, .8], Gray).Dash = "3,3";
        p = Panel(d, "(b) Exact two-quench cycle cost", "Full BKM endpoint distance L", "Cycle work / kBT");
        p.XMin = 0;
        p.XMax = 3;
        p.YMin = 0;
        p.YMax = 16;
        p.Add("Asymmetric endpoints (κ = 0.7)", length, asymmetric, Orange);
        p.Add("Sharp bound; z0 + z1 = 1", length, bound, Blue);
        p.Add("Local bound L²", length, local, Gray).Dash = "5,3";
        Save(d, "register_work");
    }
}



