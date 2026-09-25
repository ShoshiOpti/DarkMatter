using MathNet.Numerics;
using System.Text.Json;
using MathNet.Numerics.Interpolation;

namespace DarkUniverse;

internal static class LegacyPilot
{
    static double N(JsonElement value, string key) => value.GetProperty(key).GetDouble();
    static double[] A(JsonElement value, string key) => value.GetProperty(key).EnumerateArray().Select(v => v.GetDouble()).ToArray();

    public static string Run(string sourceRoot, string outputDir)
    {
        string folder = Path.Combine(sourceRoot, "numerics/coupled_pilot");
        using var savedFile = JsonDocument.Parse(File.ReadAllText(Path.Combine(folder, "results.json")));
        using var solutionFile = JsonDocument.Parse(File.ReadAllText(Path.Combine(folder, "legacy_dense_solution.json")));
        var saved = savedFile.RootElement.GetProperty("models");
        var solutions = solutionFile.RootElement.GetProperty("solutions").EnumerateArray().ToDictionary(s => (s.GetProperty("split").GetString()!, s.GetProperty("galaxy").GetString()!));
        var document = new PlotDocument("Conditional spherical-source pilot: four shared/local parameters per comparison", 2, 2) { Scope = "Conditional two-galaxy spherical pilot; saved-parameter field solutions, not a population result" };
        var checks = new List<object>();
        foreach (string name in new[] { "F583-1", "UGC00731" })
            foreach (string split in new[] { "inner", "full" })
            {
                var g = saved.GetProperty(split).GetProperty("galaxies").GetProperty(name);
                var s = solutions[(split, name)];
                var radii = A(g, "r");
                var observed = A(g, "observed");
                var errors = A(g, "error");
                var baryon = A(g, "baryon_v2");
                var mass = radii.Zip(baryon, (r, b) => r * b).ToArray();
                var source = new ControlledPlots.MonotoneCubic(radii, mass);
                double length = N(s, "length"), velocity2 = N(s, "velocity2");
                var knots = A(s, "knots");
                var coefficients = s.GetProperty("force_coefficients").EnumerateArray().Select(p => p.EnumerateArray().Select(v => v.GetDouble()).ToArray()).ToArray();
                var force = new CubicSpline(knots, coefficients[3], coefficients[2], coefficients[1], coefficients[0]);
                double Mass(double r) => r < radii[0] ? mass[0] * Math.Pow(r / radii[0], 3) : r > radii[^1] ? mass[^1] : source.At(r);
                double Predict(double r)
                {
                    double y = r / length, shape = y <= knots[^1] ? y * force.Interpolate(y) : knots[^1] * knots[^1] * force.Interpolate(knots[^1]) / y;
                    return Math.Sqrt(Mass(r) / r + velocity2 * shape);
                }
                var nf = g.GetProperty("nfw");
                double Nfw(double r) => Math.Sqrt(Mass(r) / r + Math.Pow(N(nf, "amplitude_kms"), 2) * NfwShape(r / N(nf, "length_kpc")));
                double nodeError = radii.Select(Predict).Zip(A(g, "coupled_prediction"), (a, b) => Math.Abs(a - b)).Max();
                double nfwError = radii.Select(Nfw).Zip(A(g, "nfw_prediction"), (a, b) => Math.Abs(a - b)).Max();
                if (nodeError > 1e-7 || nfwError > 1e-8)
                    throw new InvalidDataException($"Legacy pilot nodal mismatch: {split}/{name}: coupled {nodeError:G17}, NFW {nfwError:G17}");
                checks.Add(new
                {
                    split,
                    galaxy = name,
                    dense_points = 500,
                    coupled_nodal_max_absolute_error_kms = nodeError,
                    nfw_nodal_max_absolute_error_kms = nfwError
                });
                var rr = Enumerable.Range(0, 500).Select(i => radii[0] + (radii[^1] - radii[0]) * i / 499).ToArray();
                var panel = new PlotPanel(name + " | " + (split == "inner" ? "inner-only fit" : "all-radius fit"), "Radius [kpc]", "Circular speed [km/s]") { YMin = 0 };
                document.Panels.Add(panel);
                panel.Add("Coupled field + baryons", rr, rr.Select(Predict).ToArray(), "#21608b");
                panel.Add("NFW + baryons", rr, rr.Select(Nfw).ToArray(), "#c05e28");
                panel.Add("Baryons", rr, rr.Select(r => Math.Sqrt(Mass(r) / r)).ToArray(), "#777777").Dash = "3,3";
                int training = g.GetProperty("train_count").GetInt32();
                panel.ErrorBars("Training", radii[..training], observed[..training], errors[..training], "#222222");
                if (training < radii.Length)
                {
                    panel.ErrorBars("Withheld outer", radii[training..], observed[training..], errors[training..], "#222222").Hollow = true;
                    panel.Shades.Add(((radii[training - 1] + radii[training]) / 2, radii[^1], "#f0eee9"));
                }
            }
        string outFolder = Path.Combine(outputDir, "legacy/coupled_pilot");
        Directory.CreateDirectory(outFolder);
        string path = Path.Combine(outFolder, "coupled_pilot.svg");
        document.Save(path);
        File.WriteAllText(Path.Combine(outFolder, "validation.json"), JsonSerializer.Serialize(new
        {
            scope = document.Scope,
            source = "legacy_dense_solution.json",
            checks
        }, new JsonSerializerOptions { WriteIndented = true }));
        return path;
    }


    static double NfwShape(double x) => x < 1e-3 ? x / 2 - 2 * x * x / 3 + 3 * x * x * x / 4 - 4 * Math.Pow(x, 4) / 5 + 5 * Math.Pow(x, 5) / 6 : (SpecialFunctions.Log1p(x) - x / (1 + x)) / x;
}


