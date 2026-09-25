using System.Globalization;
using System.Text.Json;

namespace DarkUniverse;

public sealed record SparcCheckReport(bool Passed, int Checks, string[] Failures,
    double MaximumAbsolutePointError, int ComparedPointRows, int ComparedMetadataRows, int ComparedHoldoutGalaxies);

public static class SparcChecks
{
    public static SparcCheckReport Run(SparcDataset data, string dataRoot, string outputRoot)
    {
        var failures = new List<string>();
        int checks = 0;
        double maxError = 0;
        void Check(bool condition, string label)
        {
            checks++;
            if (!condition)
                failures.Add(label);
        }
        void Close(double actual, double expected, string label)
        {
            double error = Math.Abs(actual - expected);
            maxError = Math.Max(maxError, error);
            Check(double.IsFinite(actual) && error <= 5e-10 + 5e-13 * Math.Abs(expected), label);
        }
        double Number(string text) => double.Parse(text, NumberStyles.Float, CultureInfo.InvariantCulture);
        Check(data.CatalogueSha256 == "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243", "Raw catalogue SHA256");
        Check(data.MassModelsSha256 == "9108994b12cc401b94a1768beca61c53ec354779385c9c9cc571049f3043244c", "Raw mass-model SHA256");
        Check(data.Catalogue.Count == 175, "175 catalogue galaxies");
        Check(data.Points.Values.Sum(p => p.Length) == 3391, "3391 catalogue radii");
        Check(data.SelectedNames.Length == 153, "153 selected galaxies");
        Check(data.SelectedPoints.Count() == 3168, "3168 selected radii");
        Check(data.Points.Values.Sum(p => p.Count(r => r.GasSpeedKms < 0)) == 361, "361 signed negative gas rows");
        Check(data.SelectedPoints.Count(p => p.BaryonsSquared() < 0) == 2, "Two undefined baryon-only circular speeds");
        Check(data.SelectedNames.Count(n => Sparc.HoldoutCount(data.Points[n].Length) > 0) == 131, "131 holdout galaxies");
        Check(data.SelectedNames.Sum(n => Sparc.HoldoutCount(data.Points[n].Length)) == 659, "659 held-out radii");
        foreach (var (radius, value) in new[] { (.08, -33.0975), (.23, -27.92125) })
        {
            var point = data.Points["UGC01281"].Single(p => p.RadiusKpc == radius);
            Close(point.BaryonsSquared(), value, "UGC01281 signed B at " + radius);
        }
        string referenceRoot = "numerics/observations/results/";
        using var metadata = JsonDocument.Parse(File.ReadAllText(Sparc.FindFile(dataRoot, referenceRoot + "selected_metadata.json")));
        var savedMeta = metadata.RootElement.EnumerateArray().ToArray();
        Check(savedMeta.Select(m => m.GetProperty("name").GetString()).SequenceEqual(data.SelectedNames), "Selected metadata order");
        foreach (var expected in savedMeta)
        {
            string name = expected.GetProperty("name").GetString()!;
            var m = data.Catalogue[name];
            Check(m.Type == expected.GetProperty("type").GetInt32(), name + " type");
            Check(m.Quality == expected.GetProperty("Q").GetInt32(), name + " quality");
            Check(m.References == expected.GetProperty("refs").GetString(), name + " references");
            foreach (var (key, actual) in new[] { ("D", m.DistanceMpc), ("e_D", m.DistanceErrorMpc),
                ("inc", m.InclinationDeg), ("e_inc", m.InclinationErrorDeg), ("L36", m.Luminosity36) })
                Close(actual, expected.GetProperty(key).GetDouble(), name + " " + key);
        }
        string[] pointLines = File.ReadAllLines(Sparc.FindFile(dataRoot, referenceRoot + "observed_points_and_predictions.csv"));
        var actualPoints = data.SelectedPoints.ToArray();
        Check(pointLines.Length - 1 == actualPoints.Length, "Saved and imported point counts match");
        int comparedPoints = Math.Min(pointLines.Length - 1, actualPoints.Length);
        for (int i = 0; i < comparedPoints; i++)
        {
            var expected = pointLines[i + 1].Split(',');
            var p = actualPoints[i];
            Check(p.Galaxy == expected[0], "Point galaxy row " + i);
            Close(p.RadiusKpc, Number(expected[1]), p.Galaxy + " radius");
            Close(p.ObservedSpeedKms, Number(expected[2]), p.Galaxy + " observed speed");
            Close(p.ErrorSpeedKms, Number(expected[3]), p.Galaxy + " speed uncertainty");
            double b = p.BaryonsSquared();
            Close(b, Number(expected[4]), p.Galaxy + " signed baryons squared");
            if (b < 0)
                Check(expected[5].Length == 0, p.Galaxy + " undefined circular speed retained");
            else
                Close(Math.Sqrt(b), Number(expected[5]), p.Galaxy + " baryon circular speed");
        }
        using var summary = JsonDocument.Parse(File.ReadAllText(Sparc.FindFile(dataRoot, referenceRoot + "fit_summary.json")));
        Check(data.DisplayNames.SequenceEqual(summary.RootElement.GetProperty("displayed").EnumerateArray().Select(x => x.GetString())), "Six predeclared display galaxies");
        // Reproduce the original descriptive baryon summary over valid circular speeds.
        foreach (double ml in new[] { .5, .3, .7 })
        {
            int invalid = 0;
            var losses = new List<double>();
            foreach (string name in data.SelectedNames)
            {
                double sum = 0;
                int valid = 0;
                foreach (var p in data.Points[name])
                {
                    double b = p.BaryonsSquared(ml, 1.4 * ml);
                    if (b < 0)
                    {
                        invalid++;
                        continue;
                    }
                    double residual = (Math.Sqrt(b) - p.ObservedSpeedKms) / p.ErrorSpeedKms;
                    sum += residual * residual;
                    valid++;
                }
                losses.Add(sum / valid);
            }
            var expected = summary.RootElement.GetProperty("summary").GetProperty(ml.ToString("0.0", CultureInfo.InvariantCulture));
            Check(invalid == expected.GetProperty("baryons_no_circular_orbit_rows").GetInt32(), "Undefined baryon rows for M/L=" + ml);
            Close(Sparc.Median(losses), expected.GetProperty("baryon_median_chi2_per_valid_row").GetDouble(), "Baryon median loss for M/L=" + ml);
        }
        int comparedHoldouts = 0;
        var holdoutNames = new HashSet<string>(StringComparer.Ordinal);
        foreach (string line in File.ReadLines(Sparc.FindFile(dataRoot, "numerics/observations/uncertainty/penalized_fits.csv")).Skip(1))
        {
            string[] fields = line.Split(',');
            if (fields[1] != "core" || fields[2] != "fixed" || fields[3] != "True")
                continue;
            string name = fields[0];
            int n = data.Points[name].Length, test = Sparc.HoldoutCount(n);
            Check(holdoutNames.Add(name), "Unique saved holdout galaxy " + name);
            Check(n == Number(fields[4]) && n - test == Number(fields[5]) && test == Number(fields[6]), name + " held-out partition");
            comparedHoldouts++;
        }
        Check(comparedHoldouts == 131, "All 131 saved partitions checked");
        var report = new SparcCheckReport(failures.Count == 0, checks, failures.ToArray(), maxError, comparedPoints, savedMeta.Length, comparedHoldouts);
        string folder = Path.Combine(outputRoot, "sparc");
        Directory.CreateDirectory(folder);
        Sparc.SaveJson(Path.Combine(folder, "verification.json"), report);
        if (!report.Passed)
            throw new InvalidDataException("SPARC verification failed: " + string.Join("; ", failures.Take(10)));
        return report;
    }
}
