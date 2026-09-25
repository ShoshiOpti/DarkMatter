using System.IO.Compression;
using System.Text;
using System.Text.Json;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class PointwiseStatistics
{
    static Result Validate(string source, string destination, Result profile, Result checks)
    {
        var tables = new Result();
        string[] names = ["model_galaxy_squared_losses.csv", "pointwise_paired_differences.csv", "score_sensitivity.csv", "primary_cluster_results.csv", "galaxy_influence.csv", "galaxy_contrasts.csv", "pooled_point_sensitivity.csv", "residual_distribution_summary.csv"];
        foreach (string name in names)
            tables[name] = Statistics.ValidateCsv(Path.Combine(source, "results", name), Path.Combine(destination, name));
        int fields = 0;
        foreach (string name in new[] { "profiled_baseline_failure_sensitivity.json", "checks.json" })
            CompareJson(Data.Json(Path.Combine(source, "results", name)), Data.Json(Path.Combine(destination, name)), name, ref fields);
        var primary = Csv.Read(Path.Combine(destination, "primary_cluster_results.csv"));
        Require(primary.Select(r => (int)Csv.Number(r, "wild_exceedances")).SequenceEqual(new[] { 0, 0, 7145, 0 }), "Wild-bootstrap tail counts differ from the frozen run.");
        Require(primary.Count(r => Csv.Flag(r, "monte_carlo_floor")) == 3, "Simulation-floor comparisons must be flagged.");
        Require(primary.All(r => Csv.Number(r, "p_wild") > 0 && Csv.Number(r, "p_holm4") >= Csv.Number(r, "p_wild")), "Invalid simulation probability or multiplicity adjustment.");
        Require(profile["full_sample_p"] is null && !(bool)profile["primary_p_assigned"]!, "A failed full-sample profile must not receive a finite primary test.");
        return new Result
        {
            ["passed"] = true,
            ["reference"] = "numerics/" + Analysis + "/results",
            ["csv_tables"] = tables,
            ["json_fields_checked"] = fields,
            ["input_sha256"] = checks["input_sha256"],
            ["numeric_tolerance"] = "1e-10 absolute plus 1e-10 relative; p_* fields use 1e-9 relative with 1e-300 absolute",
            ["primary_exceedances"] = new[] { 0, 0, 7145, 0 },
            ["simulation_floor_comparisons"] = 3,
            ["numpy_random_streams"] = PointwiseRandom.SelfCheck(Path.Combine(source, "random_streams.json")),
            ["bootstrap_archive"] = ValidateArchive(Path.Combine(destination, "bootstrap_summaries.npz"))
        };
    }

    static void CompareJson(JsonElement expected, JsonElement actual, string path, ref int fields)
    {
        Require(expected.ValueKind == actual.ValueKind, "JSON value type mismatch: " + path);
        switch (expected.ValueKind)
        {
            case JsonValueKind.Object:
                Require(expected.EnumerateObject().Count() == actual.EnumerateObject().Count(), "JSON field count mismatch: " + path);
                foreach (var p in expected.EnumerateObject())
                    CompareJson(p.Value, actual.GetProperty(p.Name), path + "." + p.Name, ref fields);
                break;
            case JsonValueKind.Array:
                Require(expected.GetArrayLength() == actual.GetArrayLength(), "JSON array length mismatch: " + path);
                for (int i = 0; i < expected.GetArrayLength(); i++)
                    CompareJson(expected[i], actual[i], path + "[" + i + "]", ref fields);
                break;
            case JsonValueKind.Number:
                double a = expected.GetDouble(), b = actual.GetDouble();
                Require(Math.Abs(a - b) <= 1e-10 + 1e-10 * Math.Abs(a), $"JSON numeric mismatch: {path}, expected {a:R}, actual {b:R}.");
                fields++;
                break;
            default:
                Require(expected.GetRawText() == actual.GetRawText() || expected.ToString() == actual.ToString(), "JSON value mismatch: " + path);
                fields++;
                break;
        }
    }

    static Result ValidateArchive(string path)
    {
        using var archive = ZipFile.OpenRead(path);
        Require(archive.Entries.Count == 3, "Bootstrap NPZ array count.");
        long numbers = 0;
        foreach (var (name, columns) in new[] { ("means", 24), ("tstats", 24), ("pooled_means", 4) })
        {
            var entry = archive.GetEntry(name + ".npy") ?? throw new InvalidDataException("Missing bootstrap array: " + name);
            using var stream = entry.Open();
            using var reader = new BinaryReader(stream, Encoding.ASCII);
            Require(reader.ReadBytes(8).SequenceEqual(new byte[] { 0x93, (byte)'N', (byte)'U', (byte)'M', (byte)'P', (byte)'Y', 1, 0 }), "Invalid NPY signature.");
            int length = reader.ReadUInt16();
            string header = Encoding.ASCII.GetString(reader.ReadBytes(length));
            Require(header.Contains($"'shape': ({Draws}, {columns})") && header.Contains("'descr': '<f8'") && header.Contains("'fortran_order': False"), "Invalid NPY shape/type/order.");
            Require(entry.Length == 10 + length + (long)Draws * columns * sizeof(double), "Invalid bootstrap array byte length.");
            for (int b = 0; b < Draws; b++)
                for (int j = 0; j < columns; j++)
                {
                    Require(double.IsFinite(reader.ReadDouble()), "Non-finite retained bootstrap draw.");
                    numbers++;
                }
        }
        return new Result
        {
            ["passed"] = true,
            ["arrays"] = 3,
            ["finite_values"] = numbers,
            ["format"] = "NumPy NPZ, little-endian float64, C row order"
        };
    }

    static Result SelfChecks()
    {
        int count = 0;
        void Check(bool condition, string name)
        {
            Require(condition, "self-check: " + name);
            count++;
        }
        var (mean, se) = MeanSe([1, 2, 3]);
        Check(mean == 2 && Math.Abs(se - 1 / Math.Sqrt(3)) < 1e-15, "sample standard error uses G-1");
        Check(MeanSe([1, 1, 1]).Se == 0, "degenerate sample is detectable");
        double[] q = QuantilePair([0, 10, 20, 30, 40]);
        Check(Math.Abs(q[0] - 1) < 1e-14 && Math.Abs(q[1] - 39) < 1e-14, "NumPy linear quantile interpolation");
        Check(Correlated(3, 2, 5, 0) == 7, "rho zero gives mean squared residual");
        Check(Correlated(0, 2, 5, .5) == 4.0 / 3, "correlated score downweights common residual mode");
        Check(Correlated(3, 0, 5, .5) == 6, "correlated score upweights centered residual mode");
        Check(1.0 / (Draws + 1) == .00001 && Statistics.Holm([.00001, .00001, .07146, .00001]).SequenceEqual(new[] { .00004, .00004, .07146, .00004 }), "finite-resolution p and four-comparison Holm family");
        Check(Math.Abs(Statistics.BetaInverse(.975, 1, Draws) - (1 - Math.Pow(.025, 1.0 / Draws))) < 1e-12, "zero-exceedance Clopper-Pearson boundary");
        Check(Math.Abs(Statistics.RegularizedBeta(.5, .5, .5) - .5) < 1e-13, "Student-t Cauchy reference");
        return new Result { ["passed"] = true, ["checks"] = count };
    }
}
