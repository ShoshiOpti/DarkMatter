using System.Globalization;
using System.Text.Json;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class PublicationStatistics
{
    static Result[] Bootstrap(double[][] contrasts, int[] counts, ulong pairsSeed, ulong? wildSeed, JsonElement reference)
    {
        int K = contrasts.Length, G = counts.Length;
        Require(K > 0 && G >= 2 && contrasts.All(d => d.Length == G && d.All(double.IsFinite)), "Invalid bootstrap contrasts.");
        double[] mu = contrasts.Select(d => d.Average()).ToArray();
        double[] se = contrasts.Select((d, j) => Math.Sqrt(d.Sum(x => Square(x - mu[j])) / (G - 1) / G)).ToArray();
        bool[] zero = contrasts.Select(d => d.All(x => x == 0)).ToArray();
        Require(Enumerable.Range(0, K).All(j => se[j] > 0 || zero[j]), "Nonzero constant effect makes studentization undefined.");
        double[][] means = Enumerable.Range(0, K).Select(_ => new double[Draws]).ToArray();
        double[][] tstars = Enumerable.Range(0, K).Select(_ => new double[Draws]).ToArray();
        double[][] pooled = Enumerable.Range(0, K).Select(_ => new double[Draws]).ToArray();
        var rng = new PublicationStatisticsRandom(reference, pairsSeed);
        int[] ids = new int[G];
        double[] bm = new double[K], variances = new double[K], weighted = new double[K];
        for (int draw = 0; draw < Draws; draw++)
        {
            Array.Clear(bm); Array.Clear(variances); Array.Clear(weighted);
            int n = 0;
            for (int g = 0; g < G; g++)
            {
                int id = ids[g] = rng.NextInt(G);
                n += counts[id];
                for (int j = 0; j < K; j++) { double x = contrasts[j][id]; bm[j] += x; weighted[j] += x * counts[id]; }
            }
            for (int j = 0; j < K; j++) bm[j] /= G;
            for (int g = 0; g < G; g++)
            for (int j = 0; j < K; j++) variances[j] += Square(contrasts[j][ids[g]] - bm[j]);
            for (int j = 0; j < K; j++)
            {
                double bs = Math.Sqrt(variances[j] / (G - 1) / G);
                Require(bs > 0 || zero[j], "Degenerate nonconstant pairs bootstrap draw; draws must not be discarded.");
                means[j][draw] = bm[j]; tstars[j][draw] = bs > 0 ? (bm[j] - mu[j]) / bs : 0;
                pooled[j][draw] = weighted[j] / n;
            }
        }
        int[] hits = new int[K];
        double[] observed = Enumerable.Range(0, K).Select(j => se[j] > 0 ? mu[j] / se[j] : 0).ToArray();
        if (wildSeed is not null)
        {
            rng = new PublicationStatisticsRandom(reference, wildSeed.Value);
            double[] sumSquares = contrasts.Select(d => d.Sum(Square)).ToArray();
            for (int start = 0; start < Draws; start += 5000)
            {
                rng.StartInt8Batch();
                for (int draw = start; draw < Math.Min(start + 5000, Draws); draw++)
                {
                    Array.Clear(bm);
                    for (int g = 0; g < G; g++)
                    {
                        int sign = 2 * rng.NextSignBit() - 1;
                        for (int j = 0; j < K; j++) bm[j] += sign * contrasts[j][g];
                    }
                    for (int j = 0; j < K; j++)
                    {
                        double mean = bm[j] / G;
                        double ss = Math.Sqrt(Math.Max((sumSquares[j] - G * mean * mean) / (G - 1), 0) / G);
                        Require(ss > 0 || zero[j], "Nonzero degenerate wild bootstrap draw.");
                        double t = ss > 0 ? mean / ss : 0;
                        if (Math.Abs(t) >= Math.Abs(observed[j])) hits[j]++;
                    }
                }
            }
        }
        var results = new Result[K];
        for (int j = 0; j < K; j++)
        {
            double[] d = contrasts[j], sorted = d.Order().ToArray();
            double sum = d.Sum(), denom = d.Sum(x => Square(x - mu[j]));
            double[] shares = d.Select(x => denom > 0 ? Square(x - mu[j]) / denom : 0).ToArray();
            double[] loo = d.Select(x => (sum - x) / (G - 1)).ToArray();
            int trim = (int)Math.Floor(.1 * G);
            Array.Sort(means[j]); Array.Sort(tstars[j]); Array.Sort(pooled[j]);
            var row = new Result
            {
                ["effect"] = mu[j], ["cluster_se"] = se[j],
                ["bootstrap_t_ci_low"] = mu[j] - QuantileSorted(tstars[j], .975) * se[j],
                ["bootstrap_t_ci_high"] = mu[j] - QuantileSorted(tstars[j], .025) * se[j],
                ["percentile_ci_low"] = QuantileSorted(means[j], .025), ["percentile_ci_high"] = QuantileSorted(means[j], .975),
                ["median_paired_difference"] = QuantileSorted(sorted, .5), ["trimmed10_paired_difference"] = sorted.Skip(trim).Take(G - 2 * trim).Average(),
                ["pooled_point_effect"] = d.Select((value, g) => value * counts[g]).Sum() / counts.Sum(),
                ["pooled_percentile_ci_low"] = QuantileSorted(pooled[j], .025), ["pooled_percentile_ci_high"] = QuantileSorted(pooled[j], .975),
                ["leave_one_out_min"] = loo.Min(), ["leave_one_out_max"] = loo.Max(),
                ["leave_one_out_reversals"] = loo.Count(value => Math.Sign(value) != Math.Sign(mu[j])),
                ["max_variance_share"] = shares.Max(), ["top3_variance_share"] = shares.OrderDescending().Take(3).Sum(),
                ["variance_concentration_count"] = denom > 0 ? 1 / shares.Sum(Square) : null, ["exact_identical_scores"] = zero[j],
                ["wild_exceedances"] = null, ["wild_draws"] = null, ["p_wild_raw"] = null, ["t_statistic"] = null,
                ["monte_carlo_floor"] = null, ["raw_MC_95low"] = null, ["raw_MC_95high"] = null
            };
            if (wildSeed is not null)
            {
                int k = hits[j];
                row["wild_exceedances"] = k; row["wild_draws"] = Draws; row["p_wild_raw"] = (k + 1.0) / (Draws + 1);
                row["t_statistic"] = observed[j]; row["monte_carlo_floor"] = k == 0;
                row["raw_MC_95low"] = k == 0 ? 0 : Statistics.BetaInverse(.025, k, Draws - k + 1);
                row["raw_MC_95high"] = k == Draws ? 1 : Statistics.BetaInverse(.975, k + 1, Draws - k);
            }
            results[j] = row;
        }
        return results;
    }

    static double QuantileSorted(double[] values, double probability)
    {
        double index = probability * (values.Length - 1);
        int left = (int)index;
        return values[left] + (index - left) * (values[Math.Min(left + 1, values.Length - 1)] - values[left]);
    }

    static Result CompareReference(string path, List<Result> actual, string[] keys, string[]? fields = null)
    {
        var expected = Csv.Read(path);
        Require(expected.Count == actual.Count, $"Reference row count differs for {Path.GetFileName(path)}.");
        string KeyValue(object? value) => value is double d ? d.ToString("G17", CultureInfo.InvariantCulture) : Convert.ToString(value, CultureInfo.InvariantCulture) ?? "";
        string ExpectedKey(Dictionary<string, string> r) => string.Join('|', keys.Select(k =>
            double.TryParse(r[k], NumberStyles.Float, CultureInfo.InvariantCulture, out double d) ? d.ToString("G17", CultureInfo.InvariantCulture) : r[k]));
        var byKey = actual.ToDictionary(r => string.Join('|', keys.Select(k => KeyValue(r[k]))));
        int checkedValues = 0;
        double maximumAbsolute = 0, maximumRelative = 0;
        string largestField = "";
        foreach (var row in expected)
        {
            string key = ExpectedKey(row);
            Require(byKey.TryGetValue(key, out var current), $"Reference row absent: {key}.");
            foreach (string field in fields ?? row.Keys.ToArray())
            {
                Require(current!.ContainsKey(field), $"Recomputed table lacks reference column {field}.");
                string value = row[field];
                object? got = current[field];
                if (value.Length == 0) Require(got is null || Convert.ToString(got, CultureInfo.InvariantCulture) == "", $"Expected absent value: {key}/{field}.");
                else if (bool.TryParse(value, out bool flag)) Require(got is bool b && b == flag, $"Boolean mismatch: {key}/{field}.");
                else if (double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double number))
                {
                    Require(got is not null, $"Missing numeric value: {key}/{field}.");
                    double numberGot = Convert.ToDouble(got, CultureInfo.InvariantCulture);
                    bool integer = field is "wins" or "losses" or "ties" or "joint_failures" or "informative_n" or "informative_signs"
                        or "wild_exceedances" or "wild_draws" or "n_galaxies" or "n_points" or "pairs_seed_used" or "wild_seed_used" or "leave_one_out_reversals"
                        || field.EndsWith("_galaxies") || field.EndsWith("_points");
                    bool probability = field.EndsWith("_p") || field.StartsWith("p_");
                    double absoluteTolerance = integer ? 0 : probability ? 1e-45 : 2e-10;

                    double relativeTolerance = integer ? 0 : 2e-10;
                    Require(Near(numberGot, number, absoluteTolerance, relativeTolerance),
                        $"Publication reference mismatch {Path.GetFileName(path)} {key}/{field}: actual={numberGot:R}, expected={number:R}.");
                    double error = Math.Abs(numberGot - number);
                    if (error > maximumAbsolute) { maximumAbsolute = error; largestField = key + "/" + field; }
                    maximumRelative = Math.Max(maximumRelative, number == 0 ? 0 : error / Math.Abs(number));
                }
                else Require(Convert.ToString(got, CultureInfo.InvariantCulture) == value, $"Text mismatch: {key}/{field}.");
                checkedValues++;
            }
        }
        return new() { ["status"] = "pass", ["reference_file"] = Path.GetFileName(path), ["rows"] = expected.Count,
            ["checked_values"] = checkedValues, ["maximum_absolute_error"] = maximumAbsolute,
            ["maximum_relative_error_nonzero"] = maximumRelative, ["largest_absolute_error_field"] = largestField,
            ["tolerance"] = "Counts and flags exact; p values relative 2e-10 (absolute floor 1e-45); other values absolute 2e-10 plus relative 2e-10." };
    }
}


