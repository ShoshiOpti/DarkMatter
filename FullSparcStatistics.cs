using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;
using MathNet.Numerics;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static class FullSparcStatistics
{
    public const int Draws = 99999;

    /// <summary>Finite-sample bounds for T=max(0,Q0-Q1). No radial independence or Wilks assumption.</summary>
    public static Result Bounds(int n, double q0, double q1)
    {
        if (n < 1 || !double.IsFinite(q0) || !double.IsFinite(q1) || q0 < 0 || q1 < 0)
            throw new ArgumentException("Positive count and finite nonnegative residual sums required.");
        double delta = q0 - q1, statistic = Math.Max(0, delta);
        double logP = statistic == 0 ? 0 : Math.Min(0, Math.Log(2.0 * n) + LogNormalSurvival(Math.Sqrt(statistic / n)));
        double varianceP = statistic == 0 ? 1 : Math.Min(1, n / statistic);
        return new()
        {
            ["rows"] = n, ["chi2_baryons"] = q0, ["chi2_green"] = q1, ["delta_chi2"] = delta,
            ["gaussian_marginal_conservative_p"] = Math.Exp(logP), ["gaussian_marginal_log10_p"] = logP / Math.Log(10),
            ["gaussian_marginal_equivalent_two_sided_Z"] = NormalZFromLogP(logP),
            ["variance_only_conservative_p"] = varianceP,
            ["variance_only_equivalent_two_sided_Z"] = Statistics.NormalMagnitude(varianceP),
            ["green_Q_per_row"] = q1 / n,
            ["interpretation"] = "Conservative p upper bounds / Z lower bounds; fixed grey null and correct quoted marginal errors. Not a detection of the action."
        };
    }

    public static double LogNormalSurvival(double z)
    {
        if (!double.IsFinite(z) || z < 0) throw new ArgumentOutOfRangeException(nameof(z));
        if (z < 26) return Math.Log(.5 * SpecialFunctions.Erfc(z / Math.Sqrt(2)));
        double term = 1, sum = 1, last = double.PositiveInfinity;
        for (int j = 1; j < 40; j++)
        {
            term *= -(2.0 * j - 1) / (z * z);
            if (Math.Abs(term) >= last) break;
            sum += term; last = Math.Abs(term);
            if (Math.Abs(term) < 1e-17 * Math.Abs(sum)) break;
        }
        return -.5 * z * z - Math.Log(z) - .5 * Math.Log(2 * Math.PI) + Math.Log(sum);
    }

    public static double NormalZFromLogP(double logP)
    {
        if (!double.IsFinite(logP) || logP > 0) throw new ArgumentOutOfRangeException(nameof(logP));
        if (logP == 0) return 0;
        double target = logP - Math.Log(2), lo = 0, hi = Math.Sqrt(-2 * target) + 2;
        for (int i = 0; i < 100; i++)
        {
            double mid = (lo + hi) / 2;
            if (LogNormalSurvival(mid) > target) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }

    /// <summary>Native conditional-binomial multinomial draws, verified against complete NumPy count streams.</summary>
    public static byte[] BootstrapCounts(JsonElement contract, int groups)
    {
        var expected = contract.GetProperty("streams").EnumerateArray().Single(s => s.GetProperty("groups").GetInt32() == groups);
        var random = new Pcg64Random(UInt128.Parse(contract.GetProperty("state").GetString()!, CultureInfo.InvariantCulture),
            UInt128.Parse(contract.GetProperty("increment").GetString()!, CultureInfo.InvariantCulture));
        double Next() => (random.Next64() >> 11) * (1.0 / 9007199254740992.0);
        int Binomial(int n, double p)
        {
            if (n == 0 || p == 0) return 0;
            bool reflect = p > .5; if (reflect) p = 1 - p;
            if (n * p > 30) throw new InvalidDataException("Frozen multinomial stream left its certified inversion regime.");
            double q = 1 - p, first = Math.Exp(n * double.LogP1(-p)), mass = first, u = Next();
            int bound = Math.Min(n, (int)(n * p + 10 * Math.Sqrt(n * p * q + 1))), x = 0;
            while (u > mass)
            {
                x++;
                if (x > bound) { x = 0; mass = first; u = Next(); }
                else { u -= mass; mass = (n - x + 1) * p * mass / (x * q); }
            }
            return reflect ? n - x : x;
        }
        var counts = new byte[Draws * groups]; double each = 1.0 / groups;
        for (int draw = 0; draw < Draws; draw++)
        {
            int remaining = groups; double probability = 1;
            for (int j = 0; j < groups - 1; j++)
            {
                int k = Binomial(remaining, each / probability); counts[draw * groups + j] = checked((byte)k);
                remaining -= k; if (remaining <= 0) break;
                probability -= each;
            }
            counts[draw * groups + groups - 1] = checked((byte)remaining);
        }
        string hash = Convert.ToHexStringLower(SHA256.HashData(counts));
        if (hash != expected.GetProperty("sha256_uint8_counts").GetString())
            throw new InvalidDataException($"Native multinomial stream mismatch for {groups} groups: {hash}");
        return counts;
    }

    public static Result Inference(double[] differences, string[] groups, byte[] counts)
    {
        if (differences.Length != groups.Length || differences.Length < 2 || differences.Any(v => !double.IsFinite(v)))
            throw new ArgumentException("Paired finite galaxy differences required.");
        string[] unique = groups.Distinct().Order(StringComparer.Ordinal).ToArray(); int ng = unique.Length;
        if (counts.Length != Draws * ng) throw new InvalidDataException("Wrong bootstrap stream dimensions.");
        double[] sum = new double[ng], number = new double[ng];
        for (int i = 0; i < groups.Length; i++) { int j = Array.BinarySearch(unique, groups[i], StringComparer.Ordinal); sum[j] += differences[i]; number[j]++; }
        double obs = differences.Average(); var sampled = new double[Draws]; int exceed = 0;
        for (int b = 0; b < Draws; b++)
        {
            double total = 0, n = 0; int offset = b * ng;
            for (int g = 0; g < ng; g++) { total += counts[offset + g] * sum[g]; n += counts[offset + g] * number[g]; }
            double v = sampled[b] = total / n;
            if (Math.Abs(v - obs) >= Math.Abs(obs)) exceed++;
        }
        int wins = 0, losses = 0;
        for (int g = 0; g < ng; g++) { if (sum[g] / number[g] < -1e-7) wins++; else if (sum[g] / number[g] > 1e-7) losses++; }
        double p = (exceed + 1.0) / (Draws + 1), signP = SignP(wins, losses);
        return new()
        {
            ["mean_difference"] = obs, ["bootstrap95"] = new[] { Data.Quantile(sampled, .025), Data.Quantile(sampled, .975) },
            ["two_sided_p"] = p, ["equivalent_Z"] = Statistics.NormalMagnitude(p), ["bootstrap_exceedances"] = exceed,
            ["replications"] = Draws, ["minimum_reportable_p"] = 1.0 / (Draws + 1), ["monte_carlo_floor"] = exceed == 0,
            ["groups"] = ng, ["group_wins"] = wins, ["group_losses"] = losses, ["group_ties"] = ng - wins - losses,
            ["group_sign_p"] = signP, ["group_sign_Z"] = Statistics.NormalMagnitude(signP),
            ["galaxy_wins"] = differences.Count(v => v < -1e-7), ["galaxy_losses"] = differences.Count(v => v > 1e-7)
        };
    }

    public static double SignP(int wins, int losses)
    {
        if (wins < 0 || losses < 0) throw new ArgumentOutOfRangeException(nameof(wins));
        int n = wins + losses; if (n == 0) return 1;
        double term = Math.Pow(.5, n), sum = term;
        for (int k = 1; k <= Math.Min(wins, losses); k++) { term *= (n - k + 1.0) / k; sum += term; }
        return Math.Min(1, 2 * sum);
    }

    public static void AdjustHolm(List<Result> tests)
    {
        foreach (var (raw, adjusted, z) in new[] { ("two_sided_p", "Holm_p", "Holm_Z"), ("group_sign_p", "Holm_group_sign_p", "Holm_group_sign_Z") })
        {
            var order = tests.OrderBy(t => (double)t[raw]!).ToArray(); double current = 0;
            for (int i = 0; i < order.Length; i++)
            {
                current = Math.Min(1, Math.Max(current, (order.Length - i) * (double)order[i][raw]!));
                order[i][adjusted] = current; order[i][z] = Statistics.NormalMagnitude(current);
            }
        }
    }
}
