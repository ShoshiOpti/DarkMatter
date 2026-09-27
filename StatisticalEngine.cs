using System.Numerics;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace DarkUniverse;

/// <summary>Versioned settings for optional calibration and mechanism-resolution studies.</summary>
public sealed record AnalysisStudyProtocol
{
    public int SimulationReplicates { get; init; } = 500;
    public int BootstrapDraws { get; init; } = 999;
    public int GalaxyCount { get; init; } = 131;
    public ulong Seed { get; init; } = 2026092701;
    public double ConfidenceLevel { get; init; } = .95;
    public double TestAlpha { get; init; } = .05;
    public double StandardizedEffect { get; init; } = .3;
    public int MaxParallelism { get; init; } = 4;
    public double? MechanismSpeedToleranceKms { get; init; }
    public double? MechanismLossMargin { get; init; }
    public double? SolverSpeedErrorBoundKms { get; init; }

    public void Validate()
    {
        if (SimulationReplicates < 1 || BootstrapDraws < 99 || GalaxyCount < 8 || MaxParallelism < 1)
            throw new ArgumentException("Research studies require positive simulation count/parallelism, at least 99 bootstrap draws and at least eight galaxies.");
        if (!(ConfidenceLevel > 0 && ConfidenceLevel < 1) || !(TestAlpha > 0 && TestAlpha < 1) ||
            !double.IsFinite(StandardizedEffect) || StandardizedEffect == 0)
            throw new ArgumentException("Confidence level and test alpha must lie between zero and one; the known effect must be finite and nonzero.");
        foreach (double? margin in new[] { MechanismSpeedToleranceKms, MechanismLossMargin })
            if (margin is not null && (!double.IsFinite(margin.Value) || margin.Value <= 0))
                throw new ArgumentException("An explicitly supplied mechanism margin must be finite and positive.");
        if (SolverSpeedErrorBoundKms is double error && (!double.IsFinite(error) || error < 0))
            throw new ArgumentException("The supplied per-prediction solver speed error bound must be finite and nonnegative.");
    }
}

public sealed record PairedGalaxyObservation(string Galaxy, int PointCount, double Difference);
public sealed record StatisticalInterval(double Lower, double Upper);
public sealed record PairedBootstrapOptions
{
    public int Draws { get; init; } = 999;
    public ulong PairsSeed { get; init; } = 2026092701;
    public ulong WildSeed { get; init; } = 2026092702;
    public double ConfidenceLevel { get; init; } = .95;
}

public sealed record PairedBootstrapResult
{
    public string Status { get; init; } = "";
    public int Galaxies { get; init; }
    public long Points { get; init; }
    public double Effect { get; init; }
    public double ClusterStandardError { get; init; }
    public double PooledPointEffect { get; init; }
    public double MedianEffect { get; init; }
    public double Trimmed10Effect { get; init; }
    public double? ObservedT { get; init; }
    public StatisticalInterval? BootstrapTInterval { get; init; }
    public StatisticalInterval? PercentileInterval { get; init; }
    public StatisticalInterval? PooledPercentileInterval { get; init; }
    public double? WildPValue { get; init; }
    public int? WildExceedances { get; init; }
    public StatisticalInterval? WildMonteCarloInterval { get; init; }
    public bool MonteCarloFloor { get; init; }
    public int Draws { get; init; }
    public int DegeneratePairsDraws { get; init; }
    public int DegenerateWildDraws { get; init; }
    public double MaximumVarianceShare { get; init; }
    public double? VarianceConcentrationCount { get; init; }
    public ulong PairsSeed { get; init; }
    public ulong WildSeed { get; init; }
    public double ConfidenceLevel { get; init; }
    public string RandomAlgorithm { get; init; } = Pcg64Random.ResearchAlgorithm;
    public string Scope { get; init; } = "Conditional paired-galaxy mean inference assuming independent clusters. Marginal pairs-bootstrap-t intervals do not invert the approximate null-imposed Rademacher wild test.";
}

/// <summary>
/// Reusable research engine. The frozen-publication adapters retain their original
/// reduction order, archived protocols and NumPy seed states for exact reproduction.
/// </summary>
public static class StatisticalEngine
{
    public static PairedBootstrapResult Run(IReadOnlyList<PairedGalaxyObservation> observations, PairedBootstrapOptions? options = null)
    {
        options ??= new();
        if (observations.Count < 2 || observations.Select(o => o.Galaxy).Distinct(StringComparer.Ordinal).Count() != observations.Count ||
            observations.Any(o => string.IsNullOrWhiteSpace(o.Galaxy) || o.PointCount <= 0 || !double.IsFinite(o.Difference)))
            throw new ArgumentException("Supply at least two unique galaxy identifiers, positive point counts, and finite paired differences; unresolved failures must be handled before inference.");
        if (options.Draws < 99 || !(options.ConfidenceLevel > 0 && options.ConfidenceLevel < 1))
            throw new ArgumentException("At least 99 draws and a confidence level strictly between zero and one are required.");

        // Sorting makes a seed-defined experiment invariant to input row order.
        var rows = observations.OrderBy(o => o.Galaxy, StringComparer.Ordinal).ToArray();
        int G = rows.Length, B = options.Draws;
        double[] raw = rows.Select(o => o.Difference).ToArray();
        double scale = raw.Max(Math.Abs);
        double[] x = scale > 0 ? raw.Select(v => v / scale).ToArray() : new double[G];
        double mean = Sum(x) / G, variance = Sum(x.Select(v => Square(v - mean))) / (G - 1), se = Math.Sqrt(variance / G);
        long points = rows.Sum(o => (long)o.PointCount);
        double pooledEffect = Sum(rows.Select((o, i) => x[i] * o.PointCount)) / points * scale;
        double[] sorted = raw.Order().ToArray();
        int trim = (int)Math.Floor(.1 * G);
        double sumDeviations = variance * (G - 1);
        double[] shares = x.Select(v => sumDeviations > 0 ? Square(v - mean) / sumDeviations : 0).ToArray();
        var result = new PairedBootstrapResult
        {
            Status = "pass", Galaxies = G, Points = points, Effect = mean * scale, ClusterStandardError = se * scale,
            PooledPointEffect = pooledEffect, MedianEffect = Quantile(sorted, .5),
            Trimmed10Effect = Sum(sorted.Skip(trim).Take(G - 2 * trim).Select(v => scale > 0 ? v / scale : 0)) / (G - 2 * trim) * scale,
            MaximumVarianceShare = shares.Max(), VarianceConcentrationCount = sumDeviations > 0 ? 1 / Sum(shares.Select(Square)) : null,
            Draws = B, PairsSeed = options.PairsSeed, WildSeed = options.WildSeed, ConfidenceLevel = options.ConfidenceLevel
        };
        if (scale == 0)
            return result with { Status = "identical_zero", ObservedT = 0, BootstrapTInterval = new(0, 0),
                PercentileInterval = new(0, 0), PooledPercentileInterval = new(0, 0), WildPValue = 1,
                WildExceedances = B, WildMonteCarloInterval = BinomialInterval(B, B) };
        if (se == 0 || !double.IsFinite(se))
            return result with { Status = "undefined_nonzero_constant_effect" };

        var pairsRandom = Pcg64Random.FromSeed(options.PairsSeed);
        var wildRandom = Pcg64Random.FromSeed(options.WildSeed);
        var means = new double[B]; var tstars = new double[B]; var pooled = new double[B];
        int[] ids = new int[G]; double[] signed = new double[G];
        int degeneratePairs = 0, degenerateWild = 0, hits = 0;
        double observed = mean / se;
        for (int draw = 0; draw < B; draw++)
        {
            double bm = 0, weighted = 0; long denominator = 0;
            for (int g = 0; g < G; g++)
            {
                int id = ids[g] = pairsRandom.NextInt(G);
                bm += x[id]; weighted += x[id] * rows[id].PointCount; denominator += rows[id].PointCount;
            }
            bm /= G;
            double vs = 0;
            for (int g = 0; g < G; g++) vs += Square(x[ids[g]] - bm);
            double bs = Math.Sqrt(vs / (G - 1) / G);
            means[draw] = bm * scale; pooled[draw] = weighted / denominator * scale;
            if (!(bs > 0) || !double.IsFinite(bs)) degeneratePairs++;
            else tstars[draw] = (bm - mean) / bs;

            // Direct centered moments avoid Q-G*mean^2 cancellation in nearly
            // constant sign draws. Research streams use one uninterrupted int32 stream.
            double wm = 0;
            for (int g = 0; g < G; g++) { signed[g] = (2 * wildRandom.NextInt(2) - 1) * x[g]; wm += signed[g]; }
            wm /= G; double wv = 0;
            for (int g = 0; g < G; g++) wv += Square(signed[g] - wm);
            double ws = Math.Sqrt(wv / (G - 1) / G);
            if (!(ws > 0) || !double.IsFinite(ws)) degenerateWild++;
            else if (Math.Abs(wm / ws) >= Math.Abs(observed)) hits++;
        }
        double alpha = (1 - options.ConfidenceLevel) / 2;
        Array.Sort(means); Array.Sort(tstars); Array.Sort(pooled);
        StatisticalInterval? interval = degeneratePairs == 0
            ? new((mean - Quantile(tstars, 1 - alpha) * se) * scale, (mean - Quantile(tstars, alpha) * se) * scale) : null;
        return result with
        {
            Status = degeneratePairs + degenerateWild == 0 ? "pass" : "undefined_degenerate_resampling",
            ObservedT = observed, BootstrapTInterval = interval,
            PercentileInterval = new(Quantile(means, alpha), Quantile(means, 1 - alpha)),
            PooledPercentileInterval = new(Quantile(pooled, alpha), Quantile(pooled, 1 - alpha)),
            WildPValue = degenerateWild == 0 ? (hits + 1.0) / (B + 1) : null,
            WildExceedances = degenerateWild == 0 ? hits : null,
            WildMonteCarloInterval = degenerateWild == 0 ? BinomialInterval(hits, B) : null,
            MonteCarloFloor = degenerateWild == 0 && hits == 0,
            DegeneratePairsDraws = degeneratePairs, DegenerateWildDraws = degenerateWild
        };
    }

    public static StatisticalInterval BinomialInterval(int successes, int trials, double confidenceLevel = .95)
    {
        if (trials < 1 || successes < 0 || successes > trials || !(confidenceLevel > 0 && confidenceLevel < 1))
            throw new ArgumentException("Invalid binomial interval arguments.");
        double tail = (1 - confidenceLevel) / 2;
        return new(successes == 0 ? 0 : Statistics.BetaInverse(tail, successes, trials - successes + 1),
            successes == trials ? 1 : Statistics.BetaInverse(1 - tail, successes + 1, trials - successes));
    }

    public static ulong DeriveSeed(ulong seed, ulong stream) => Pcg64Random.MixSeed(unchecked(seed + stream * 0x9e3779b97f4a7c15UL));
    internal static double Square(double value) => value * value;
    internal static double Sum(IEnumerable<double> values)
    {
        double sum = 0, correction = 0;
        foreach (double value in values) { double nextValue = value - correction, next = sum + nextValue; correction = (next - sum) - nextValue; sum = next; }
        return sum;
    }
    internal static double Quantile(double[] sorted, double level)
    {
        double position = level * (sorted.Length - 1); int left = (int)position;
        // Convex interpolation avoids overflow for finite values of opposite signs.
        double fraction = position - left;
        return sorted[left] * (1 - fraction) + sorted[Math.Min(left + 1, sorted.Length - 1)] * fraction;
    }

    public static Dictionary<string, object?> SelfCheck()
    {
        int checks = 0;
        void Check(bool condition, string label) { if (!condition) throw new InvalidDataException("Statistical engine self-check: " + label); checks++; }
        var frozen = new Pcg64Random(UInt128.Parse("114815172112769181402702086940346429578"), UInt128.Parse("147885377315751686733782345246761456639"));
        ulong[] known = [12614079201176063507UL, 8874092783013893688UL, 189220602283836686UL, 338392960798599543UL];
        Check(known.All(value => frozen.Next64() == value), "frozen NumPy PCG64 raw golden stream");
        var first = Pcg64Random.FromSeed(12345); var second = Pcg64Random.FromSeed(12345); var other = Pcg64Random.FromSeed(12346);
        Check(Enumerable.Range(0, 32).All(_ => first.Next64() == second.Next64()), "repeat seed determinism");
        Check(Pcg64Random.FromSeed(12345).Next64() != other.Next64(), "distinct seed");
        var researchGolden = Pcg64Random.FromSeed(12345);
        ulong[] researchKnown = [4268611720711217064UL, 14502665754970336712UL, 7585861177730345277UL, 3939900106928982441UL];
        Check(researchKnown.All(value => researchGolden.Next64() == value), "research seed expansion matches independent Python integer reference");
        var bounds = Pcg64Random.FromSeed(12);
        Check(Enumerable.Range(0, 10000).All(_ => bounds.NextInt(131) is >= 0 and < 131), "bounded integers");
        var input = Enumerable.Range(0, 64).Select(i => new PairedGalaxyObservation($"g{i:D3}", 2 + i % 7, Math.Sin(i * 1.3) + (i % 5 - 2) * .1)).ToArray();
        var options = new PairedBootstrapOptions { Draws = 199, PairsSeed = 100, WildSeed = 200 };
        var analytical = Run([new("a", 1, 1), new("b", 1, 2), new("c", 1, 3)], options);
        Check(Math.Abs(analytical.Effect - 2) < 1e-14 && Math.Abs(analytical.ClusterStandardError - 1 / Math.Sqrt(3)) < 1e-14, "analytic sample mean and G-1 standard error");
        Check(analytical.DegeneratePairsDraws > 0 && analytical.BootstrapTInterval is null && analytical.Status == "undefined_degenerate_resampling", "degenerate pairs draws invalidate the interval instead of being dropped");
        var baseline = Run(input, options);
        var rescaled = Run(input.Select(o => o with { Difference = 7 * o.Difference }).ToArray(), options);
        var shifted = Run(input.Select(o => o with { Difference = o.Difference + 2 }).ToArray(), options);
        bool Close(double a, double b) => Math.Abs(a - b) <= 2e-12 * (1 + Math.Abs(b));
        Check(Close(baseline.Effect, Sum(input.Select(o => o.Difference)) / input.Length), "known equal-galaxy mean");
        Check(rescaled.WildExceedances == baseline.WildExceedances && Close(rescaled.BootstrapTInterval!.Lower, 7 * baseline.BootstrapTInterval!.Lower), "positive-scale invariant test and interval");
        Check(Close(shifted.BootstrapTInterval!.Lower, baseline.BootstrapTInterval!.Lower + 2) && Close(shifted.BootstrapTInterval.Upper, baseline.BootstrapTInterval.Upper + 2), "translation equivariance of pairs interval");
        Check(Run(input.Reverse().ToArray(), options) == baseline, "row-order independent result");
        var zeros = Run(input.Select(o => o with { Difference = 0 }).ToArray(), options);
        Check(zeros.Status == "identical_zero" && zeros.WildPValue == 1 && zeros.BootstrapTInterval == new StatisticalInterval(0, 0), "exact zero control");
        var constant = Run(input.Select(o => o with { Difference = 2 }).ToArray(), options);
        Check(constant.Status == "undefined_nonzero_constant_effect" && constant.WildPValue is null && constant.BootstrapTInterval is null, "undefined constant effect is explicit");
        var tiny = Run(input.Select(o => o with { Difference = o.Difference * 1e-200 }).ToArray(), options);
        Check(tiny.Status == "pass" && tiny.WildExceedances == baseline.WildExceedances, "scaled moments avoid tiny-value underflow");
        var huge = Run(input.Select(o => o with { Difference = o.Difference * 1e200 }).ToArray(), options);
        Check(huge.Status == "pass" && huge.WildExceedances == baseline.WildExceedances, "scaled moments avoid large-value overflow");
        bool invalid = false;
        try { _ = Run([new("duplicate", 1, 1), new("duplicate", 1, 2)], options); } catch (ArgumentException) { invalid = true; }
        Check(invalid, "duplicate galaxy rejected");
        Check(BinomialInterval(0, 100).Lower == 0 && BinomialInterval(100, 100).Upper == 1, "binomial boundary coverage");
        Check(Math.Abs(BinomialInterval(0, 100).Upper - (1 - Math.Pow(.025, .01))) < 1e-12, "analytic zero-tail Clopper-Pearson endpoint");
        return new() { ["status"] = "pass", ["checks"] = checks, ["algorithm"] = Pcg64Random.ResearchAlgorithm };
    }
}

/// <summary>Shared PCG64 XSL-RR core. Explicit-state construction preserves NumPy frozen streams.</summary>
internal sealed class Pcg64Random
{
    internal const string ResearchAlgorithm = "PCG64-XSL-RR with SplitMix64 seed expansion v1; research seed mapping is explicitly distinct from NumPy SeedSequence";
    static readonly UInt128 Multiplier = ((UInt128)2549297995355413924UL << 64) | 4865540595714422341UL;
    UInt128 state;
    readonly UInt128 increment;
    uint cached32, byteBuffer;
    bool hasCached32;
    int bytesRemaining;
    bool hasNormal;
    double cachedNormal;

    internal Pcg64Random(UInt128 initialState, UInt128 streamIncrement)
    {
        if ((streamIncrement & 1) == 0) throw new ArgumentException("PCG64 stream increment must be odd.");
        state = initialState; increment = streamIncrement;
    }
    internal static Pcg64Random FromSeed(ulong seed)
    {
        ulong a = MixSeed(seed), b = MixSeed(unchecked(seed + 0x9e3779b97f4a7c15UL));
        ulong c = MixSeed(unchecked(seed + 2 * 0x9e3779b97f4a7c15UL)), d = MixSeed(unchecked(seed + 3 * 0x9e3779b97f4a7c15UL));
        UInt128 initialState = ((UInt128)a << 64) | b;
        UInt128 stream = unchecked(((((UInt128)c << 64) | d) << 1) | 1);
        var rng = new Pcg64Random(0, stream);
        _ = rng.Next64(); rng.state = unchecked(rng.state + initialState); _ = rng.Next64();
        return rng;
    }
    internal static ulong MixSeed(ulong value)
    {
        value = unchecked(value + 0x9e3779b97f4a7c15UL);
        value = unchecked((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9UL);
        value = unchecked((value ^ (value >> 27)) * 0x94d049bb133111ebUL);
        return value ^ (value >> 31);
    }
    internal ulong Next64()
    {
        state = unchecked(state * Multiplier + increment);
        return BitOperations.RotateRight((ulong)(state >> 64) ^ (ulong)state, (int)(state >> 122));
    }
    uint Next32()
    {
        if (hasCached32) { hasCached32 = false; return cached32; }
        ulong value = Next64(); cached32 = (uint)(value >> 32); hasCached32 = true; return (uint)value;
    }
    internal int NextInt(int upper)
    {
        if (upper <= 0) throw new ArgumentOutOfRangeException(nameof(upper));
        if (upper == 1) return 0;
        uint bound = (uint)upper, threshold = unchecked(0u - bound) % bound; ulong product;
        do { product = (ulong)Next32() * bound; } while ((uint)product < threshold);
        return (int)(product >> 32);
    }
    internal void StartInt8Batch() => bytesRemaining = 0;
    internal int NextInt8(int upper)
    {
        if (upper <= 0 || upper > 128) throw new ArgumentOutOfRangeException(nameof(upper));
        if (upper == 1) return 0;
        int threshold = 256 % upper, product;
        do
        {
            if (bytesRemaining == 0) { byteBuffer = Next32(); bytesRemaining = 4; }
            int value = (int)(byteBuffer & 255); byteBuffer >>= 8; bytesRemaining--;
            product = value * upper;
        } while ((product & 255) < threshold);
        return product >> 8;
    }
    internal double NextOpenDouble() => ((Next64() >> 12) + .5) / 4503599627370496.0;
    internal double NextNormal()
    {
        if (hasNormal) { hasNormal = false; return cachedNormal; }
        double radius = Math.Sqrt(-2 * Math.Log(NextOpenDouble())), angle = 2 * Math.PI * NextOpenDouble();
        cachedNormal = radius * Math.Sin(angle); hasNormal = true;
        return radius * Math.Cos(angle);
    }
}


