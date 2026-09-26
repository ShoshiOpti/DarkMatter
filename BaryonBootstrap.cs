using System.Globalization;
using System.Numerics;
using Result = System.Collections.Generic.Dictionary<string, object?>;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

public static class BaryonBootstrap
{
    public const int Seed = 20260924, Resamples = 20000;
    const string InitialState = "87917040079043241584574165612868195315", Increment = "77762782066293002654591666826778350581";
    static readonly string[] Effects = ["mean_loss_reduction", "median_loss_reduction", "mean_rms_reduction", "median_rms_reduction",
        "ratio_of_medians_reduction", "mean_fractional_loss_reduction", "median_fractional_loss_reduction"];
    static readonly (string Augmented, string Baseline)[] Pairs = [("A_core", "bar_fixed"), ("A_nfw", "bar_fixed"),
        ("B_core", "bar_profiled"), ("C_core", "bar_profiled"), ("B_nfw", "bar_profiled")];

    public static Result Run(string dataRoot, string outputRoot)
    {
        SelfCheck();
        string numerics = Directory.Exists(Path.Combine(dataRoot, "numerics")) ? Path.Combine(dataRoot, "numerics") : dataRoot;
        string folder = Path.Combine(outputRoot, "baryons_only");
        var data = Csv.Read(Path.Combine(folder, "paired_galaxy_scores.csv"));
        var results = new List<Result>();
        foreach (var (mode, region) in new[] { ("full", "all"), ("inner", "outer") })
        {
            foreach (var (augmented, baseline) in Pairs)
            {
                var all = data.Where(r => r["augmented"] == augmented && r["baseline"] == baseline && r["fit_mode"] == mode && r["region"] == region)
                    .OrderBy(r => r["galaxy"], StringComparer.Ordinal).ToArray();
                if (all.Length != (mode == "full" ? 153 : 131) || all.Select(r => r["galaxy"]).Distinct().Count() != all.Length)
                    throw new InvalidDataException("Missing or repeated paired bootstrap galaxies.");
                int baselineFailures = all.Count(r => !double.IsFinite(Csv.Number(r, "baseline_loss"))), augmentedFailures = all.Count(r => !double.IsFinite(Csv.Number(r, "augmented_loss")));
                foreach (string sample in new[] { "all_selected", "common_finite" })
                {
                    var selected = sample == "all_selected" ? all : all.Where(r => double.IsFinite(Csv.Number(r, "augmented_loss")) && double.IsFinite(Csv.Number(r, "baseline_loss"))).ToArray();
                    double[] augmentedLoss = selected.Column("augmented_loss"), baselineLoss = selected.Column("baseline_loss"), augmentedRms = selected.Column("augmented_rms_kms"), baselineRms = selected.Column("baseline_rms_kms");
                    var row = new Result
                    {
                        ["augmented"] = augmented,
                        ["baseline"] = baseline,
                        ["fit_mode"] = mode,
                        ["region"] = region,
                        ["sample"] = sample,
                        ["n_total"] = all.Length,
                        ["n_galaxies"] = selected.Length,
                        ["original_baseline_failures"] = baselineFailures,
                        ["original_augmented_failures"] = augmentedFailures,
                        ["augmented_median_loss"] = Median((double[])augmentedLoss.Clone()),
                        ["baseline_median_loss"] = Median((double[])baselineLoss.Clone()),
                        ["augmented_mean_loss"] = augmentedLoss.Average(),
                        ["baseline_mean_loss"] = baselineLoss.Average(),
                        ["augmented_median_rms_kms"] = Median((double[])augmentedRms.Clone()),
                        ["baseline_median_rms_kms"] = Median((double[])baselineRms.Clone())
                    };
                    foreach (var item in Sign(augmentedLoss, baselineLoss))
                        row[item.Key] = item.Value;
                    foreach (var item in Bootstrap(augmentedLoss, baselineLoss, augmentedRms, baselineRms))
                        row[item.Key] = item.Value;
                    results.Add(row);
                }
            }
        }

        string output = Path.Combine(folder, "paired_summary.csv");
        Csv.Write(output, results);
        var verification = Statistics.ValidateCsv(Path.Combine(numerics, "baryons_only", "paired_summary.csv"), output);
        verification["seed"] = Seed;
        verification["resamples_per_pair"] = Resamples;
        verification["bootstrap_runs"] = results.Count;
        verification["pcg64_initial_state"] = InitialState;
        verification["pcg64_increment"] = Increment;
        verification["integer_generator"] = "NumPy PCG64 XSL-RR 128/64 with cached uint32 and Lemire bounded integers";
        verification["interval_method"] = "Paired galaxy percentile intervals, inverse empirical CDF";
        Data.SaveJson(Path.Combine(folder, "bootstrap_verification.json"), verification);
        return verification;
    }

    static Result Sign(double[] augmentedLoss, double[] baselineLoss)
    {
        var summary = Statistics.BaryonSign(augmentedLoss, baselineLoss);
        return new()
        {
            ["n_informative"] = summary["n_informative"],
            ["augmented_wins"] = summary["augmented_wins"],
            ["baseline_wins"] = summary["baseline_wins"],
            ["finite_ties"] = summary["finite_ties"],
            ["augmented_failures"] = summary["augmented_failures"],
            ["baseline_failures"] = summary["baseline_failures"],
            ["joint_failures"] = summary["joint_failures"],
            ["direction"] = summary["direction"],
            ["win_fraction"] = summary["win_fraction"] ?? double.NaN,
            ["win_ci_low"] = summary["win_ci_low"] ?? double.NaN,
            ["win_ci_high"] = summary["win_ci_high"] ?? double.NaN
        };
    }

    static Result Bootstrap(double[] augmentedLoss, double[] baselineLoss, double[] augmentedRms, double[] baselineRms)
    {
        int galaxyCount = augmentedLoss.Length;
        var indices = Enumerable.Range(0, galaxyCount).ToArray();
        var work = Enumerable.Range(0, 5).Select(_ => new double[galaxyCount]).ToArray();
        var estimate = Evaluate(augmentedLoss, baselineLoss, augmentedRms, baselineRms, indices, work);
        var samples = Effects.Select(_ => new double[Resamples]).ToArray();
        var random = new Pcg64();
        for (int draw = 0; draw < Resamples; draw++)
        {
            // Resample paired galaxies, keeping both models' losses and RMS values together.
            for (int i = 0; i < galaxyCount; i++)
                indices[i] = random.Next(galaxyCount);
            var values = Evaluate(augmentedLoss, baselineLoss, augmentedRms, baselineRms, indices, work);
            for (int effect = 0; effect < Effects.Length; effect++)
                samples[effect][draw] = values[effect];
        }
        var result = new Result();
        for (int effect = 0; effect < Effects.Length; effect++)
            result[Effects[effect]] = estimate[effect];
        for (int effect = 0; effect < Effects.Length; effect++)
        {
            bool invalid = samples[effect].Any(double.IsNaN);
            Array.Sort(samples[effect]);
            // This historical analysis uses inverse empirical-CDF ranks, not interpolated percentiles.
            result[Effects[effect] + "_ci_low"] = invalid ? double.NaN : samples[effect][(int)Math.Ceiling(.025 * Resamples) - 1];
            result[Effects[effect] + "_ci_high"] = invalid ? double.NaN : samples[effect][(int)Math.Ceiling(.975 * Resamples) - 1];
        }
        return result;
    }

    static double[] Evaluate(double[] augmentedLoss, double[] baselineLoss, double[] augmentedRms, double[] baselineRms, int[] indices, double[][] work)
    {
        var sampledAugmented = work[0];
        var sampledBaseline = work[1];
        var lossReductions = work[2];
        var rmsReductions = work[3];
        var fractionalReductions = work[4];
        double lossSum = 0, rmsSum = 0, fractionSum = 0;
        for (int i = 0; i < indices.Length; i++)
        {
            int galaxy = indices[i];
            double loss = baselineLoss[galaxy] - augmentedLoss[galaxy];
            double rms = baselineRms[galaxy] - augmentedRms[galaxy];
            double fraction = double.IsFinite(augmentedLoss[galaxy]) && double.IsFinite(baselineLoss[galaxy]) && baselineLoss[galaxy] > 0
                ? 1 - augmentedLoss[galaxy] / baselineLoss[galaxy] : double.NaN;
            sampledAugmented[i] = augmentedLoss[galaxy];
            sampledBaseline[i] = baselineLoss[galaxy];
            lossReductions[i] = loss;
            rmsReductions[i] = rms;
            fractionalReductions[i] = fraction;
            lossSum += loss;
            rmsSum += rms;
            fractionSum += fraction;
        }
        return [lossSum / indices.Length, Median(lossReductions), rmsSum / indices.Length, Median(rmsReductions),
            1 - Median(sampledAugmented) / Median(sampledBaseline), fractionSum / indices.Length, Median(fractionalReductions)];
    }

    static double Median(double[] values)
    {
        if (values.Length == 0 || values.Any(double.IsNaN))
            return double.NaN;
        Array.Sort(values);
        int middle = values.Length / 2;
        return values.Length % 2 == 1 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
    }

    public static void SelfCheck()
    {
        ulong[] reference = [1216218152932875416,2883288408839688800,11960178178589183471,9591654829523265231,770372349272430564,11566435092094653015,
            9296526231677078352,13211656360047194462,4258977636004141579,6273694152402553786,8013760359359257596,15422889583280318740];
        var raw = new Pcg64();
        if (!reference.All(value => raw.Next64() == value))
            throw new InvalidDataException("PCG64 raw sequence differs from NumPy.");
        int[] indices = [56, 10, 10, 23, 81, 99, 101, 79, 142, 6, 41, 95];
        var bounded = new Pcg64();
        if (!indices.All(value => bounded.Next(153) == value))
            throw new InvalidDataException("PCG64 bounded sequence differs from NumPy.");
    }

    sealed class Pcg64
    {
        UInt128 state = UInt128.Parse(InitialState, CultureInfo.InvariantCulture);
        readonly UInt128 increment = UInt128.Parse(Increment, CultureInfo.InvariantCulture);
        static readonly UInt128 Multiplier = ((UInt128)2549297995355413924UL << 64) | 4865540595714422341UL;
        uint cached;
        bool hasCached;
        public ulong Next64()
        {
            state = unchecked(state * Multiplier + increment);
            return BitOperations.RotateRight((ulong)(state >> 64) ^ (ulong)state, (int)(state >> 122));
        }
        uint Next32()
        {
            if (hasCached)
            {
                hasCached = false;
                return cached;
            }
            ulong value = Next64();
            cached = (uint)(value >> 32);
            hasCached = true;
            return (uint)value;
        }
        public int Next(int upper)
        {
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(upper);
            if (upper == 1)
                return 0;
            // Lemire bounded integers; uint32 caching matches NumPy PCG64.
            uint bound = (uint)upper, threshold = unchecked(0u - bound) % bound;
            ulong product;
            do
            {
                product = (ulong)Next32() * bound;
            } while ((uint)product < threshold);
            return (int)(product >> 32);
        }
    }
}
