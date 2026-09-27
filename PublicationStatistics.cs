using System.Globalization;
using System.Numerics;
using System.Text.Json;
using Result = System.Collections.Generic.Dictionary<string, object?>;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

/// <summary>
/// Recomputes the consolidated manuscript's H4/H10 sign families and pressure/charged H6
/// inference from frozen predictions. This does not fit or evolve physical fields.
/// </summary>
public static partial class PublicationStatistics
{
    const int Draws = 99999;
    static readonly JsonSerializerOptions JsonOptions = new() { WriteIndented = true };
    sealed record Observation(string Galaxy, int Index, double Radius, double Observed, double Sigma, bool Outer);
    sealed record Prediction(Observation Observation, string Model, double Value, bool Failed);
    sealed record Sample(string[] Models, string[] Galaxies, int[] Counts, double[,] Loss, double[,] AbsZ,
        double[,] AbsKms, double[,] MeanZ, double[,] VarianceZ, bool[,] Failed, Prediction[] Predictions);

    public static Result Run(string dataRoot, string outputRoot)
    {
        string data = Directory.Exists(Path.Combine(dataRoot, "publication")) ? dataRoot : Path.GetDirectoryName(dataRoot)!;
        string source = Path.Combine(data, "publication", "statistics");
        string destination = Path.Combine(outputRoot, "publication_statistics");
        Directory.CreateDirectory(destination);
        var integrity = VerifyInputs(source);
        string original = Path.Combine(data, "numerics", "pointwise_comparison_2026_09_25", "data", "pointwise_all.csv");
        Require(Data.FileSha256(original) == "2bf5d6323cf92b232dfed60d04d2e7b4ed95f8eba0005f4865a4f883bd396052", "Original observation reference changed.");
        var observations = ReadObservations(original);
        using var rng = JsonDocument.Parse(File.ReadAllText(Path.Combine(source, "numpy_rng_reference.json")));
        var randomChecks = PublicationStatisticsRandom.SelfCheck(rng.RootElement);
        var policyChecks = SelfCheck(rng.RootElement);
        var envelope = RunEnvelope(source, destination, original, observations);
        var pressure = RunPhysical("pressure", source, destination, observations, rng.RootElement);
        var charged = RunPhysical("charged", source, destination, observations, rng.RootElement);
        var report = new Result
        {
            ["status"] = "pass", ["analysis"] = "consolidated_manuscript_2026_09_26",
            ["statistical_unit"] = "one galaxy; original 131 galaxies and 659 outer measurements",
            ["scope"] = "Retrospective inference conditional on frozen predictions. No fit, population formation, physical stability, nuisance propagation, survey covariance or model-history correction is implied.",
            ["failure_policy"] = "Finite prediction wins over failure in signs. Undefined full-sample continuous contrasts have no p value; only jointly feasible subset intervals are reported. Conservative p=1 placeholders are internal to the fixed six-slot Holm adjustment.",
            ["monte_carlo"] = "99999 whole-galaxy pairs draws and 99999 Rademacher wild draws; plus-one p values and explicit zero-exceedance floor flags",
            ["input_integrity"] = integrity, ["random_stream_validation"] = randomChecks, ["policy_validation"] = policyChecks,
            ["envelope"] = envelope, ["pressure"] = pressure, ["charged"] = charged
        };
        WriteJson(Path.Combine(destination, "validation.json"), report);
        File.WriteAllText(Path.Combine(destination, "README.md"),
            "# Consolidated publication statistics\n\nNative .NET reconstruction from frozen predictions; no fits or field evolutions are performed. " +
            "Each family directory contains recomputed tables and a comparison against the retained independent Python outputs.\n\n" +
            "Envelope H4/H10 includes the declared 1e-4 tie rule and the 1e-8 sensitivity. " +
            "The historical exact-zero sensitivity is recomputed from the preserved per-galaxy score serialization, after auditing those scores against prediction rows; exact-zero signs are sensitive to floating-point serialization.\n\n" +
            "Pressure and charged H6 tables retain all six declared slots, including the quartic/zero-contact controls. " +
            "The failed profiled-baryon fit UGC01281 remains a failure in signs; its continuous comparison is explicitly conditional on 130 feasible galaxies and has no primary p value. " +
            "Intervals are marginal bootstrap-t intervals, not simultaneous intervals or inversions of Holm-adjusted wild tests.\n\n" +
            "Probabilities and unsigned Gaussian equivalents are nominal retrospective diagnostics; neither measures detection of the cubic mechanism.\n");
        return report;
    }

    static Result VerifyInputs(string source)
    {
        using var manifest = JsonDocument.Parse(File.ReadAllText(Path.Combine(source, "input_manifest.json")));
        int count = 0;
        foreach (var item in manifest.RootElement.GetProperty("files").EnumerateArray())
        {
            string relative = item.GetProperty("path").GetString()!;
            string path = Path.GetFullPath(Path.Combine(source, relative));
            Require(path.StartsWith(Path.GetFullPath(source) + Path.DirectorySeparatorChar, StringComparison.OrdinalIgnoreCase), "Input manifest path escapes statistics data.");
            Require(File.Exists(path) && Data.FileSha256(path) == item.GetProperty("sha256").GetString(), $"Frozen publication input changed or missing: {relative}");
            count++;
        }
        return new() { ["status"] = "pass", ["verified_files"] = count };
    }

    static Dictionary<(string, int), Observation> ReadObservations(string path)
    {
        var result = new Dictionary<(string, int), Observation>();
        foreach (var row in Csv.Read(path))
        {
            var obs = ParseObservation(row);
            if (result.TryGetValue((obs.Galaxy, obs.Index), out var old)) Require(EqualObservation(old, obs), "Original models disagree on observations.");
            else result.Add((obs.Galaxy, obs.Index), obs);
        }
        Require(result.Count == 3034 && result.Values.Count(o => o.Outer) == 659 && result.Values.Select(o => o.Galaxy).Distinct().Count() == 131, "Incorrect original sample.");
        return result;
    }

    static Sample ReadSample(string path, string[] models, Dictionary<(string, int), Observation> reference)
    {
        var predictions = new List<Prediction>();
        var keys = new HashSet<(string, int, string)>();
        foreach (var row in Csv.Read(path))
        {
            string model = row["model"];
            if (!models.Contains(model)) continue;
            var obs = ParseObservation(row);
            Require(reference.TryGetValue((obs.Galaxy, obs.Index), out var expected) && EqualObservation(obs, expected), $"Changed observation metadata or split: {obs.Galaxy}/{obs.Index}.");
            Require(keys.Add((obs.Galaxy, obs.Index, model)), $"Duplicate prediction: {obs.Galaxy}/{obs.Index}/{model}.");
            double prediction = Number(row, "predicted", "predicted_kms");
            bool failed = row.ContainsKey("fit_failed") && Flag(row["fit_failed"]);
            predictions.Add(new(obs, model, prediction, failed));
        }
        var supplied = predictions.Select(p => (p.Observation.Galaxy, p.Observation.Index)).ToHashSet();
        Require(supplied.Count is 3034 or 659, "Input must supply all 3034 observations or all 659 outer observations.");
        var expectedKeys = reference.Where(x => supplied.Count == 3034 || x.Value.Outer).Select(x => x.Key).ToHashSet();
        Require(supplied.SetEquals(expectedKeys), "Incomplete original observation coverage.");
        foreach (var model in models)
            Require(predictions.Count(p => p.Model == model) == supplied.Count, $"Incomplete prediction coverage for {model}.");
        string[] galaxies = reference.Values.Select(p => p.Galaxy).Distinct().Order(StringComparer.Ordinal).ToArray();
        var groups = predictions.GroupBy(p => (p.Observation.Galaxy, p.Model)).ToDictionary(g => g.Key, g => g.OrderBy(p => p.Observation.Index).ToArray());
        int G = galaxies.Length, M = models.Length;
        int[] counts = galaxies.Select(g => reference.Values.Count(o => o.Galaxy == g && o.Outer)).ToArray();
        double[,] loss = new double[G, M], absZ = new double[G, M], absKms = new double[G, M], meanZ = new double[G, M], varianceZ = new double[G, M];
        bool[,] failedBlocks = new bool[G, M];
        for (int g = 0; g < G; g++)
        for (int m = 0; m < M; m++)
        {
            var all = groups[(galaxies[g], models[m])];
            var outer = all.Where(p => p.Observation.Outer).ToArray();
            Require(outer.Length == counts[g], "Outer sample mismatch.");
            bool failed = BlockFailed(all);
            failedBlocks[g, m] = failed;
            double[] z = outer.Select(p => (p.Value - p.Observation.Observed) / p.Observation.Sigma).ToArray();
            loss[g, m] = failed ? double.PositiveInfinity : Kahan(z.Select(Square)) / z.Length;
            absZ[g, m] = failed ? double.PositiveInfinity : Kahan(z.Select(Math.Abs)) / z.Length;
            absKms[g, m] = failed ? double.PositiveInfinity : Kahan(outer.Select(p => Math.Abs(p.Value - p.Observation.Observed))) / z.Length;
            meanZ[g, m] = Kahan(z) / z.Length;
            double mu = meanZ[g, m];
            varianceZ[g, m] = Kahan(z.Select(v => Square(v - mu))) / z.Length;
        }
        return new(models, galaxies, counts, loss, absZ, absKms, meanZ, varianceZ, failedBlocks, predictions.ToArray());
    }

    static Observation ParseObservation(Row row)
    {
        double index = Number(row, "index", "point_index");
        Require(double.IsFinite(index) && index >= 0 && index == Math.Floor(index), "Invalid observation index.");
        double radius = Number(row, "radius", "radius_kpc"), obs = Number(row, "observed", "observed_kms"), sigma = Number(row, "sigma", "sigma_kms");
        Require(double.IsFinite(radius) && radius > 0 && double.IsFinite(obs) && double.IsFinite(sigma) && sigma > 0, "Invalid observation coordinate.");
        return new(row["galaxy"], checked((int)index), radius, obs, sigma, Flag(row["is_outer"]));
    }

    static bool EqualObservation(Observation a, Observation b) => a.Galaxy == b.Galaxy && a.Index == b.Index && a.Outer == b.Outer
        && Near(a.Radius, b.Radius, 1e-12, 1e-12) && Near(a.Observed, b.Observed, 1e-12, 1e-12) && Near(a.Sigma, b.Sigma, 1e-12, 1e-12);
    static bool Flag(string value) => value.Trim().ToLowerInvariant() switch
    {
        "true" or "1" or "1.0" => true, "false" or "0" or "0.0" => false,
        _ => throw new InvalidDataException($"Invalid boolean {value}.")
    };
    static double Number(Row row, string key, string alias) => Csv.Number(row, row.ContainsKey(key) ? key : alias);
    static bool Near(double a, double b, double absolute = 2e-10, double relative = 2e-10) =>
        a == b || double.IsFinite(a) && double.IsFinite(b) && Math.Abs(a - b) <= absolute + relative * Math.Abs(b);
    static double Square(double value) => value * value;
    static double Kahan(IEnumerable<double> values)
    {
        double sum = 0, correction = 0;
        foreach (double value in values) { double adjusted = value - correction, next = sum + adjusted; correction = next - sum - adjusted; sum = next; }
        return sum;
    }
    static void Require(bool condition, string message) { if (!condition) throw new InvalidDataException(message); }
    static void WriteJson(string path, object value) => File.WriteAllText(path, JsonSerializer.Serialize(value, JsonOptions) + Environment.NewLine);

    static Result Sign(double[] left, double[] right, double tolerance)
    {
        int wins = 0, losses = 0, ties = 0, joint = 0;
        for (int i = 0; i < left.Length; i++)
        {
            bool af = double.IsFinite(left[i]), bf = double.IsFinite(right[i]);
            if (!af && !bf) joint++;
            else if (!af) losses++;
            else if (!bf) wins++;
            else if (Math.Abs(left[i] - right[i]) <= tolerance * (1 + Math.Max(left[i], right[i]))) ties++;
            else if (left[i] < right[i]) wins++;
            else losses++;
        }
        int n = wins + losses;
        BigInteger term = BigInteger.One, numerator = BigInteger.One;
        for (int k = 1; k <= Math.Min(wins, losses); k++) { term = term * (n - k + 1) / k; numerator += term; }
        double p = n == 0 ? 1 : Math.Min(1, (double)numerator / Math.Pow(2, n - 1));
        return new() { ["wins"] = wins, ["losses"] = losses, ["ties"] = ties, ["joint_failures"] = joint,
            ["informative_signs"] = n, ["p_sign_raw"] = p };
    }

    static Result RunEnvelope(string source, string destination, string original, Dictionary<(string, int), Observation> observations)
    {
        string[] models = ["extended", "adaptive_primary", "baryon_only", "scalar_core", "compact_plummer", "nfw", "profiled_baryons"];
        var sample = ReadSample(original, models, observations);
        var saved = Csv.Read(Path.Combine(source, "envelope", "saved_galaxy_scores.csv"))
            .Where(r => models.Contains(r["model"])).ToDictionary(r => (r["galaxy"], r["model"]), r => Csv.Number(r, "outer_loss"));
        double maxScoreError = 0;
        for (int g = 0; g < sample.Galaxies.Length; g++)
        for (int m = 0; m < models.Length - 1; m++)
        {
            double expected = saved[(sample.Galaxies[g], models[m])];
            Require(Near(sample.Loss[g, m], expected), "Saved envelope galaxy loss disagrees with raw predictions.");
            maxScoreError = Math.Max(maxScoreError, Math.Abs(sample.Loss[g, m] - expected));
        }
        (string Left, string Right, string Family)[] pairs =
        [
            ("extended","baryon_only","fixed_envelope"), ("extended","scalar_core","fixed_envelope"),
            ("extended","compact_plummer","fixed_envelope"), ("extended","nfw","fixed_envelope"),
            ("adaptive_primary","baryon_only","adaptive"), ("adaptive_primary","scalar_core","adaptive"),
            ("adaptive_primary","extended","adaptive"), ("adaptive_primary","nfw","adaptive"),
            ("extended","profiled_baryons","added_profiled_baryon_sensitivity"), ("adaptive_primary","profiled_baryons","added_profiled_baryon_sensitivity")
        ];
        var rows = new List<Result>();
        foreach (var (variant, tolerance) in new[] { ("primary_symmetric_1e-4", 1e-4), ("sensitivity_symmetric_1e-8", 1e-8), ("sensitivity_exact_zero", 0.0) })
        {
            var familyRows = new List<Result>();
            foreach (var pair in pairs)
            {
                int left = Array.IndexOf(models, pair.Left), right = Array.IndexOf(models, pair.Right);
                // Historical zero-tolerance signs depend on saved score roundoff. Keep that input
                // distinct from freshly recomputed squared losses, which are audited above.
                double[] a = sample.Galaxies.Select((g, i) => tolerance == 0 ? saved[(g, pair.Left)] : sample.Loss[i, left]).ToArray();
                double[] b = sample.Galaxies.Select((g, i) => tolerance == 0 && pair.Right != "profiled_baryons" ? saved[(g, pair.Right)] : sample.Loss[i, right]).ToArray();
                var sign = Sign(a, b, tolerance);
                var row = new Result { ["variant"] = variant, ["family_source"] = pair.Family, ["left"] = pair.Left, ["right"] = pair.Right,
                    ["wins"] = sign["wins"], ["losses"] = sign["losses"], ["ties"] = sign["ties"], ["joint_failures"] = sign["joint_failures"],
                    ["informative_n"] = sign["informative_signs"], ["raw_p"] = sign["p_sign_raw"], ["raw_Z"] = Statistics.NormalMagnitude((double)sign["p_sign_raw"]!),
                    ["previous_holm4_p"] = null, ["previous_holm4_Z"] = null,
                    ["comparator_infeasible_galaxies"] = b.Count(x => !double.IsFinite(x)), ["tie_relative_scale"] = tolerance,
                    ["holm10_p"] = null, ["holm10_Z"] = null,
                    ["inference_scope"] = "nominal_retrospective_family_sensitivity_not_history_adjusted" };
                familyRows.Add(row);
            }
            foreach (int start in new[] { 0, 4 })
            {
                double[] adjusted = Statistics.Holm(familyRows.Skip(start).Take(4).Select(r => (double)r["raw_p"]!).ToArray());
                for (int j = 0; j < 4; j++) { familyRows[start + j]["previous_holm4_p"] = adjusted[j]; familyRows[start + j]["previous_holm4_Z"] = Statistics.NormalMagnitude(adjusted[j]); }
            }
            double[] h10 = Statistics.Holm(familyRows.Select(r => (double)r["raw_p"]!).ToArray());
            for (int j = 0; j < 10; j++) { familyRows[j]["holm10_p"] = h10[j]; familyRows[j]["holm10_Z"] = Statistics.NormalMagnitude(h10[j]); }
            rows.AddRange(familyRows);
        }
        string output = Path.Combine(destination, "envelope");
        Directory.CreateDirectory(output);
        Csv.Write(Path.Combine(output, "sign_family_h10.csv"), rows);
        var originalRows = new List<Result>();
        foreach (var variant in rows.GroupBy(r => (string)r["variant"]!))
        {
            var originalFamily = variant.Where(r => (string)r["family_source"]! != "added_profiled_baryon_sensitivity").ToArray();
            double[] h8 = Statistics.Holm(originalFamily.Select(r => (double)r["raw_p"]!).ToArray());
            for (int j = 0; j < originalFamily.Length; j++)
            {
                var row = originalFamily[j];
                string left = (string)row["left"]!, right = (string)row["right"]!;
                // Mean differences remain fresh reconstructions; frozen scores are used only
                // by the explicitly labeled exact-zero sign diagnostic above.
                int lm = Array.IndexOf(models, left), rm = Array.IndexOf(models, right);
                double[] delta = Enumerable.Range(0, sample.Galaxies.Length).Select(g => sample.Loss[g, lm] - sample.Loss[g, rm]).ToArray();
                originalRows.Add(new() { ["family"] = row["family_source"], ["variant"] = row["variant"],
                    ["left"] = left, ["right"] = right, ["total_galaxies"] = sample.Galaxies.Length, ["outer_rows"] = sample.Counts.Sum(),
                    ["wins"] = row["wins"], ["losses"] = row["losses"], ["ties"] = row["ties"], ["informative_n"] = row["informative_n"],
                    ["win_fraction_all"] = (int)row["wins"]! / (double)sample.Galaxies.Length,
                    ["win_fraction_non_ties"] = (int)row["wins"]! / (double)(int)row["informative_n"]!,
                    ["raw_p"] = row["raw_p"], ["raw_Z"] = row["raw_Z"], ["tie_relative_scale"] = row["tie_relative_scale"],
                    ["mean_delta"] = delta.Average(), ["median_delta"] = Data.Median(delta),
                    ["holm4_p"] = row["previous_holm4_p"], ["holm4_Z"] = row["previous_holm4_Z"],
                    ["holm8_p_sensitivity"] = h8[j], ["holm8_Z_sensitivity"] = Statistics.NormalMagnitude(h8[j]) });
            }
        }
        Csv.Write(Path.Combine(output, "sign_family_h4.csv"), originalRows);
        Csv.Write(Path.Combine(output, "model_galaxy_squared_losses.csv"), Enumerable.Range(0, sample.Galaxies.Length).Select(g =>
        {
            var row = new Result { ["galaxy"] = sample.Galaxies[g] };
            for (int m = 0; m < models.Length; m++) row[models[m]] = double.IsFinite(sample.Loss[g, m]) ? sample.Loss[g, m] : "inf";
            return row;
        }));
        var audit = CompareReference(Path.Combine(source, "envelope", "reference", "sign_family_h10.csv"), rows,
            ["variant", "left", "right"], ["wins", "losses", "ties", "informative_n", "raw_p", "previous_holm4_p", "previous_holm4_Z", "holm10_p", "holm10_Z"]);
        var originalAudit = CompareReference(Path.Combine(source, "envelope", "reference", "sign_family_original.csv"), originalRows,
            ["variant", "left", "right"]);
        var report = new Result { ["status"] = "pass", ["comparisons"] = rows.Count, ["maximum_reconstructed_score_error"] = maxScoreError,
            ["reference_comparison"] = audit, ["original_h4_h8_reference_comparison"] = originalAudit,
            ["zero_tolerance_source"] = "Preserved galaxy scores, independently checked against raw prediction rows" };
        WriteJson(Path.Combine(output, "checks.json"), report);
        return report;
    }

    static Result RunPhysical(string family, string source, string destination, Dictionary<(string, int), Observation> observations, JsonElement rng)
    {
        string input = Path.Combine(source, family), output = Path.Combine(destination, family);
        Directory.CreateDirectory(output);
        using var document = JsonDocument.Parse(File.ReadAllText(Path.Combine(input, "PROTOCOL.json")));
        var protocol = document.RootElement;
        string primary = protocol.GetProperty("primary_model").GetString()!;
        string[] comparators = protocol.GetProperty("primary_comparators").EnumerateArray().Select(v => v.GetString()!).ToArray();
        Require(comparators.Length == 6 && protocol.GetProperty("bootstrap_draws").GetInt32() == Draws, "Changed publication protocol.");
        ulong pairsSeed = protocol.GetProperty("pairs_seed").GetUInt64(), wildSeed = protocol.GetProperty("wild_seed").GetUInt64(), subsetSeed = protocol.GetProperty("finite_subset_seed").GetUInt64();
        string[] models = [primary, .. comparators];
        var s = ReadSample(Path.Combine(input, "predicted_rows.csv"), models, observations);
        int G = s.Galaxies.Length;
        var modelRows = new List<Result>();
        for (int m = 0; m < models.Length; m++)
        {
            int[] finite = Enumerable.Range(0, G).Where(g => !s.Failed[g, m]).ToArray();
            double[] losses = finite.Select(g => s.Loss[g, m]).ToArray();
            modelRows.Add(new() { ["model"] = models[m], ["total_galaxies"] = G, ["total_outer_points"] = s.Counts.Sum(),
                ["failed_galaxies"] = G - finite.Length, ["failed_outer_points"] = Enumerable.Range(0, G).Where(g => s.Failed[g, m]).Sum(g => s.Counts[g]),
                ["failed_galaxy_names"] = string.Join(';', Enumerable.Range(0, G).Where(g => s.Failed[g, m]).Select(g => s.Galaxies[g])),
                ["finite_galaxies"] = finite.Length, ["finite_outer_points"] = finite.Sum(g => s.Counts[g]),
                ["full_sample_mean_squared_loss"] = finite.Length == G ? losses.Average() : null,
                ["finite_subset_mean_squared_loss"] = losses.Length > 0 ? losses.Average() : null,
                ["finite_subset_median_squared_loss"] = losses.Length > 0 ? Data.Median(losses) : null,
                ["finite_subset_mean_absolute_standardized"] = finite.Length > 0 ? finite.Average(g => s.AbsZ[g, m]) : null,
                ["finite_subset_mean_absolute_kms"] = finite.Length > 0 ? finite.Average(g => s.AbsKms[g, m]) : null,
                ["finite_subset_galaxies_rms_z_above3"] = losses.Count(x => x > 9),
                ["conditioning"] = finite.Length == G ? "all131" : "explicitly_conditional_on_feasible_galaxies" });
        }
        Csv.Write(Path.Combine(output, "model_summary.csv"), modelRows);
        Csv.Write(Path.Combine(output, "model_galaxy_squared_losses.csv"), Enumerable.Range(0, G).Select(g =>
        {
            var row = new Result { ["galaxy"] = s.Galaxies[g] };
            for (int m = 0; m < models.Length; m++) row[models[m]] = double.IsFinite(s.Loss[g, m]) ? s.Loss[g, m] : "inf";
            return row;
        }));
        Csv.Write(Path.Combine(output, "model_galaxy_failures.csv"), Enumerable.Range(0, G).Select(g =>
        {
            var row = new Result { ["galaxy"] = s.Galaxies[g] };
            for (int m = 0; m < models.Length; m++) row[models[m]] = s.Failed[g, m];
            return row;
        }));
        var signs = new List<Result>();
        for (int j = 0; j < comparators.Length; j++)
        {
            var row = new Result { ["primary"] = primary, ["comparator"] = comparators[j] };
            foreach (var item in Sign(Enumerable.Range(0, G).Select(g => s.Loss[g, 0]).ToArray(), Enumerable.Range(0, G).Select(g => s.Loss[g, j + 1]).ToArray(), 1e-4)) row[item.Key] = item.Value;
            signs.Add(row);
        }
        double[] h6 = Statistics.Holm(signs.Select(r => (double)r["p_sign_raw"]!).ToArray());
        for (int j = 0; j < comparators.Length; j++)
        {
            signs[j]["p_sign_holm6"] = h6[j]; signs[j]["Z_sign_holm6"] = Statistics.NormalMagnitude(h6[j]);
            signs[j]["inference"] = "nominal_retrospective_independent_fair_sign_diagnostic";
        }
        Csv.Write(Path.Combine(output, "sign_comparisons.csv"), signs);
        int[] full = Enumerable.Range(0, comparators.Length).Where(j => Enumerable.Range(0, G).All(g => !s.Failed[g, 0] && !s.Failed[g, j + 1])).ToArray();
        Console.WriteLine($"Publication {family}: {G} galaxies, {s.Counts.Sum()} outer rows; {Draws:N0} pairs and wild draws.");
        double[][] fullD = full.Select(j => Enumerable.Range(0, G).Select(g => s.Loss[g, 0] - s.Loss[g, j + 1]).ToArray()).ToArray();
        var fullResults = full.Length > 0 ? Bootstrap(fullD, s.Counts, pairsSeed, wildSeed, rng) : [];
        var continuous = new List<Result>(); var influences = new List<Result>(); var contrasts = new List<Result>(); var covariance = new List<Result>();
        for (int j = 0; j < comparators.Length; j++)
        {
            int m = j + 1;
            int[] feasible = Enumerable.Range(0, G).Where(g => !s.Failed[g, 0] && !s.Failed[g, m]).ToArray();
            double[] d = feasible.Select(g => s.Loss[g, 0] - s.Loss[g, m]).ToArray();
            bool all = feasible.Length == G;
            Result row;
            if (all) row = fullResults[Array.IndexOf(full, j)];
            else
            {
                Require(feasible.Length >= 2, "Too few jointly feasible galaxies for saved publication comparison.");
                row = Bootstrap([d], feasible.Select(g => s.Counts[g]).ToArray(), subsetSeed + (ulong)j, null, rng)[0];
            }
            row["full_sample_continuous_defined"] = all; row["primary_p_assigned"] = all;
            row["primary"] = primary; row["comparator"] = comparators[j]; row["n_galaxies"] = feasible.Length;
            row["n_points"] = feasible.Sum(g => s.Counts[g]);
            row["excluded_galaxies"] = string.Join(';', Enumerable.Range(0, G).Except(feasible).Select(g => s.Galaxies[g]));
            row["conditioning"] = all ? "all131" : "jointly_feasible_subset_no_primary_p";
            row["mean_absolute_kms_difference"] = feasible.Average(g => s.AbsKms[g, 0] - s.AbsKms[g, m]);
            row["mean_absolute_standardized_difference"] = feasible.Average(g => s.AbsZ[g, 0] - s.AbsZ[g, m]);
            row["pairs_seed_used"] = all ? pairsSeed : subsetSeed + (ulong)j; row["wild_seed_used"] = all ? wildSeed : null;
            row["p_wild_holm6"] = null; row["Z_wild_holm6"] = null;
            row["inference"] = "nominal_retrospective_fixed_prediction_cluster_mean_loss";
            row["reason"] = all ? null : "Failed fit makes full-sample continuous mean inference undefined; finite-pair sensitivity only";
            continuous.Add(row);
            double mu = d.Average(), sum = d.Sum(), denom = d.Sum(x => Square(x - mu));
            for (int i = 0; i < feasible.Length; i++)
                influences.Add(new() { ["comparator"] = comparators[j], ["galaxy"] = s.Galaxies[feasible[i]], ["delta_galaxy_loss"] = d[i],
                    ["variance_share"] = denom == 0 ? 0 : Square(d[i] - mu) / denom,
                    ["mean_without"] = (sum - d[i]) / (d.Length - 1), ["conditioning"] = row["conditioning"] });
            for (int g = 0; g < G; g++)
                contrasts.Add(new() { ["galaxy"] = s.Galaxies[g], ["comparator"] = comparators[j],
                    ["left_loss"] = double.IsFinite(s.Loss[g, 0]) ? s.Loss[g, 0] : "inf", ["right_loss"] = double.IsFinite(s.Loss[g, m]) ? s.Loss[g, m] : "inf",
                    ["delta"] = s.Failed[g, 0] || s.Failed[g, m] ? null : s.Loss[g, 0] - s.Loss[g, m],
                    ["left_failed"] = s.Failed[g, 0], ["right_failed"] = s.Failed[g, m] });
            foreach (double rho in new[] { 0.0, .25, .5, .75 })
            {
                double Score(int g, int k) => s.VarianceZ[g, k] / (1 - rho) + Square(s.MeanZ[g, k]) / (1 + (s.Counts[g] - 1) * rho);
                double[] delta = feasible.Select(g => Score(g, 0) - Score(g, m)).ToArray();
                covariance.Add(new() { ["comparator"] = comparators[j], ["assumed_rho"] = rho, ["n_galaxies"] = feasible.Length,
                    ["mean_difference"] = delta.Average(), ["median_difference"] = Data.Median(delta),
                    ["conditioning"] = row["conditioning"], ["new_p_assigned"] = false });
            }
        }
        double[] adjusted = Statistics.Holm(continuous.Select(r => (bool)r["primary_p_assigned"]! ? (double)r["p_wild_raw"]! : 1).ToArray());
        for (int j = 0; j < comparators.Length; j++)
            if ((bool)continuous[j]["primary_p_assigned"]!) { continuous[j]["p_wild_holm6"] = adjusted[j]; continuous[j]["Z_wild_holm6"] = Statistics.NormalMagnitude(adjusted[j]); }
        Csv.Write(Path.Combine(output, "continuous_comparisons.csv"), continuous);
        Csv.Write(Path.Combine(output, "galaxy_contrasts.csv"), contrasts);
        Csv.Write(Path.Combine(output, "galaxy_influence.csv"), influences);
        Csv.Write(Path.Combine(output, "covariance_sensitivity.csv"), covariance);
        WritePointDifferences(output, s);
        var audits = new Result
        {
            ["sign"] = CompareReference(Path.Combine(input, "reference", "sign_comparisons.csv"), signs, ["comparator"], ["wins","losses","ties","joint_failures","informative_signs","p_sign_raw","p_sign_holm6","Z_sign_holm6"]),
            ["continuous"] = CompareReference(Path.Combine(input, "reference", "continuous_comparisons.csv"), continuous, ["comparator"]),
            ["models"] = CompareReference(Path.Combine(input, "reference", "model_summary.csv"), modelRows, ["model"]),
            ["covariance"] = CompareReference(Path.Combine(input, "reference", "covariance_sensitivity.csv"), covariance, ["comparator","assumed_rho"])
        };
        var result = new Result { ["status"] = "pass", ["primary"] = primary, ["comparators"] = comparators, ["galaxies"] = G,
            ["outer_rows"] = s.Counts.Sum(), ["supplied_observation_rows"] = s.Predictions.Length / models.Length,
            ["metadata_unchanged"] = true, ["split_unchanged"] = true, ["bootstrap_draws"] = Draws,
            ["pairs_seed"] = pairsSeed, ["wild_seed"] = wildSeed, ["finite_subset_seed_base"] = subsetSeed,
            ["reference_comparisons"] = audits };
        WriteJson(Path.Combine(output, "checks.json"), result);
        File.Copy(Path.Combine(input, "PROTOCOL.json"), Path.Combine(output, "PROTOCOL_USED.json"), true);
        return result;
    }

    static void WritePointDifferences(string output, Sample s)
    {
        var predictions = s.Predictions.Where(p => p.Observation.Outer).ToDictionary(p => (p.Observation.Galaxy, p.Observation.Index, p.Model));
        var rows = new List<Result>();
        for (int m = 1; m < s.Models.Length; m++)
        foreach (var primary in s.Predictions.Where(p => p.Observation.Outer && p.Model == s.Models[0]))
        {
            var obs = primary.Observation;
            var comp = predictions[(obs.Galaxy, obs.Index, s.Models[m])];
            int g = Array.IndexOf(s.Galaxies, obs.Galaxy);
            bool valid = !s.Failed[g, 0] && !s.Failed[g, m];
            double a = (primary.Value - obs.Observed) / obs.Sigma, b = (comp.Value - obs.Observed) / obs.Sigma;
            rows.Add(new() { ["galaxy"] = obs.Galaxy, ["index"] = obs.Index, ["radius"] = obs.Radius, ["comparator"] = s.Models[m],
                ["z_primary"] = double.IsFinite(a) ? a : null, ["z_comparator"] = double.IsFinite(b) ? b : null,
                ["delta_squared"] = valid ? a * a - b * b : null, ["delta_absolute_kms"] = valid ? (Math.Abs(a) - Math.Abs(b)) * obs.Sigma : null, ["comparison_valid"] = valid });
        }
        Csv.Write(Path.Combine(output, "pointwise_paired_differences.csv"), rows);
    }
}



