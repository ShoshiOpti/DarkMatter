using System.Globalization;
using System.Text.Json;
using Row = System.Collections.Generic.Dictionary<string, string>;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class Statistics
{
    static Result Validate(string numerics, string destination, Result results)
    {
        string referencePath = Path.Combine(numerics, "significance", "results.json");
        if (!File.Exists(referencePath))
            throw new FileNotFoundException("Saved significance reference is required for validation.", referencePath);
        using var expected = JsonDocument.Parse(File.ReadAllText(referencePath));
        using var actual = JsonDocument.Parse(JsonSerializer.Serialize(JsonSafe(results)));
        int fields = 0;
        foreach (var property in expected.RootElement.EnumerateObject())
        {
            if (property.Name is "software_versions" or "implementation_sha256")
                continue;
            CompareJson(property.Value, actual.RootElement.GetProperty(property.Name), property.Name, ref fields);
        }
        var tables = new Result();
        foreach (string filename in new[] { "primary_summary.csv", "paired_galaxy_losses.csv", "prediction_reconstruction.csv" })
            tables[filename] = ValidateCsv(Path.Combine(numerics, "significance", filename), Path.Combine(destination, filename));
        return new Result { ["passed"] = true, ["reference"] = "numerics/significance/results.json", ["json_fields_checked"] = fields, ["csv_tables"] = tables, ["numeric_tolerance"] = "1e-10 absolute plus 1e-10 relative; probabilities use 1e-9 relative with 1e-300 absolute", ["excluded_metadata"] = new[] { "software_versions", "implementation_sha256" } };
    }

    static void CompareJson(JsonElement expected, JsonElement actual, string path, ref int fields)
    {
        if (expected.ValueKind == JsonValueKind.Number && actual.ValueKind == JsonValueKind.Number)
        {
            CheckNumber(expected.GetDouble(), actual.GetDouble(), path);
            fields++;
            return;
        }
        if (expected.ValueKind != actual.ValueKind)
            throw new InvalidDataException("Reference type mismatch: " + path);
        switch (expected.ValueKind)
        {
            case JsonValueKind.Object:
                if (expected.EnumerateObject().Count() != actual.EnumerateObject().Count())
                    throw new InvalidDataException("Reference property count mismatch: " + path);
                foreach (var item in expected.EnumerateObject())
                    CompareJson(item.Value, actual.GetProperty(item.Name), path + "." + item.Name, ref fields);
                break;
            case JsonValueKind.Array:
                if (expected.GetArrayLength() != actual.GetArrayLength())
                    throw new InvalidDataException("Reference array length mismatch: " + path);
                for (int i = 0; i < expected.GetArrayLength(); i++)
                    CompareJson(expected[i], actual[i], path + "[" + i + "]", ref fields);
                break;
            default:
                if (expected.GetRawText() != actual.GetRawText() && expected.ToString() != actual.ToString())
                    throw new InvalidDataException("Reference value mismatch: " + path);
                fields++;
                break;
        }
    }

    static void CheckNumber(double expected, double actual, string path)
    {
        if (double.IsNaN(expected) && double.IsNaN(actual) || expected == actual)
            return;
        bool probability = path.Split('.').Last().StartsWith("p_", StringComparison.Ordinal);
        double tolerance = probability ? 1e-300 + 1e-9 * Math.Abs(expected) : 1e-10 + 1e-10 * Math.Abs(expected);
        if (!double.IsFinite(actual) || !double.IsFinite(expected) || Math.Abs(expected - actual) > tolerance)
            throw new InvalidDataException($"Reference mismatch at {path}: expected {expected:R}, actual {actual:R}.");
    }

    internal static Result ValidateCsv(string reference, string produced)
    {
        var expected = Csv.Read(reference);
        var actual = Csv.Read(produced);
        if (expected.Count != actual.Count)
            throw new InvalidDataException("Reference CSV row count mismatch: " + reference);
        int fields = 0;
        for (int i = 0; i < expected.Count; i++)
        {
            if (!expected[i].Keys.ToHashSet().SetEquals(actual[i].Keys))
                throw new InvalidDataException("Reference CSV schema mismatch: " + reference);
            foreach (var item in expected[i])
            {
                string expectedText = item.Value;
                string actualText = actual[i][item.Key];
                string path = Path.GetFileName(reference) + "[" + i + "]." + item.Key;
                if (expectedText == actualText || string.Equals(expectedText, actualText, StringComparison.OrdinalIgnoreCase)
                    || (string.IsNullOrEmpty(expectedText) && actualText == "undefined"))
                {
                    fields++;
                    continue;
                }
                if (TryNumber(expectedText, out double expectedNumber) && TryNumber(actualText, out double actualNumber))
                    CheckNumber(expectedNumber, actualNumber, path);
                else
                    throw new InvalidDataException($"Reference CSV mismatch: {path}, '{expectedText}' != '{actualText}'.");
                fields++;
            }
        }
        return new Result { ["passed"] = true, ["rows"] = expected.Count, ["fields"] = fields };
    }

    static bool TryNumber(string text, out double value)
    {
        if (text is "inf" or "Infinity" or "positive_infinity")
        {
            value = double.PositiveInfinity;
            return true;
        }
        if (text is "-inf" or "-Infinity" or "negative_infinity")
        {
            value = double.NegativeInfinity;
            return true;
        }
        if (text is "nan" or "NaN" or "undefined" or "")
        {
            value = double.NaN;
            return true;
        }
        return double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
    }

    public static Result SelfChecks()
    {
        int checks = 0;
        void Require(bool condition, string name)
        {
            if (!condition)
                throw new InvalidDataException("Statistics self-check failed: " + name);
            checks++;
        }
        void Reject(Action action, string name)
        {
            bool rejected = false;
            try
            {
                action();
            }
            catch (Exception e) when (e is InvalidDataException or ArgumentException) { rejected = true; }
            Require(rejected, name);
        }
        static Row Loss(string galaxy, double loss, bool valid = true, int n = 2) => new()
        {
            ["galaxy"] = galaxy,
            ["n_test"] = n.ToString(CultureInfo.InvariantCulture),
            ["holdout_chi2_per_row"] = loss.ToString("R", CultureInfo.InvariantCulture),
            ["holdout_admissible"] = valid.ToString()
        };
        var a = ExactSign(0, 10);
        var b = ExactSign(10, 0);
        Require((double)a["p_two_sided"]! == 2.0 / 1024 && Equals(a["p_two_sided"], b["p_two_sided"]), "exact sign symmetry and boundary");
        Require(Math.Abs(((double[])a["left_win_fraction_ci95"]!)[1] - (1 - Math.Pow(.025, .1))) < 1e-13, "Clopper-Pearson boundary");
        Require(ExactSign(0, 0)["p_two_sided"] is null && ExactSign(0, 0)["left_win_fraction_ci95"] is null, "no information remains null");
        Require((double)ExactSign(50, 50)["p_two_sided"]! == 1 && (double)ExactSign(50, 50)["z_two_sided"]! == 0, "balanced signs");
        Reject(() => ExactSign(-1, 2), "negative counts rejected");
        Require(Holm([.03, .001, .02]).Zip(new[] { .04, .003, .04 }, (x, y) => Math.Abs(x - y) < 1e-15).All(x => x), "Holm multiplicity and order");
        Require(Holm([.8, .8, 1]).SequenceEqual(new[] { 1.0, 1.0, 1.0 }), "Holm caps at one");
        Reject(() => Holm([double.NaN, .1]), "invalid probability rejected");
        var interval = MedianInterval(Enumerable.Range(0, 131).Select(i => (double)i));
        Require((int)interval["k_1based"]! == 54 && ((double[])interval["ci95"]!).SequenceEqual(new double[] { 53, 77 }), "exact median interval ranks");
        Require(((double[])MedianInterval([1, 2])["ci95"]!)[0] == double.NegativeInfinity && MedianInterval([])["ci95"] is null, "small sample intervals");
        Require(double.IsNaN(Median([double.NegativeInfinity, double.PositiveInfinity])), "opposite infinite median is undefined");
        var left = new[] { Loss("a", 1), Loss("b", 0, false), Loss("c", 1), Loss("d", 0, false), Loss("e", 2) };
        var right = new[] { Loss("a", 2), Loss("b", 1), Loss("c", 1), Loss("d", 0, false), Loss("e", 0, false) };
        var (edge, rows) = Compare(left, right, "a", "b", "edge");
        Require((int)edge["left_wins"]! == 2 && (int)edge["right_wins"]! == 1 && (int)edge["finite_ties"]! == 1 && (int)edge["joint_failures"]! == 1 && (int)edge["n_informative"]! == 3, "finite ties and one-sided and joint failures");
        Require((string)rows[3]["outcome"]! == "joint_failure" && double.IsNaN((double)rows[3]["delta_left_minus_right"]!), "undefined joint difference");
        Reject(() => Compare([Loss("a", 1)], [Loss("b", 2)], "a", "b", "bad"), "unmatched galaxy rejected");
        Reject(() => Compare([Loss("a", 1), Loss("a", 2)], [Loss("a", 2)], "a", "b", "bad"), "duplicate galaxy rejected");
        Reject(() => Compare([Loss("a", 1)], [Loss("a", 2, n: 3)], "a", "b", "bad"), "unequal counts rejected");
        var fits = new[] { new Row { ["galaxy"] = "g", ["holdout"] = "True", ["n"] = "2", ["n_train"] = "1", ["n_test"] = "1", ["chi2_data_per_row"] = "1", ["holdout_chi2_per_row"] = "4", ["holdout_admissible"] = "True" } };
        var points = new[] { true, false }.Select((train, i) => new Row { ["galaxy"] = "g", ["holdout"] = "True", ["training_row"] = train.ToString(), ["r_catalogue_kpc"] = (i + 1).ToString(), ["Vobs_kms"] = "10", ["e_Vobs_kms"] = "2", ["Vpred_catalogue_kms"] = (12 + 2 * i).ToString() }).ToList();
        Require((int)AuditPredictions(fits, points).Report["partition_checks"]! == 2, "prediction loss reconstruction");
        var changed = points.Select(r => new Row(r)).ToList();
        changed[1]["Vpred_catalogue_kms"] = "15";
        Reject(() => AuditPredictions(fits, changed), "corrupted predictions rejected");
        changed = points.Select(r => new Row(r)).ToList();
        changed[1]["training_row"] = "True";
        Reject(() => AuditPredictions(fits, changed), "corrupted partition rejected");
        changed = points.Select(r => new Row(r)).ToList();
        changed[1]["e_Vobs_kms"] = "3";
        Reject(() => CheckPairedObservations(points, changed), "changed observations rejected");
        var full = new Row { ["galaxy"] = "g", ["n_train"] = "2", ["n_test"] = "0", ["chi2_data_per_row"] = "1" };
        var fullResult = Compare([full], [new Row(full) { ["chi2_data_per_row"] = "2" }], "a", "b", "full", false, "full").Summary;
        Require(!fullResult.ContainsKey("p_two_sided") && full["n_test"] == "0", "full fits descriptive and inputs immutable");
        Require(Math.Abs(NormalMagnitude(.05) - 1.959963984540054) < 1e-13, "normal-equivalent magnitude");
        return new Result { ["passed"] = true, ["checks"] = checks };
    }
}

public static partial class Statistics
{
    static readonly Dictionary<string, string> TierIds = new() { ["fixed_inputs"] = "A", ["matched_nuisance"] = "B", ["expanded_family"] = "C" };
    static readonly string[] Morphologies = ["S0", "Sa", "Sab", "Sb", "Sbc", "Sc", "Scd", "Sd", "Sdm", "Sm", "Im", "BCD"];
    static readonly (string Left, string Right)[] ClassPairs = [("A_core", "A_nfw"), ("B_core", "B_nfw"), ("C_core", "C_nfw"), ("C_core", "B_core")];
    static readonly (string Left, string Right)[] BaryonPairs = [("A_core", "bar_fixed"), ("A_nfw", "bar_fixed"), ("B_core", "bar_profiled"), ("C_core", "bar_profiled"), ("B_nfw", "bar_profiled")];

    static Result AdditionalAnalyses(string numerics, string outputRoot)
    {
        var predictions = Csv.Read(Path.Combine(numerics, "fit_tiers/tier_predictions.csv"));
        var fits = Csv.Read(Path.Combine(numerics, "fit_tiers/tier_fits.csv"));
        static (string Tier, string Model, string Mode, string Galaxy) Key(Row r) => (Text(r, "tier"), Text(r, "model"), Text(r, "fit_mode"), Text(r, "galaxy"));
        var fitByKey = fits.ToDictionary(Key);
        var groups = predictions.GroupBy(Key).OrderBy(g => g.Key.Tier, StringComparer.Ordinal).ThenBy(g => g.Key.Model, StringComparer.Ordinal).ThenBy(g => g.Key.Mode, StringComparer.Ordinal).ThenBy(g => g.Key.Galaxy, StringComparer.Ordinal).ToList();
        if (!groups.Select(g => g.Key).ToHashSet().SetEquals(fitByKey.Keys))
            throw new InvalidDataException("Tier prediction and fit keys differ.");
        var reference = groups.Where(g => g.Key.Tier == "fixed_inputs" && g.Key.Model == "core").ToDictionary(g => (g.Key.Mode, g.Key.Galaxy), g => g.ToList());
        var classScores = new List<Result>();
        var regions = new List<Result>();
        foreach (var group in groups)
        {
            var points = group.OrderBy(r => Number(r, "radius_kpc")).ToList();
            if (points.Select(r => Number(r, "radius_kpc")).Distinct().Count() != points.Count)
                throw new InvalidDataException("Duplicate tier prediction radius.");
            MatchPoints(points, reference[(group.Key.Mode, group.Key.Galaxy)]);
            var fit = fitByKey[group.Key];
            if (Count(fit, "n_rows") != points.Count || Count(fit, "n_train") != points.Count(r => Flag(r, "is_training")) || Count(fit, "n_outer") != points.Count(r => Flag(r, "is_outer")))
                throw new InvalidDataException("Tier prediction partition count mismatch.");
            bool full = group.Key.Mode == "full";
            var picked = full ? points : points.Where(r => Flag(r, "is_outer")).ToList();
            var score = Score(picked, Flag(fit, full ? "train_admissible" : "outer_admissible"));
            var row = new Result { ["galaxy"] = group.Key.Galaxy, ["model_id"] = TierIds[group.Key.Tier] + "_" + group.Key.Model, ["fit_mode"] = group.Key.Mode, ["region"] = full ? "all" : "outer" };
            foreach (var item in score)
                row[item.Key] = item.Value;
            foreach (string field in new[] { "halo_bound_hit", "eta_bound_hit", "nuisance_bound_hit", "guard_active", "optimizer_success" })
                row[field] = Flag(fit, field);
            row["eta"] = Number(fit, "eta");
            classScores.Add(row);
            if (full)
            {
                int offset = 0;
                for (int third = 0; third < 3; third++)
                {
                    int length = points.Count / 3 + (third < points.Count % 3 ? 1 : 0);
                    var region = new Result { ["galaxy"] = group.Key.Galaxy, ["model_id"] = row["model_id"], ["fit_mode"] = "full", ["region"] = new[] { "inner_third", "middle_third", "outer_third" }[third] };
                    foreach (var item in Score(points.Skip(offset).Take(length).ToList(), Flag(fit, "train_admissible")))
                        region[item.Key] = item.Value;
                    regions.Add(region);
                    offset += length;
                }
            }
            else
            {
                var region = new Result { ["galaxy"] = group.Key.Galaxy, ["model_id"] = row["model_id"], ["fit_mode"] = "inner", ["region"] = "withheld_outer" };
                foreach (var item in score)
                    region[item.Key] = item.Value;
                regions.Add(region);
            }
        }
        var bNfw = predictions.Where(r => Text(r, "tier") == "matched_nuisance" && Text(r, "model") == "nfw").OrderBy(r => Text(r, "fit_mode"), StringComparer.Ordinal).ThenBy(r => Text(r, "galaxy"), StringComparer.Ordinal).ThenBy(r => Number(r, "radius_kpc")).ToList();
        var cNfw = predictions.Where(r => Text(r, "tier") == "expanded_family" && Text(r, "model") == "nfw").OrderBy(r => Text(r, "fit_mode"), StringComparer.Ordinal).ThenBy(r => Text(r, "galaxy"), StringComparer.Ordinal).ThenBy(r => Number(r, "radius_kpc")).ToList();
        if (bNfw.Count != cNfw.Count || bNfw.Where((r, i) => r.Any(v => v.Key != "tier" && v.Value != cNfw[i][v.Key])).Any())
            throw new InvalidDataException("B and C NFW predictions differ.");
        string classOut = Path.Combine(outputRoot, "galaxy_classes");
        Csv.Write(Path.Combine(classOut, "per_galaxy_scores.csv"), classScores);
        Csv.Write(Path.Combine(classOut, "per_galaxy_regions.csv"), regions);
        var validation = new Result
        {
            ["class_score_reconstruction"] = ValidateCsv(Path.Combine(numerics, "galaxy_classes/per_galaxy_scores.csv"), Path.Combine(classOut, "per_galaxy_scores.csv")),
            ["class_region_reconstruction"] = ValidateCsv(Path.Combine(numerics, "galaxy_classes/per_galaxy_regions.csv"), Path.Combine(classOut, "per_galaxy_regions.csv"))
        };
        var morphology = new Dictionary<string, int>();
        foreach (string line in File.ReadLines(Path.Combine(numerics, "observations/data/SPARC_Lelli2016c.mrt")))
        {
            string[] p = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (p.Length == 19 && int.TryParse(p[1], out int type) && double.TryParse(p[2], NumberStyles.Float, CultureInfo.InvariantCulture, out _))
                morphology.Add(p[0], type);
        }
        if (morphology.Count != 175 || morphology.Values.Any(v => v < 0 || v > 11))
            throw new InvalidDataException("Invalid SPARC morphology catalogue.");
        var comparisons = new List<Result>();
        var classSummary = new List<Result>();
        var residualSummary = new List<Result>();
        foreach (bool broad in new[] { false, true })
            for (int classId = 0; classId < (broad ? 3 : 12); classId++)
            {
                var label = new Result { ["classification"] = broad ? "broad" : "detailed", ["class_id"] = classId, ["class_name"] = broad ? new[] { "S0", "spirals", "Im/BCD" }[classId] : Morphologies[classId], ["class_order"] = classId };
                var names = morphology.Where(kv => (broad ? kv.Value == 0 ? 0 : kv.Value <= 9 ? 1 : 2 : kv.Value) == classId).Select(kv => kv.Key).ToHashSet();
                var scores = classScores.Where(r => names.Contains((string)r["galaxy"]!)).ToList();
                var regional = regions.Where(r => names.Contains((string)r["galaxy"]!)).ToList();
                foreach (var (left, right) in ClassPairs)
                {
                    var row = new Result(label) { ["pair_id"] = left + "_vs_" + right, ["left_model"] = left, ["right_model"] = right };
                    var pair = Align(scores, left, right, "inner", "outer");
                    foreach (var item in ClassPaired(pair.Left, pair.Right))
                        row[item.Key] = item.Value;
                    comparisons.Add(row);
                }
                foreach (string model in new[] { "A_core", "A_nfw", "B_core", "B_nfw", "C_core", "C_nfw" })
                {
                    foreach (string mode in new[] { "full", "inner" })
                    {
                        var subset = scores.Where(r => (string)r["model_id"]! == model && (string)r["fit_mode"]! == mode).ToList();
                        classSummary.Add(new Result(label) { ["model_id"] = model, ["fit_mode"] = mode, ["region"] = mode == "full" ? "all" : "outer", ["n_galaxies"] = subset.Count, ["n_failures"] = subset.Count(r => (bool)r["failed"]!), ["median_loss"] = Median(subset.Select(r => (double)r["loss"]!)), ["mean_loss"] = Mean(subset.Select(r => (double)r["loss"]!)), ["median_rms_kms"] = Median(subset.Select(r => (double)r["rms_kms"]!)), ["mean_rms_kms"] = Mean(subset.Select(r => (double)r["rms_kms"]!)) });
                    }
                    foreach (string region in new[] { "inner_third", "middle_third", "outer_third", "withheld_outer" })
                    {
                        var subset = regional.Where(r => (string)r["model_id"]! == model && (string)r["region"]! == region).ToList();
                        var valid = subset.Where(r => !(bool)r["failed"]!).ToList();
                        residualSummary.Add(new Result(label) { ["model_id"] = model, ["fit_mode"] = region == "withheld_outer" ? "inner" : "full", ["region"] = region, ["n_galaxies"] = subset.Count, ["n_failures"] = subset.Count - valid.Count, ["n_signed_defined"] = valid.Count, ["median_signed_standardized"] = Median(valid.Select(r => (double)r["signed_standardized"]!)), ["mean_signed_standardized"] = Mean(valid.Select(r => (double)r["signed_standardized"]!)), ["median_signed_kms"] = Median(valid.Select(r => (double)r["signed_kms"]!)), ["mean_signed_kms"] = Mean(valid.Select(r => (double)r["signed_kms"]!)), ["median_loss"] = Median(subset.Select(r => (double)r["loss"]!)), ["mean_loss"] = Mean(subset.Select(r => (double)r["loss"]!)) });
                    }
                }
            }
        if (comparisons.Count != 60)
            throw new InvalidDataException("Class Holm family must contain 60 comparisons.");
        var adjusted = Holm(comparisons.Select(r => r["p_raw"] is double p ? p : 1).ToArray());
        for (int i = 0; i < comparisons.Count; i++)
        {
            bool informative = (int)comparisons[i]["n_informative"]! > 0;
            comparisons[i]["holm_family_size"] = 60;
            comparisons[i]["p_holm"] = informative ? adjusted[i] : null;
            comparisons[i]["z_holm"] = informative ? NormalMagnitude(adjusted[i]) : null;
        }
        foreach (var (filename, rows) in new[] { ("class_comparisons.csv", comparisons), ("class_model_scores.csv", classSummary), ("residual_regions.csv", residualSummary) })
        {
            Csv.Write(Path.Combine(classOut, filename), rows);
            validation[filename] = ValidateCsv(Path.Combine(numerics, "galaxy_classes", filename), Path.Combine(classOut, filename));
        }
        var baseline = Csv.Read(Path.Combine(numerics, "baryons_only/baseline_predictions.csv"));
        foreach (var group in baseline.GroupBy(r => (Mode: Text(r, "fit_mode"), Galaxy: Text(r, "galaxy"), Baseline: Text(r, "baseline"))))
            MatchPoints(group.ToList(), reference[(group.Key.Mode, group.Key.Galaxy)]);
        var allPoints = predictions.Where(r => !(Text(r, "tier") == "expanded_family" && Text(r, "model") == "nfw")).Select(r => new Row(r) { ["model_id"] = TierIds[Text(r, "tier")] + "_" + Text(r, "model"), ["pipeline_failure"] = "False" })
            .Concat(baseline.Select(r => new Row(r) { ["model_id"] = "bar_" + Text(r, "baseline") }));
        var baryonScores = new List<Result>();
        foreach (var group in allPoints.GroupBy(r => (Model: Text(r, "model_id"), Mode: Text(r, "fit_mode"), Galaxy: Text(r, "galaxy"))).OrderBy(g => g.Key.Model, StringComparer.Ordinal).ThenBy(g => g.Key.Mode, StringComparer.Ordinal).ThenBy(g => g.Key.Galaxy, StringComparer.Ordinal))
            foreach (string region in group.Key.Mode == "full" ? new[] { "all" } : new[] { "training", "outer" })
            {
                var points = group.Where(r => region == "all" || Flag(r, region == "training" ? "is_training" : "is_outer")).ToList();
                var score = Score(points, true, false);
                baryonScores.Add(new Result { ["model_id"] = group.Key.Model, ["fit_mode"] = group.Key.Mode, ["region"] = region, ["galaxy"] = group.Key.Galaxy, ["n_rows"] = points.Count, ["failed"] = score["failed"], ["pipeline_failure"] = group.Any(r => Flag(r, "pipeline_failure")), ["invalid_rows"] = points.Count(r => !double.IsFinite(Number(r, "prediction_kms"))), ["loss"] = score["loss"], ["rms_kms"] = score["rms_kms"], ["signed_mean_residual"] = score["signed_kms"] });
            }
        string baryonOut = Path.Combine(outputRoot, "baryons_only");
        Csv.Write(Path.Combine(baryonOut, "galaxy_scores.csv"), baryonScores);
        validation["baryon_score_reconstruction"] = ValidateCsv(Path.Combine(numerics, "baryons_only/galaxy_scores.csv"), Path.Combine(baryonOut, "galaxy_scores.csv"));
        var tests = new List<Result>();
        var paired = new List<Result>();
        foreach (string mode in new[] { "full", "inner" })
            foreach (var (left, right) in BaryonPairs)
            {
                string region = mode == "full" ? "all" : "outer";
                var pair = Align(baryonScores, left, right, mode, region);
                if (pair.Left.Count != (mode == "full" ? 153 : 131))
                    throw new InvalidDataException("Incomplete baryon comparison sample.");
                for (int i = 0; i < pair.Left.Count; i++)
                {
                    var a = pair.Left[i];
                    var b = pair.Right[i];
                    paired.Add(new Result { ["augmented"] = left, ["baseline"] = right, ["fit_mode"] = mode, ["region"] = region, ["galaxy"] = a["galaxy"], ["augmented_loss"] = a["loss"], ["baseline_loss"] = b["loss"], ["augmented_rms_kms"] = a["rms_kms"], ["baseline_rms_kms"] = b["rms_kms"] });
                }
                if (mode == "inner")
                {
                    var row = new Result { ["augmented"] = left, ["baseline"] = right };
                    foreach (var item in BaryonSign(pair.Left.Select(r => (double)r["loss"]!).ToArray(), pair.Right.Select(r => (double)r["loss"]!).ToArray()))
                        row[item.Key] = item.Value;
                    tests.Add(row);
                }
            }
        adjusted = Holm(tests.Select(r => (double)r["p_raw"]!).ToArray());
        for (int i = 0; i < tests.Count; i++)
        {
            tests[i]["p_holm"] = adjusted[i];
            tests[i]["z_holm"] = NormalMagnitude(adjusted[i]);
            tests[i]["family_size"] = 5;
        }
        var modelSummary = new List<Result>();
        foreach (var group in baryonScores.Where(r => (string)r["fit_mode"]! == "full" && (string)r["region"]! == "all" || (string)r["fit_mode"]! == "inner" && (string)r["region"]! == "outer").GroupBy(r => (Model: (string)r["model_id"]!, Mode: (string)r["fit_mode"]!, Region: (string)r["region"]!)))
            foreach (string sample in new[] { "all_selected", "finite_predictions" })
            {
                var subset = group.Where(r => sample == "all_selected" || !(bool)r["failed"]!).ToList();
                int n = subset.Sum(r => (int)r["n_rows"]!);
                modelSummary.Add(new Result { ["model_id"] = group.Key.Model, ["fit_mode"] = group.Key.Mode, ["region"] = group.Key.Region, ["sample"] = sample, ["n_total"] = group.Count(), ["n_galaxies"] = subset.Count, ["n_failures"] = group.Count(r => (bool)r["failed"]!), ["n_rows"] = n, ["median_loss"] = Median(subset.Select(r => (double)r["loss"]!)), ["mean_loss"] = Mean(subset.Select(r => (double)r["loss"]!)), ["median_rms_kms"] = Median(subset.Select(r => (double)r["rms_kms"]!)), ["mean_rms_kms"] = Mean(subset.Select(r => (double)r["rms_kms"]!)), ["pooled_loss"] = subset.Sum(r => (double)r["loss"]! * (int)r["n_rows"]!) / n });
            }
        foreach (var (filename, rows) in new[] { ("sign_tests.csv", tests), ("paired_galaxy_scores.csv", paired), ("model_summary.csv", modelSummary) })
        {
            Csv.Write(Path.Combine(baryonOut, filename), rows);
            validation[filename] = ValidateCsv(Path.Combine(numerics, "baryons_only", filename), Path.Combine(baryonOut, filename));
        }
        validation["passed"] = true;
        Console.WriteLine("Statistics: reconstructed 60 morphology comparisons and 5 baryons-only comparisons from prediction rows.");
        return validation;
    }

    static void MatchPoints(List<Row> left, List<Row> right)
    {
        var a = left.OrderBy(r => Number(r, "radius_kpc")).ToList();
        var b = right.OrderBy(r => Number(r, "radius_kpc")).ToList();
        if (a.Count != b.Count)
            throw new InvalidDataException("Paired point count differs.");
        for (int i = 0; i < a.Count; i++)
        {
            foreach (string key in new[] { "radius_kpc", "observed_kms", "error_kms", "catalogue_baryon_v2" })
                if (Math.Abs(Number(a[i], key) - Number(b[i], key)) > 2e-12 + 2e-14 * Math.Abs(Number(b[i], key)))
                    throw new InvalidDataException("Paired point observations differ: " + key);
            foreach (string key in new[] { "is_training", "is_outer", "eligible_holdout" })
                if (Flag(a[i], key) != Flag(b[i], key))
                    throw new InvalidDataException("Paired point partition differs: " + key);
        }
    }
    static Result Score(List<Row> rows, bool admissible, bool rejectNegative = true)
    {
        if (rows.Count == 0)
            throw new InvalidDataException("Empty scored region.");
        double chiSquared = 0, squaredResidualSum = 0, residualSum = 0, standardizedSum = 0;
        bool valid = admissible;
        foreach (var row in rows)
        {
            double observed = Number(row, "observed_kms");
            double error = Number(row, "error_kms");
            double prediction = Number(row, "prediction_kms");
            if (!double.IsFinite(observed) || !double.IsFinite(error) || error <= 0)
                throw new InvalidDataException("Malformed observations/errors.");
            valid &= double.IsFinite(prediction) && (!rejectNegative || prediction >= 0);
            double residual = prediction - observed;
            double standardizedResidual = residual / error;
            chiSquared += standardizedResidual * standardizedResidual;
            squaredResidualSum += residual * residual;
            residualSum += residual;
            standardizedSum += standardizedResidual;
        }
        double rms = Math.Sqrt(squaredResidualSum / rows.Count);
        valid &= double.IsFinite(chiSquared) && double.IsFinite(rms);
        // Failed fits retain infinite loss and remain in the paired galaxy denominator.
        return new Result
        {
            ["n_rows"] = rows.Count,
            ["failed"] = !valid,
            ["chi2"] = valid ? chiSquared : double.PositiveInfinity,
            ["loss"] = valid ? chiSquared / rows.Count : double.PositiveInfinity,
            ["rms_kms"] = valid ? rms : double.PositiveInfinity,
            ["signed_kms"] = valid ? residualSum / rows.Count : double.NaN,
            ["signed_standardized"] = valid ? standardizedSum / rows.Count : double.NaN
        };
    }

    static (List<Result> Left, List<Result> Right) Align(List<Result> rows, string left, string right, string mode, string region)
    {
        List<Result> Pick(string model) => rows.Where(r => (string)r["model_id"]! == model && (string)r["fit_mode"]! == mode && (string)r["region"]! == region).OrderBy(r => (string)r["galaxy"]!, StringComparer.Ordinal).ToList();
        var leftRows = Pick(left);
        var rightRows = Pick(right);
        if (leftRows.Count != rightRows.Count
            || leftRows.Select(row => row["galaxy"]).Distinct().Count() != leftRows.Count
            || rightRows.Select(row => row["galaxy"]).Distinct().Count() != rightRows.Count
            || leftRows.Where((row, index) => !Equals(row["galaxy"], rightRows[index]["galaxy"])
                || !Equals(row["n_rows"], rightRows[index]["n_rows"])).Any())
            throw new InvalidDataException("Unmatched paired observations.");
        return (leftRows, rightRows);
    }
    static Result ClassPaired(List<Result> left, List<Result> right)
    {
        double[] x = left.Select(r => (bool)r["failed"]! ? double.PositiveInfinity : (double)r["loss"]!).ToArray(), y = right.Select(r => (bool)r["failed"]! ? double.PositiveInfinity : (double)r["loss"]!).ToArray();
        var common = x.Zip(y, (a, b) => (a, b)).ToArray();
        int wins = common.Count(v => v.a < v.b), losses = common.Count(v => v.a > v.b);
        var inference = ExactSign(wins, losses);
        var ci = inference["left_win_fraction_ci95"] as double[];
        double[] delta = common.Where(v => !(double.IsPositiveInfinity(v.a) && double.IsPositiveInfinity(v.b))).Select(v => v.a - v.b).ToArray();
        var interval = MedianInterval(delta);
        var dci = interval["ci95"] as double[];
        double[] ratio = common.Where(v => double.IsFinite(v.a) && double.IsFinite(v.b) && v.a > 0 && v.b > 0).Select(v => Math.Log10(v.a) - Math.Log10(v.b)).ToArray();
        return new Result { ["n_total"] = left.Count, ["n_informative"] = wins + losses, ["left_wins"] = wins, ["right_wins"] = losses, ["finite_ties"] = common.Count(v => double.IsFinite(v.a) && v.a == v.b), ["left_failures"] = x.Count(v => !double.IsFinite(v)), ["right_failures"] = y.Count(v => !double.IsFinite(v)), ["joint_failures"] = common.Count(v => !double.IsFinite(v.a) && !double.IsFinite(v.b)), ["direction"] = wins > losses ? "left" : losses > wins ? "right" : "none", ["left_win_fraction"] = inference["left_win_fraction"], ["win_ci_low"] = ci?[0], ["win_ci_high"] = ci?[1], ["p_raw"] = inference["p_two_sided"], ["z_raw"] = inference["z_two_sided"], ["median_delta"] = Median(delta), ["mean_delta"] = Mean(delta), ["delta_ci_low"] = dci?[0], ["delta_ci_high"] = dci?[1], ["delta_interval_n"] = interval["n"], ["delta_interval_k"] = interval["k_1based"], ["delta_interval_coverage"] = interval["coverage_continuous_null"], ["median_log10_ratio"] = Median(ratio), ["n_positive_finite_pairs"] = ratio.Length };
    }
    public static Result BaryonSign(double[] x, double[] y)
    {
        if (x.Length != y.Length)
            throw new InvalidDataException("Unequal paired arrays.");
        var pairs = x.Zip(y, (a, b) => (Augmented: a, Baseline: b)).ToArray();
        // A finite fit beats a failed fit; joint failures are neither wins nor finite ties.
        int wins = pairs.Count(pair => double.IsFinite(pair.Augmented)
            && (!double.IsFinite(pair.Baseline) || pair.Augmented < pair.Baseline));
        int losses = pairs.Count(pair => double.IsFinite(pair.Baseline)
            && (!double.IsFinite(pair.Augmented) || pair.Baseline < pair.Augmented));
        var exact = ExactSign(wins, losses);
        var interval = exact["left_win_fraction_ci95"] as double[];
        return new Result
        {
            ["n_total"] = x.Length,
            ["n_informative"] = wins + losses,
            ["augmented_wins"] = wins,
            ["baseline_wins"] = losses,
            ["finite_ties"] = pairs.Count(pair => double.IsFinite(pair.Augmented) && double.IsFinite(pair.Baseline) && pair.Augmented == pair.Baseline),
            ["augmented_failures"] = x.Count(value => !double.IsFinite(value)),
            ["baseline_failures"] = y.Count(value => !double.IsFinite(value)),
            ["joint_failures"] = pairs.Count(pair => !double.IsFinite(pair.Augmented) && !double.IsFinite(pair.Baseline)),
            ["direction"] = wins > losses ? "augmented" : losses > wins ? "baseline" : "tie",
            ["win_fraction"] = exact["left_win_fraction"],
            ["win_ci_low"] = interval?[0],
            ["win_ci_high"] = interval?[1],
            ["p_raw"] = exact["p_two_sided"],
            ["z_raw"] = exact["z_two_sided"]
        };
    }

    static double Mean(IEnumerable<double> values)
    {
        var sample = values.ToArray();
        return sample.Length == 0 ? double.NaN : sample.Average();
    }
}
