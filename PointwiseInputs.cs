using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

public static class PointwiseInputs
{
    static readonly string[] Models = ["extended", "scalar_core", "compact_plummer", "nfw", "baryon_only", "profiled_baryons", "adaptive_primary"];
    public static Dictionary<string, object?> Verify(string dataRoot, string outputRoot)
    {
        string folder = Path.Combine(dataRoot, "numerics/pointwise_comparison_2026_09_25/data");
        var allRows = Csv.Read(Path.Combine(folder, "pointwise_all.csv"));
        var outerRows = Csv.Read(Path.Combine(folder, "pointwise_outer.csv"));
        var currentPredictions = Csv.Read(Path.Combine(folder, "inputs/current_predictions.csv"))
            .Where(row => Models.Contains(row["model"])).GroupBy(row => (row["galaxy"], row["model"]))
            .ToDictionary(g => g.Key, g => g.OrderBy(row => Csv.Number(row, "r_catalogue_kpc")).ToArray());
        var profiledPredictions = Csv.Read(Path.Combine(folder, "inputs/profiled_baseline_predictions.csv"))
            .Where(row => row["baseline"] == "profiled" && row["fit_mode"] == "inner").GroupBy(row => row["galaxy"])
            .ToDictionary(g => g.Key, g => g.OrderBy(row => Csv.Number(row, "radius_kpc")).ToArray());
        var profiledFits = Csv.Read(Path.Combine(folder, "inputs/profiled_baseline_fits.csv"))
            .Where(row => row["baseline"] == "profiled" && row["fit_mode"] == "inner").ToDictionary(row => row["galaxy"]);
        var sparc = Sparc.Load(dataRoot);
        var eligibleGalaxies = sparc.SelectedNames.Where(g => sparc.Points[g].Length >= 8).ToHashSet(StringComparer.Ordinal);
        int checks = 0;
        double maximumDifference = 0;
        void Require(bool pass, string label)
        {
            checks++;
            if (!pass)
                throw new InvalidDataException("Pointwise input mismatch: " + label);
        }
        void Close(double a, double b, string label)
        {
            checks++;
            if (double.IsNaN(a) && double.IsNaN(b) || a == b)
                return;
            double error = Math.Abs(a - b);
            maximumDifference = Math.Max(maximumDifference, error);
            if (!double.IsFinite(error) || error > 1e-12 + 2e-12 * Math.Abs(b))
                throw new InvalidDataException($"Pointwise input {label}: {a:R} != {b:R}");
        }
        bool Flag(Row row, string key)
        {
            string value = row[key].ToLowerInvariant();
            if (value != "true" && value != "false")
                throw new InvalidDataException("Invalid pointwise boolean: " + key);
            return value == "true";
        }
        var predictionsByKey = new Dictionary<(string, int, string), Row>();
        var numericFields = new[] { "radius_kpc", "observed_kms", "sigma_kms", "training_edge_kpc", "predicted_kms", "inclination_factor", "distance_ratio", "residual_kms", "standardized_residual", "squared_standardized_error" };
        Require(eligibleGalaxies.Count == 131, "eligible galaxy count");
        Require(allRows.Count == 3034 * 7 && outerRows.Count == 659 * 7, "complete model-point rows");
        // Compare every saved prediction and split with its SPARC observation and source fit.
        foreach (var row in allRows)
        {
            string galaxy = row["galaxy"], model = row["model"];
            double indexValue = Csv.Number(row, "point_index");
            int index = (int)indexValue;
            Require(indexValue == index && index >= 0, "integer point index");
            Require(eligibleGalaxies.Contains(galaxy) && Models.Contains(model), "eligible galaxy/model");
            var points = sparc.Points[galaxy];
            Require(index < points.Length, "point index range");
            Require(predictionsByKey.TryAdd((galaxy, index, model), row), "duplicate model-point key");
            var observation = points[index];
            int trainingCount = points.Length - Sparc.HoldoutCount(points.Length);
            Close(Csv.Number(row, "radius_kpc"), observation.RadiusKpc, "raw SPARC radius");
            Close(Csv.Number(row, "observed_kms"), observation.ObservedSpeedKms, "raw SPARC speed");
            Close(Csv.Number(row, "sigma_kms"), observation.ErrorSpeedKms, "raw SPARC error");
            Require(Csv.Number(row, "sigma_kms") > 0, "positive quoted error");
            Require(Flag(row, "is_outer") == (index >= trainingCount) && Flag(row, "is_training") == (index < trainingCount), "original outer split");
            Close(Csv.Number(row, "training_edge_kpc"), points[trainingCount - 1].RadiusKpc, "training edge");
            double predicted = Csv.Number(row, "predicted_kms"), sigma = Csv.Number(row, "sigma_kms"), residual = predicted - Csv.Number(row, "observed_kms");
            Close(Csv.Number(row, "residual_kms"), residual, "velocity residual");
            Close(Csv.Number(row, "standardized_residual"), residual / sigma, "standardized residual");
            Close(Csv.Number(row, "squared_standardized_error"), Math.Pow(residual / sigma, 2), "squared score");
            bool valid = Flag(row, "valid_prediction"), failed = Flag(row, "fit_failed");
            Require(valid == (double.IsFinite(predicted) && !failed), "validity and failed-fit state");
            if (model == "profiled_baryons")
            {
                var original = profiledPredictions[galaxy][index];
                var fit = profiledFits[galaxy];
                Close(predicted, Csv.Number(original, "prediction_kms"), "profiled prediction snapshot");
                Require(failed == Csv.Flag(original, "pipeline_failure"), "profiled training failure");
                Require(row["failure_reason"] == (valid ? "" : original["fit_status"]), "profiled failure reason");
                Require(row["conditioning"] == "independently_profiled_baryons_inner_only_matched_priors", "profiled conditioning");
                Close(Csv.Number(row, "inclination_factor"), Csv.Number(fit, "inclination_factor"), "profiled inclination");
                Close(Csv.Number(row, "distance_ratio"), Csv.Number(fit, "distance_ratio"), "profiled distance");
            }
            else
            {
                var original = currentPredictions[(galaxy, model)][index];
                Close(predicted, Csv.Number(original, "predicted_catalogue_kms"), "current prediction snapshot");
                Close(Csv.Number(original, "r_catalogue_kpc"), observation.RadiusKpc, "snapshot radius");
                Close(Csv.Number(original, "observed_catalogue_kms"), observation.ObservedSpeedKms, "snapshot observed speed");
                Close(Csv.Number(original, "error_catalogue_kms"), observation.ErrorSpeedKms, "snapshot error");
                Require(Flag(row, "is_outer") == Csv.Flag(original, "is_outer") && Flag(row, "is_training") == Csv.Flag(original, "is_training"), "snapshot split");
                Require(!failed, "finite-primary fit status");
                Require(row["conditioning"] == original["conditioning"], "prediction conditioning");
                Require(valid == (double.IsFinite(predicted) && !Csv.Flag(original, "invalid_prediction")), "declared prediction validity");
                Close(Csv.Number(row, "inclination_factor"), Csv.Number(original, "inclination_factor"), "prediction inclination");
                Close(Csv.Number(row, "distance_ratio"), Csv.Number(original, "distance_ratio"), "prediction distance");
            }
        }
        Require(allRows.Select(row => row["galaxy"]).ToHashSet().SetEquals(eligibleGalaxies), "complete galaxy membership");
        foreach (string galaxy in eligibleGalaxies)
            foreach (string model in Models)
                Require(allRows.Count(row => row["galaxy"] == galaxy && row["model"] == model) == sparc.Points[galaxy].Length, "no missing model radii");

        // The outer file must be the exact withheld subset, including failed predictions.
        var outerKeys = new HashSet<(string, int, string)>();
        foreach (var row in outerRows)
        {
            double indexValue = Csv.Number(row, "point_index");
            int index = (int)indexValue;
            Require(indexValue == index && index >= 0, "integer outer point index");
            var key = (row["galaxy"], index, row["model"]);
            Require(outerKeys.Add(key), "duplicate outer row");
            Require(Flag(row, "is_outer") && !Flag(row, "is_training"), "only withheld points");
            Require(predictionsByKey.TryGetValue(key, out var original), "outer row present in full input");
            foreach (string field in numericFields)
                Close(Csv.Number(row, field), Csv.Number(original!, field), "outer/full " + field);
            foreach (string field in new[] { "valid_prediction", "fit_failed", "failure_reason", "conditioning" })
                Require(row[field] == original![field], "outer/full " + field);
        }
        Require(outerKeys.SetEquals(predictionsByKey.Where(kv => Flag(kv.Value, "is_outer")).Select(kv => kv.Key)), "complete outer block");
        var failedOuter = outerRows.Where(row => Flag(row, "fit_failed")).ToArray();
        Require(failedOuter.Length == 5 && failedOuter.All(row => row["galaxy"] == "UGC01281" && row["model"] == "profiled_baryons" && !Flag(row, "valid_prediction") && double.IsNaN(Csv.Number(row, "predicted_kms"))), "retained profiled failure");
        Require(outerRows.Where(row => row["model"] != "profiled_baryons").All(row => Flag(row, "valid_prediction")), "finite outer comparison predictions");
        int snapshotHashes = 0;
        foreach (var item in Data.Json(Path.Combine(folder, "input_manifest.json")).EnumerateArray())
        {
            string name = item.GetProperty("snapshot").GetString()!;
            if (!name.EndsWith(".csv", StringComparison.Ordinal))
                continue;
            string digest = Data.FileSha256(Path.Combine(folder, "inputs", name));
            Require(digest == item.GetProperty("sha256").GetString(), "immutable source snapshot " + name);
            snapshotHashes++;
        }
        var report = new Dictionary<string, object?>
        {
            ["passed"] = true,
            ["checks"] = checks,
            ["snapshot_hashes"] = snapshotHashes,
            ["models"] = Models,
            ["galaxies"] = eligibleGalaxies.Count,
            ["all_measurements"] = allRows.Count / 7,
            ["outer_measurements"] = outerRows.Count / 7,
            ["outer_model_rows"] = outerRows.Count,
            ["maximum_absolute_difference"] = maximumDifference,
            ["raw_sparc_coordinates_and_splits_match"] = true,
            ["predictions_match_frozen_sources"] = true,
            ["failed_profiled_galaxy"] = "UGC01281",
            ["failed_outer_measurements"] = 5,
            ["fitted_parameters_or_observations_changed"] = false
        };
        Data.SaveJson(Path.Combine(outputRoot, "pointwise_comparison/input_validation.json"), report);
        return report;
    }
}
