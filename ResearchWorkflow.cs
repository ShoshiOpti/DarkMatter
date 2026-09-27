using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;

namespace DarkUniverse;

/// <summary>Locked development protocols and auditable fitting/prediction contracts, separate from publication replay.</summary>
public static class ResearchWorkflow
{
    static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true
    };

    static string Hash(string path) => Data.FileSha256(path);

    static void Check(bool condition, string message)
    {
        if (!condition) throw new InvalidDataException(message);
    }

    static JsonObject Read(string path) => JsonNode.Parse(File.ReadAllText(path))!.AsObject();

    static string[] Strings(JsonNode? node) => node!.AsArray().Select(x => x!.GetValue<string>()).ToArray();

    static string Key(IReadOnlyDictionary<string, string> row) => row["galaxy"] + ":" + row["index"];

    static void Save(string path, object value)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(path))!);
        File.WriteAllText(path, JsonSerializer.Serialize(value, JsonOptions));
    }

    static string Local(string root, string name)
    {
        string full = Path.GetFullPath(Path.Combine(root, name));
        Check(!Path.IsPathRooted(name) && Inside(root, full), "Protocol paths must remain local: " + name);
        Check(ExecutionSafety.CanonicalizePath(full).Equals(full, OperatingSystem.IsWindows() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal), "Research inputs cannot contain filesystem links.");
        if (File.Exists(full)) ExecutionSafety.AssertUnaliasedInputFile(full);
        return full;
    }

    static bool Inside(string root, string path)
    {
        string rel = Path.GetRelativePath(Path.GetFullPath(root), Path.GetFullPath(path));
        return rel == "." || (!Path.IsPathRooted(rel) && rel != ".." && !rel.StartsWith(".." + Path.DirectorySeparatorChar));
    }

    public static object CreateProtocol(string dataRoot, string outputRoot)
    {
        ExecutionSafety.ValidateRoots(dataRoot, outputRoot);
        Directory.CreateDirectory(outputRoot);
        Check(!File.Exists(Path.Combine(outputRoot, "protocol.json")), "Protocol already exists; use a fresh output directory.");

        var data = Sparc.Load(dataRoot);
        var galaxies = data.SelectedNames.Where(n => data.Points[n].Length >= 8).ToArray();

        var predictors = new List<Dictionary<string, object?>>();
        var train = new List<Dictionary<string, object?>>();
        var evaluation = new List<Dictionary<string, object?>>();
        var split = new List<Dictionary<string, object?>>();
        var metadata = new Dictionary<string, object?>();

        for (int g = 0; g < galaxies.Length; g++)
        {
            string name = galaxies[g];
            var m = data.Catalogue[name];
            var rows = data.Points[name];
            int nt = rows.Length - Math.Max(2, (int)Math.Ceiling(.2 * rows.Length));

            metadata[name] = new
            {
                D = m.DistanceMpc,
                e_D = m.DistanceErrorMpc,
                inc = m.InclinationDeg,
                e_inc = m.InclinationErrorDeg,
                Rdisk = m.DiskScaleKpc,
                distance_method = m.DistanceMethod,
                catalogue_D = m.DistanceMpc,
                catalogue_inc = m.InclinationDeg
            };

            for (int i = 0; i < rows.Length; i++)
            {
                var p = rows[i];
                string role = i < nt ? "inner" : "outer", population = g % 5 == 0 ? "transfer" : "train";

                predictors.Add(new()
                {
                    ["galaxy"] = name,
                    ["index"] = i,
                    ["radius_kpc"] = p.RadiusKpc,
                    ["sigma_kms"] = p.ErrorSpeedKms,
                    ["gas_kms"] = p.GasSpeedKms,
                    ["disk_kms"] = p.DiskSpeedKms,
                    ["bulge_kms"] = p.BulgeSpeedKms,
                    ["partition"] = role,
                    ["population_role"] = population
                });

                var obs = new Dictionary<string, object?>()
                {
                    ["galaxy"] = name,
                    ["index"] = i,
                    ["observed_kms"] = p.ObservedSpeedKms,
                    ["sigma_kms"] = p.ErrorSpeedKms
                };

                (i < nt ? train : evaluation).Add(obs);

                split.Add(new()
                {
                    ["galaxy"] = name,
                    ["index"] = i,
                    ["partition"] = role,
                    ["population_role"] = population
                });
            }
        }

        Csv.Write(Path.Combine(outputRoot, "predictors.csv"), predictors);
        Csv.Write(Path.Combine(outputRoot, "training.csv"), train);
        Csv.Write(Path.Combine(outputRoot, "evaluation.csv"), evaluation);

        Save(Path.Combine(outputRoot, "metadata.json"), metadata);
        Save(Path.Combine(outputRoot, "split.json"), split);

        var protocol = new
        {
            schema_version = 1,
            study_status = "retrospective_development",
            dataset_source = "Retained SPARC catalogue; previously examined in the consolidated manuscript",
            acquisition_statement = "Existing explored data; no independently acquired or blinded observations",
            previously_inspected = true,
            files = new Dictionary<string, string>{
{
"predictors","predictors.csv"}
,{
"training","training.csv"}
,{
"evaluation","evaluation.csv"}
,{
"metadata","metadata.json"}
,{
"split","split.json"}
}
,
            models = new[]{
"core","nfw","profiled_baryons","envelope"}
,
            primary_model = "envelope",
            evaluation_galaxies = new[]{
"DDO154","NGC0300"}
,
            shared_parameter_models = Array.Empty<string>(),
            shared_training_galaxies = Array.Empty<string>(),
            fitting_rules = new
            {
                stellar_centres = new[]{
.5,.7}
,
                stellar_width_dex = .1,
                pull_bound = 6,
                distance_ratio_floor = .001,
                core_depth_guard = .001,
                envelope_taper = 20,
                training_weight = "pooled_points",
                failure_policy = "retain_failed_fit_no_finite_penalty",
                calibration = "same_prior_centres_and_bounds; envelope conditional on newly fitted guarded core",
                optimizer_budget = "original_multistart_and_tolerances",
                shared_parameter_policy = "learn_only_population_train_galaxies_then_freeze_for_transfer",
                local_state_policy = "inner_rows_only_in_both_train_and_transfer_galaxies"
            }
,
            pilot = new
            {
                backend = "empirical_original",
                perturbation_replicates = 2,
                seed = 2026092801L,
                timeout_seconds = 3600,
                measurement_model = "independent Gaussian training velocity perturbations about observed velocities; Gaussian D/i prior-centre perturbations within physical domain; fixed stellar priors; sensitivity ensemble, not posterior or calibrated confidence interval"
            }
,
            inference = new
            {
                enabled = true,
                minimum_galaxies = 20,
                bootstrap_draws = 999,
                confidence_level = .95,
                seed = 2026092802UL
            }
,
            ensemble_contract = new
            {
                kind = "nuisance_sensitivity",
                generative_model = "The declared pilot measurement model; no population resampling or shared physical recalibration",
                required_refit_stages = new[]{
"core","nfw","profiled_baryons","conditional_envelope"}
            }
,
            evaluation_policy = "equal galaxy mean outer squared standardized loss and MAE; failures retained; no pilot p-values",
            supported_backend = "empirical_original",
            external_backends = new[]{
"pressure_shared_scale","charged_shared_mass"}
,
            catalogue_sha256 = data.CatalogueSha256,
            mass_models_sha256 = data.MassModelsSha256
        };

        string path = Path.Combine(outputRoot, "protocol.json");
        Save(path, protocol);
        return new
        {
            status = "draft",
            protocol = path,
            galaxies = galaxies.Length,
            training_rows = train.Count,
            evaluation_rows = evaluation.Count,
            next = "Review/edit protocol and data, then research-lock. Existing data remain retrospective."
        };
    }

    static (JsonObject Protocol, List<Dictionary<string, string>> Predictors) ValidateDraft(string protocolPath)
    {
        var p = Read(protocolPath);
        string root = Path.GetDirectoryName(Path.GetFullPath(protocolPath))!;

        Check(p["schema_version"]!.GetValue<int>() == 1, "Unsupported research protocol version.");

        string status = p["study_status"]!.GetValue<string>();
        Check(status is "retrospective_development" or "prospective_candidate", "Declare retrospective_development or prospective_candidate.");

        Check(!string.IsNullOrWhiteSpace(p["dataset_source"]!.GetValue<string>()) && !string.IsNullOrWhiteSpace(p["acquisition_statement"]!.GetValue<string>()), "Dataset source and acquisition statement required.");

        Check(status != "prospective_candidate" || !p["previously_inspected"]!.GetValue<bool>(), "Previously inspected data cannot be declared prospective.");

        var files = p["files"]!.AsObject();
        Check(files.Select(f => f.Key).ToHashSet().SetEquals(new[]{
"predictors","training","evaluation","metadata","split"}) && files.Select(f => f.Value!.GetValue<string>()).Distinct().Count() == 5, "Exactly five separate input file roles required.");
        Local(root, Path.GetFileName(protocolPath));
        string FilePath(string n) => Local(root, files[n]!.GetValue<string>());

        var predictors = Csv.Read(FilePath("predictors"));
        var tr = Csv.Read(FilePath("training"));
        var ev = Csv.Read(FilePath("evaluation"));

        var keys = predictors.Select(Key).ToHashSet();
        Check(keys.Count == predictors.Count && keys.Count > 0, "Duplicate/empty predictor keys.");

        Check(predictors.All(r => !r.Keys.Any(k => k.Contains("observed", StringComparison.OrdinalIgnoreCase))), "Predictor table must not contain outcomes.");

        var tk = tr.Select(Key).ToHashSet();
        var ek = ev.Select(Key).ToHashSet();
        Check(tk.Count == tr.Count && ek.Count == ev.Count && !tk.Overlaps(ek), "Duplicate or overlapping training/evaluation keys.");

        Check(keys.SetEquals(tk.Concat(ek)), "Training/evaluation must partition every predictor key.");

        foreach (var r in predictors)
        {
            Check(r["partition"] is "inner" or "outer", "Unknown radial partition.");
            Check(r["population_role"] is "train" or "transfer", "Unknown population role.");
            Check((r["partition"] == "inner") == tk.Contains(Key(r)), "Split and outcome partitions disagree.");
            foreach (string field in new[]{
"radius_kpc","sigma_kms","gas_kms","disk_kms","bulge_kms"}) Check(double.IsFinite(Csv.Number(r, field)), "Nonfinite predictor: " + field);
            Check(Csv.Number(r, "radius_kpc") > 0 && Csv.Number(r, "sigma_kms") > 0, "Invalid radius or error.");
        }

        var observations = tr.Concat(ev).ToDictionary(Key);
        foreach (var r in predictors)
        {
            var o = observations[Key(r)];
            Check(double.IsFinite(Csv.Number(o, "observed_kms")) && Csv.Number(o, "sigma_kms") == Csv.Number(r, "sigma_kms"), "Invalid outcome or differing observational error.");
        }

        var metadata = Read(FilePath("metadata"));
        foreach (var group in predictors.GroupBy(r => r["galaxy"]))
        {
            var rows = group.OrderBy(r => int.Parse(r["index"], CultureInfo.InvariantCulture)).ToArray();
            Check(rows.Select(r => r["population_role"]).Distinct().Count() == 1, "Galaxy straddles population training and transfer.");
            Check(rows.Select((r, i) => int.Parse(r["index"], CultureInfo.InvariantCulture) == i).All(x => x), "Indices must be contiguous from zero.");
            Check(rows.Count(r => r["partition"] == "inner") >= 3 && rows.Any(r => r["partition"] == "outer"), "Each galaxy needs inner calibration and outer evaluation.");
            Check(rows.Select(r => r["partition"]).SequenceEqual(rows.OrderBy(r => r["partition"] == "inner" ? 0 : 1).Select(r => r["partition"])), "Inner rows must precede outer rows.");
            Check(rows.Select(r => Csv.Number(r, "radius_kpc")).SequenceEqual(rows.Select(r => Csv.Number(r, "radius_kpc")).Order()), "Radii must be nondecreasing with index.");
            Check(metadata.ContainsKey(group.Key), "Missing galaxy metadata.");
            var m = metadata[group.Key]!;
            foreach (string field in new[]{
"D","e_D","inc","e_inc","Rdisk","catalogue_D","catalogue_inc"}) Check(double.IsFinite(m[field]!.GetValue<double>()), "Nonfinite metadata: " + field);
            Check(m["D"]!.GetValue<double>() > 0 && m["catalogue_D"]!.GetValue<double>() > 0 && m["e_D"]!.GetValue<double>() >= 0 && m["Rdisk"]!.GetValue<double>() > 0, "Invalid distance/error/disk scale.");
            Check(m["inc"]!.GetValue<double>() > 0 && m["inc"]!.GetValue<double>() <= 90 && m["catalogue_inc"]!.GetValue<double>() > 0 && m["catalogue_inc"]!.GetValue<double>() <= 90 && m["e_inc"]!.GetValue<double>() >= 0, "Invalid inclination/error.");
        }

        var split = JsonNode.Parse(System.IO.File.ReadAllText(FilePath("split")))!.AsArray();
        Check(split.Count == predictors.Count, "Split manifest coverage mismatch.");
        var byKey = predictors.ToDictionary(Key);
        var sk = new HashSet<string>();
        foreach (var r in split)
        {
            string k = r!["galaxy"]!.GetValue<string>() + ":" + r["index"]!.ToString();
            Check(sk.Add(k) && byKey.ContainsKey(k), "Invalid split key.");
            Check(r["partition"]!.GetValue<string>() == byKey[k]["partition"] && r["population_role"]!.GetValue<string>() == byKey[k]["population_role"], "Split role mismatch.");
        }

        var models = Strings(p["models"]);
        Check(models.Length > 1 && models.Distinct().Count() == models.Length && models.Contains(p["primary_model"]!.GetValue<string>()), "Invalid model family.");
        var chosen = Strings(p["evaluation_galaxies"]);
        Check(chosen.Length > 0 && chosen.Distinct().Count() == chosen.Length && chosen.All(g => predictors.Any(r => r["galaxy"] == g)), "Invalid declared evaluation galaxies.");

        Check(p["fitting_rules"]!["training_weight"]!.GetValue<string>() is "pooled_points" or "equal_galaxy", "Unsupported training weighting.");

        var sharedModels = Strings(p["shared_parameter_models"]);
        var sharedGalaxies = Strings(p["shared_training_galaxies"]);
        Check(sharedModels.Distinct().Count() == sharedModels.Length && sharedModels.All(models.Contains), "Invalid shared model declarations.");
        Check(sharedGalaxies.Distinct().Count() == sharedGalaxies.Length && sharedGalaxies.All(g => predictors.Any(r => r["galaxy"] == g && r["population_role"] == "train")), "Shared training galaxies must belong to population training.");
        Check((sharedModels.Length == 0) == (sharedGalaxies.Length == 0), "Shared models and training galaxies must both be declared.");

        Check(p["ensemble_contract"]!["kind"]!.GetValue<string>() is "baseline_only" or "nuisance_sensitivity" or "full_fitted_procedure", "Unknown ensemble contract.");
        Check(!string.IsNullOrWhiteSpace(p["ensemble_contract"]!["generative_model"]!.GetValue<string>()), "Explicit generative/sensitivity model required.");

        var inference = p["inference"]!;
        Check(inference["minimum_galaxies"]!.GetValue<int>() >= 2 && inference["bootstrap_draws"]!.GetValue<int>() >= 99 && inference["confidence_level"]!.GetValue<double>() > 0 && inference["confidence_level"]!.GetValue<double>() < 1, "Invalid inference protocol.");

        return (p, predictors);
    }

    public static object LockProtocol(string protocolPath)
    {
        var (p, rows) = ValidateDraft(protocolPath);
        string root = Path.GetDirectoryName(Path.GetFullPath(protocolPath))!;
        string target = Path.Combine(root, "protocol.lock.json");
        Check(!File.Exists(target), "Protocol already locked; create a new study for changes.");

        var files = p["files"]!.AsObject().ToDictionary(k => k.Value!.GetValue<string>(), v => Hash(Local(root, v.Value!.GetValue<string>())));

        Save(target, new
        {
            schema_version = 1,
            protocol_sha256 = Hash(protocolPath),
            files,
            locked_at_utc = DateTimeOffset.UtcNow,
            fitting_rules_json = p["fitting_rules"]!.ToJsonString(),
            independence_verified = false
        });

        return new
        {
            status = "locked",
            path = target,
            rows = rows.Count,
            independence_verified = false
        };
    }

    public static object ValidateProtocol(string protocolPath)
    {
        var (p, rows) = ValidateDraft(protocolPath);
        string root = Path.GetDirectoryName(Path.GetFullPath(protocolPath))!;
        var l = Read(Local(root, "protocol.lock.json"));
        Check(l["files"]!.AsObject().Select(f => f.Key).ToHashSet().SetEquals(p["files"]!.AsObject().Select(f => f.Value!.GetValue<string>())), "Lock must cover every declared file exactly.");
        Check(l["fitting_rules_json"]!.GetValue<string>() == p["fitting_rules"]!.ToJsonString(), "Lock fitting rules changed.");
        Check(l["protocol_sha256"]!.GetValue<string>() == Hash(protocolPath), "Locked protocol changed.");
        foreach (var f in l["files"]!.AsObject()) Check(Hash(Local(root, f.Key)) == f.Value!.GetValue<string>(), "Locked input changed: " + f.Key);

        return new
        {
            status = "pass",
            protocol_sha256 = Hash(protocolPath),
            rows = rows.Count,
            study_status = p["study_status"]!.GetValue<string>(),
            independence_verified = false
        };
    }

    public static object Evaluate(string protocolPath, string predictionsPath, string fitArtifactPath, string outputRoot)
    {
        ValidateProtocol(protocolPath);
        var (p, inputs) = ValidateDraft(protocolPath);
        ExecutionSafety.ValidateRoots(Path.GetDirectoryName(Path.GetFullPath(protocolPath))!, outputRoot);
        Check(!File.Exists(Path.Combine(outputRoot, "research_evaluation.json")), "Use fresh evaluation output.");
        string root = Path.GetDirectoryName(Path.GetFullPath(protocolPath))!;
        var a = Read(fitArtifactPath);
        var chosen = Strings(p["evaluation_galaxies"]).ToHashSet();
        var models = Strings(p["models"]);
        Check(a["protocol_sha256"]!.GetValue<string>() == Hash(protocolPath) && a["predictions_sha256"]!.GetValue<string>() == Hash(predictionsPath), "Fit artifact protocol/prediction identity mismatch.");
        Check(a["fitting_rules_json"]!.GetValue<string>() == p["fitting_rules"]!.ToJsonString(), "Fit used different priors, bounds, budgets or weighting.");
        Check(a["read_evaluation_outcomes"]!.GetValue<bool>() == false, "Fitter reports reading evaluation outcomes.");

        var train = inputs.Where(r => r["partition"] == "inner" && chosen.Contains(r["galaxy"])).Select(Key).ToHashSet();
        Check(train.SetEquals(Strings(a["training_keys"])) && Strings(a["training_keys"]).Distinct().Count() == Strings(a["training_keys"]).Length, "Fitter training keys differ from declared inner data.");

        ValidateEnsembleContract(p, inputs, a, fitArtifactPath);

        int reps = a["replicates"]!.GetValue<int>();
        Check(reps > 0, "No fit replicates.");
        var predictions = Csv.Read(predictionsPath);
        var expected = inputs.Where(r => r["partition"] == "outer" && chosen.Contains(r["galaxy"])).Select(Key).ToHashSet();
        var seen = new HashSet<string>();
        foreach (var r in predictions)
        {
            int rep = CanonicalReplicate(r["replicate"]);
            Check(rep >= 0 && rep < reps && models.Contains(r["model"]) && expected.Contains(Key(r)), "Unexpected prediction key.");
            Check(seen.Add(rep + ":" + r["model"] + ":" + Key(r)), "Duplicate prediction.");
            Check(r["fit_failed"].ToLowerInvariant() is "true" or "false", "Explicit fit failure flag required.");
            if (!Csv.Flag(r, "fit_failed")) Check(double.IsFinite(Csv.Number(r, "predicted_kms")) && Csv.Number(r, "predicted_kms") >= 0, "Invalid successful prediction.");
        }

        Check(seen.Count == reps * models.Length * expected.Count, "Prediction coverage incomplete.");
        var obs = Csv.Read(Local(root, p["files"]!["evaluation"]!.GetValue<string>())).ToDictionary(Key);
        var scores = new List<Dictionary<string, object?>>();
        foreach (var group in predictions.GroupBy(r => (r["replicate"], r["model"], r["galaxy"])))
        {
            bool failed = group.Any(r => Csv.Flag(r, "fit_failed"));
            double[] e = failed ? [] : group.Select(r => Csv.Number(r, "predicted_kms") - Csv.Number(obs[Key(r)], "observed_kms")).ToArray();
            scores.Add(new()
            {
                ["replicate"] = group.Key.Item1,
                ["model"] = group.Key.Item2,
                ["galaxy"] = group.Key.Item3,
                ["fit_failed"] = failed,
                ["outer_points"] = group.Count(),
                ["squared_loss"] = failed ? null : group.Average(r => Math.Pow((Csv.Number(r, "predicted_kms") - Csv.Number(obs[Key(r)], "observed_kms")) / Csv.Number(obs[Key(r)], "sigma_kms"), 2)),
                ["mae_kms"] = failed ? null : e.Average(Math.Abs)
            });
        }

        foreach (var score in scores)
            if (!(bool)score["fit_failed"]!)
                Check(double.IsFinite((double)score["squared_loss"]!) && double.IsFinite((double)score["mae_kms"]!), "Nonfinite derived score: prediction/outcome/error scale overflows for " + score["model"] + "/" + score["galaxy"]);
        Directory.CreateDirectory(outputRoot);
        Csv.Write(Path.Combine(outputRoot, "research_galaxy_scores.csv"), scores);
        var summary = scores.GroupBy(r => (r["replicate"], r["model"])).Select(g => new
        {
            replicate = g.Key.Item1,
            model = g.Key.Item2,
            failed_galaxies = g.Count(r => (bool)r["fit_failed"]!),
            full_mean_loss = g.Any(r => (bool)r["fit_failed"]!) ? (double?)null : g.Average(r => (double)r["squared_loss"]!),
            finite_subset_mean_loss = g.Any(r => !(bool)r["fit_failed"]!) ? g.Where(r => !(bool)r["fit_failed"]!).Average(r => (double)r["squared_loss"]!) : (double?)null
        }).ToArray();

        foreach (var row in summary)
            Check((row.full_mean_loss is not double full || double.IsFinite(full)) && (row.finite_subset_mean_loss is not double subset || double.IsFinite(subset)), "Nonfinite aggregate score; numerical scale is unsupported.");
        var report = new
        {
            status = "pass",
            study_status = p["study_status"]!.GetValue<string>(),
            independence_verified = false,
            scope = "Declared development/transfer predictions; nuisance perturbation pilot is not posterior inference, a calibrated confidence interval or new independent evidence.",
            replicates = reps,
            galaxies = chosen.Count,
            models,
            summary,
            source_artifact_sha256 = Hash(fitArtifactPath),
            comparisons = CompareBaseline(p, scores, chosen.Count),
            contract_scope = "External provenance is self-attested consistency, not proof of scientific correctness or outcome blindness. Conditional baseline inference assumes independent galaxy clusters. Minimum galaxy count is an operational gate, not calibration certification. No calibrated full-procedure interval is produced."
        };
        Save(Path.Combine(outputRoot, "research_evaluation.json"), report);
        return report;
    }

    static int CanonicalReplicate(string text)
    {
        Check(int.TryParse(text, NumberStyles.None, CultureInfo.InvariantCulture, out int value) && value >= 0 && text == value.ToString(CultureInfo.InvariantCulture), "Replicate IDs must be canonical nonnegative integers.");
        return value;
    }

    // External evidence is checked for consistency, not treated as proof of a fitter's equations.
    static void ValidateEnsembleContract(JsonObject p, List<Dictionary<string, string>> inputs, JsonObject a, string artifactPath)
    {
        string root = Path.GetDirectoryName(Path.GetFullPath(artifactPath))!;

        string kind = p["ensemble_contract"]!["kind"]!.GetValue<string>();

        Check(a["ensemble_kind"]!.GetValue<string>() == kind, "Different ensemble design.");

        int reps = a["replicates"]!.GetValue<int>();

        Check(reps > 0 && (kind != "baseline_only" || reps == 1), "Invalid replicate count.");

        var models = Strings(p["models"]);
        var sharedModels = Strings(p["shared_parameter_models"]);
        var galaxies = Strings(p["shared_training_galaxies"]).ToHashSet();

        var expected = inputs.Where(r => r["partition"] == "inner" && galaxies.Contains(r["galaxy"])).Select(Key).ToHashSet();

        var keys = Strings(a["shared_training_keys"]);

        Check(a["has_shared_parameters"]!.GetValue<bool>() == (sharedModels.Length > 0), "Shared parameter flag contradicts declared models.");

        Check(keys.Distinct().Count() == keys.Length && expected.SetEquals(keys), "Shared training keys differ from declared population inner data.");

        if (sharedModels.Length > 0)
        {
            var seen = new HashSet<string>();

            foreach (var row in a["shared_parameter_fits"]!.AsArray())
            {
                int rep = row!["replicate"]!.GetValue<int>();
                string model = row["model"]!.GetValue<string>();

                Check(rep >= 0 && rep < reps && sharedModels.Contains(model) && seen.Add(rep + ":" + model), "Invalid shared fit record.");

                Check(expected.SetEquals(Strings(row["training_keys"])), "Shared fit used undeclared outcomes.");

                string file = Local(root, row["parameters_file"]!.GetValue<string>());

                Check(Hash(file) == row["parameters_sha256"]!.GetValue<string>() && Hash(file) == row["transfer_parameters_sha256"]!.GetValue<string>(), "Shared parameters changed during transfer.");
            }

            Check(seen.Count == reps * sharedModels.Length, "Each replicate/shared model needs a calibration and frozen transfer artifact.");
        }

        if (kind == "full_fitted_procedure")
        {
            Check(reps >= 2 && galaxies.Count > 0, "Full-fit ensemble needs population training and refits.");

            var records = a["replicate_provenance"]!.AsArray();
            var seen = new HashSet<int>();

            Check(records.Count == reps, "Every full-fit replicate needs provenance.");

            var stages = Strings(p["ensemble_contract"]!["required_refit_stages"]).ToHashSet();
            Check(stages.Count > 0, "Declare every affected refit stage.");

            foreach (var row in records)
            {
                int rep = row!["replicate"]!.GetValue<int>();
                Check(rep >= 0 && rep < reps && seen.Add(rep), "Invalid replicate provenance.");

                var counts = row["population_multiplicities"]!.AsObject();

                Check(counts.Select(x => x.Key).ToHashSet().SetEquals(galaxies) && counts.All(x => x.Value!.GetValue<int>() >= 0) && counts.Sum(x => x.Value!.GetValue<int>()) == galaxies.Count, "Invalid population bootstrap multiplicities.");

                Check(rep != 0 || counts.All(x => x.Value!.GetValue<int>() == 1), "Baseline must use every training galaxy once.");

                Check(stages.SetEquals(Strings(row["refitted_stages"])) && models.ToHashSet().SetEquals(Strings(row["paired_models"])), "All affected stages/models must share each declared replicate.");

                foreach (string field in new[]{
"measurement_draws","fit_records"}) Check(Hash(Local(root, row[field + "_file"]!.GetValue<string>())) == row[field + "_sha256"]!.GetValue<string>(), "Full-fit evidence changed.");
            }
        }
    }

    static object[] CompareBaseline(JsonObject p, List<Dictionary<string, object?>> scores, int galaxyCount)
    {
        string primary = p["primary_model"]!.GetValue<string>();
        var comparators = Strings(p["models"]).Where(m => m != primary).ToArray();
        var rules = p["inference"]!;

        var baseline = scores.Where(r => (string)r["replicate"]! == "0").ToDictionary(r => ((string)r["model"]!, (string)r["galaxy"]!));

        bool sharedTrainingReused = Strings(p["shared_parameter_models"]).Length > 0 && Strings(p["evaluation_galaxies"]).Intersect(Strings(p["shared_training_galaxies"])).Any();
        var result = new List<Dictionary<string, object?>>();

        for (int m = 0; m < comparators.Length; m++)
        {
            string comparator = comparators[m];

            var paired = baseline.Values.Where(r => (string)r["model"]! == primary && !(bool)r["fit_failed"]! && !(bool)baseline[(comparator, (string)r["galaxy"]!)]["fit_failed"]!).Select(r => new PairedGalaxyObservation((string)r["galaxy"]!, (int)r["outer_points"]!, (double)baseline[(comparator, (string)r["galaxy"]!)]["squared_loss"]! - (double)r["squared_loss"]!)).ToArray();

            bool complete = paired.Length == galaxyCount;
            bool allowed = complete && !sharedTrainingReused && rules["enabled"]!.GetValue<bool>() && galaxyCount >= rules["minimum_galaxies"]!.GetValue<int>();

            PairedBootstrapResult? inference = allowed ? StatisticalEngine.Run(paired, new()
            {
                Draws = rules["bootstrap_draws"]!.GetValue<int>(),
                ConfidenceLevel = rules["confidence_level"]!.GetValue<double>(),
                PairsSeed = StatisticalEngine.DeriveSeed(rules["seed"]!.GetValue<ulong>(), (ulong)(2 * m)),
                WildSeed = StatisticalEngine.DeriveSeed(rules["seed"]!.GetValue<ulong>(), (ulong)(2 * m + 1))
            }) : null;

            Check(paired.Length == 0 || double.IsFinite(paired.Average(r => r.Difference)), "Nonfinite paired mean; numerical scale is unsupported.");
            result.Add(new()
            {
                ["primary"] = primary,
                ["comparator"] = comparator,
                ["difference_direction"] = "comparator minus primary; positive favors primary",
                ["status"] = !complete ? "incomplete_family_member_no_full_sample_inference" : sharedTrainingReused ? "descriptive_shared_training_reuse" : !allowed ? "descriptive_operational_gate" : inference!.Status,
                ["jointly_feasible_galaxies"] = paired.Length,
                ["full_sample_effect"] = complete ? paired.Average(r => r.Difference) : (double?)null,
                ["joint_feasible_subset_effect"] = paired.Length > 0 ? paired.Average(r => r.Difference) : (double?)null,
                ["conditional_baseline_inference"] = inference,
                ["wild_p_value"] = inference?.WildPValue,
                ["holm_p_value"] = null
            });
        }

        // Failed/undefined slots contribute one only internally, preserving the declared family.
        var ordered = result.Select((r, i) => (i, p: r["wild_p_value"] is double v ? v : 1d)).OrderBy(x => x.p).ToArray();
        double previous = 0;

        for (int rank = 0; rank < ordered.Length; rank++)
        {
            previous = Math.Max(previous, Math.Min(1, (ordered.Length - rank) * ordered[rank].p));
            if (result[ordered[rank].i]["wild_p_value"] is double) result[ordered[rank].i]["holm_p_value"] = previous;
        }

        return result.Cast<object>().ToArray();
    }

    /// <summary>Adversarial protocol/fit-contract checks on temporary synthetic data.</summary>
    public static object SelfCheck()
    {
        string root = Path.Combine(Path.GetTempPath(), "du-research-selfcheck-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(root); var passed = new List<string>();
        string Setup(string name)
        {
            string folder = Path.Combine(root, name); Directory.CreateDirectory(folder);
            File.WriteAllText(Path.Combine(folder, "predictors.csv"), "galaxy,index,radius_kpc,sigma_kms,gas_kms,disk_kms,bulge_kms,partition,population_role\n" + string.Join("\n", new[] { "A", "B" }.SelectMany(g => Enumerable.Range(0, 4).Select(i => $"{g},{i},{i + 1},2,3,4,0,{(i < 3 ? "inner" : "outer")},{(g == "A" ? "train" : "transfer")}"))) + "\n");
            foreach (string role in new[] { "training", "evaluation" }) File.WriteAllText(Path.Combine(folder, role + ".csv"), "galaxy,index,observed_kms,sigma_kms\n" + string.Join("\n", new[] { "A", "B" }.SelectMany(g => (role == "training" ? new[] { 0, 1, 2 } : new[] { 3 }).Select(i => $"{g},{i},10,2"))) + "\n");
            Save(Path.Combine(folder, "metadata.json"), new[] { "A", "B" }.ToDictionary(g => g, g => (object)new { D = 5, e_D = .1, inc = 45, e_inc = 2, Rdisk = 1, catalogue_D = 5, catalogue_inc = 45 }));
            Save(Path.Combine(folder, "split.json"), new[] { "A", "B" }.SelectMany(g => Enumerable.Range(0, 4).Select(i => new { galaxy = g, index = i, partition = i < 3 ? "inner" : "outer", population_role = g == "A" ? "train" : "transfer" })).ToArray());
            Save(Path.Combine(folder, "protocol.json"), new { schema_version = 1, study_status = "retrospective_development", dataset_source = "synthetic test", acquisition_statement = "synthetic, not observational evidence", previously_inspected = true, files = new Dictionary<string, string> { { "predictors", "predictors.csv" }, { "training", "training.csv" }, { "evaluation", "evaluation.csv" }, { "metadata", "metadata.json" }, { "split", "split.json" } }, models = new[] { "core", "nfw" }, primary_model = "core", evaluation_galaxies = new[] { "A", "B" }, shared_parameter_models = Array.Empty<string>(), shared_training_galaxies = Array.Empty<string>(), fitting_rules = new { training_weight = "pooled_points" }, ensemble_contract = new { kind = "baseline_only", generative_model = "synthetic baseline", required_refit_stages = new[] { "core", "nfw" } }, inference = new { enabled = true, minimum_galaxies = 20, bootstrap_draws = 99, confidence_level = .95, seed = 1UL } });
            return Path.Combine(folder, "protocol.json");
        }
        void Test(string name, Action<string> action, bool rejects)
        {
            string protocol = Setup(name); bool rejected = false;
            try { action(protocol); } catch (InvalidDataException) { rejected = true; }
            Check(rejected == rejects, "Research self-check failed: " + name); passed.Add(name);
        }
        void Edit(string path, Action<JsonObject> change) { var node = Read(path); change(node); Save(path, node); }
        Test("valid-lock", p => { LockProtocol(p); ValidateProtocol(p); }, false);
        foreach (string field in new[] { "radius_kpc", "sigma_kms", "gas_kms" })
            Test("nonfinite-" + field, p => { string f = Path.Combine(Path.GetDirectoryName(p)!, "predictors.csv"); var rows = Csv.Read(f); rows[0][field] = "Infinity"; Csv.Write(f, rows.Select(r => r.ToDictionary(k => k.Key, k => (object?)k.Value))); LockProtocol(p); }, true);
        Test("overlap", p => { string f = Path.Combine(Path.GetDirectoryName(p)!, "evaluation.csv"); File.WriteAllText(f, File.ReadAllText(f).Replace("A,3,", "A,0,")); LockProtocol(p); }, true);
        Test("invalid-inclination", p => { Edit(Path.Combine(Path.GetDirectoryName(p)!, "metadata.json"), x => x["A"]!["inc"] = 0); LockProtocol(p); }, true);
        Test("radius-order", p => { string f = Path.Combine(Path.GetDirectoryName(p)!, "predictors.csv"); File.WriteAllText(f, File.ReadAllText(f).Replace("A,0,1,", "A,0,8,")); LockProtocol(p); }, true);
        Test("omitted-lock-entry", p => { LockProtocol(p); Edit(Path.Combine(Path.GetDirectoryName(p)!, "protocol.lock.json"), x => x["files"]!.AsObject().Remove("evaluation.csv")); ValidateProtocol(p); }, true);
        Test("changed-input", p => { LockProtocol(p); string f = Path.Combine(Path.GetDirectoryName(p)!, "evaluation.csv"); File.WriteAllText(f, File.ReadAllText(f).Replace(",10,", ",11,")); ValidateProtocol(p); }, true);
        Test("false-prospective", p => { Edit(p, x => x["study_status"] = "prospective_candidate"); LockProtocol(p); }, true);
        foreach (string scenario in new[] { "shared-flag", "shared-outer-leak", "shared-transfer-change", "valid-shared-transfer" })
            Test(scenario, p =>
            {
                Edit(p, x => { x["shared_parameter_models"] = JsonSerializer.SerializeToNode(new[] { "core" }); x["shared_training_galaxies"] = JsonSerializer.SerializeToNode(new[] { "A" }); }); LockProtocol(p);
                var (protocol, inputs) = ValidateDraft(p); string folder = Path.GetDirectoryName(p)!; string state = Path.Combine(folder, "parameters.json"); File.WriteAllText(state, "{}");
                var keys = new[] { "A:0", "A:1", "A:2" }; var artifact = JsonSerializer.SerializeToNode(new { ensemble_kind = "baseline_only", replicates = 1, has_shared_parameters = scenario != "shared-flag", shared_training_keys = scenario == "shared-outer-leak" ? new[] { "A:0", "A:1", "A:3" } : keys, shared_parameter_fits = new[] { new { replicate = 0, model = "core", training_keys = keys, parameters_file = "parameters.json", parameters_sha256 = Hash(state), transfer_parameters_sha256 = scenario == "shared-transfer-change" ? "bad" : Hash(state) } } })!.AsObject();
                ValidateEnsembleContract(protocol, inputs, artifact, Path.Combine(folder, "artifact.json"));
            }, scenario != "valid-shared-transfer");
        Test("descriptive-family", p =>
        {
            var protocol = Read(p); var scores = new List<Dictionary<string, object?>>(); foreach (string model in new[] { "core", "nfw" }) foreach (string galaxy in new[] { "A", "B" }) scores.Add(new() { ["replicate"] = "0", ["model"] = model, ["galaxy"] = galaxy, ["outer_points"] = 1, ["fit_failed"] = false, ["squared_loss"] = model == "core" ? 1d : 2d });
            var result = JsonSerializer.SerializeToNode(CompareBaseline(protocol, scores, 2))!.AsArray(); Check(result.Count == 1 && result[0]!["wild_p_value"] == null && result[0]!["full_sample_effect"]!.GetValue<double>() == 1, "Pilot gate/effect incorrect.");
            scores[3]["fit_failed"] = true; result = JsonSerializer.SerializeToNode(CompareBaseline(protocol, scores, 2))!.AsArray(); Check(result.Count == 1 && result[0]!["full_sample_effect"] == null && result[0]!["wild_p_value"] == null, "Failure was dropped from family.");
        }, false);
        foreach (string invalid in new[] { "00", "+0", " 0", "-1" })
            Test("replicate-" + invalid.Replace(' ', '_'), _ => CanonicalReplicate(invalid), true);
        Test("derived-score-overflow", path =>
        {
            LockProtocol(path); string folder = Path.GetDirectoryName(path)!; var protocol = Read(path);
            string predictions = Path.Combine(folder, "predictions.csv");
            File.WriteAllText(predictions, "replicate,galaxy,index,model,predicted_kms,fit_failed\n0,A,3,core,1e308,False\n0,B,3,core,1e308,False\n0,A,3,nfw,10,False\n0,B,3,nfw,10,False\n");
            string artifact = Path.Combine(folder, "artifact.json");
            Save(artifact, new { protocol_sha256 = Hash(path), predictions_sha256 = Hash(predictions), fitting_rules_json = protocol["fitting_rules"]!.ToJsonString(), read_evaluation_outcomes = false, training_keys = new[] { "A:0", "A:1", "A:2", "B:0", "B:1", "B:2" }, shared_training_keys = Array.Empty<string>(), has_shared_parameters = false, replicates = 1, ensemble_kind = "baseline_only" });
            Evaluate(path, predictions, artifact, Path.Combine(root, "overflow-output"));
        }, true);
        Test("shared-training-inference-gate", path =>
        {
            var protocol = Read(path); protocol["shared_parameter_models"] = JsonSerializer.SerializeToNode(new[] { "core" }); protocol["shared_training_galaxies"] = JsonSerializer.SerializeToNode(new[] { "A" }); protocol["inference"]!["minimum_galaxies"] = 2;
            var scores = new List<Dictionary<string, object?>>(); foreach (string model in new[] { "core", "nfw" }) foreach (string galaxy in new[] { "A", "B" }) scores.Add(new() { ["replicate"] = "0", ["model"] = model, ["galaxy"] = galaxy, ["outer_points"] = 1, ["fit_failed"] = false, ["squared_loss"] = model == "core" ? 1d : 2d });
            var result = JsonSerializer.SerializeToNode(CompareBaseline(protocol, scores, 2))!.AsArray(); Check(result[0]!["status"]!.GetValue<string>() == "descriptive_shared_training_reuse" && result[0]!["wild_p_value"] == null, "Shared training reuse must disable inference.");
        }, false);
        return new { status = "pass", checks = passed.Count, passed, temporary_fixture_root = root };
    }

    public static object RunPilot(string protocolPath, string referenceRoot, string outputRoot, string pythonExecutable, string? packagesPath = null)
    {
        ValidateProtocol(protocolPath);
        ExecutionSafety.ValidateRoots(referenceRoot, outputRoot);
        ExecutionSafety.ValidateRoots(Path.GetDirectoryName(Path.GetFullPath(protocolPath))!, outputRoot);
        Check(!Inside(referenceRoot, outputRoot) && !Inside(Path.GetDirectoryName(Path.GetFullPath(protocolPath))!, outputRoot), "Pilot output must be outside frozen reference and locked study.");

        string runner = Path.Combine(Path.GetDirectoryName(Path.GetFullPath(referenceRoot))!, "DarkMatter", "tools", "research_python_runner.py");
        if (!File.Exists(runner)) runner = Path.Combine(AppContext.BaseDirectory, "tools", "research_python_runner.py");
        Check(File.Exists(runner), "Research Python adapter not found.");
        Directory.CreateDirectory(outputRoot);

        var start = new ProcessStartInfo(pythonExecutable)
        {
            UseShellExecute = false,
            WorkingDirectory = outputRoot,
            RedirectStandardOutput = true,
            RedirectStandardError = true
        };
        foreach (string arg in new[]{
"-B",runner,"--protocol",Path.GetFullPath(protocolPath),"--reference",Path.GetFullPath(referenceRoot),"--out",Path.GetFullPath(outputRoot)}) start.ArgumentList.Add(arg);
        if (packagesPath != null)
        {
            Check(Directory.Exists(packagesPath), "Python package path missing.");
            start.Environment["PYTHONPATH"] = Path.GetFullPath(packagesPath);
        }
        start.Environment["PYTHONDONTWRITEBYTECODE"] = "1";
        start.Environment["OPENBLAS_NUM_THREADS"] = "1";
        start.Environment["OMP_NUM_THREADS"] = "1";

        int timeout = Read(protocolPath)["pilot"]!["timeout_seconds"]!.GetValue<int>();
        Check(timeout >= 1 && timeout <= 86400, "Pilot timeout must be 1..86400 seconds.");

        using var process = Process.Start(start) ?? throw new InvalidOperationException("Python did not start.");

        var stdout = process.StandardOutput.ReadToEndAsync();
        var stderr = process.StandardError.ReadToEndAsync();
        bool cancelled = false;

        ConsoleCancelEventHandler cancel = (_, e) =>
        {
            e.Cancel = true;
            cancelled = true;
            if (!process.HasExited) process.Kill(true);
        };
        Console.CancelKeyPress += cancel;

        bool completed;
        try
        {
            completed = process.WaitForExit(timeout * 1000);
            if (!completed)
            {
                process.Kill(true);
                process.WaitForExit();
            }
        }
        finally
        {
            Console.CancelKeyPress -= cancel;
        }

        File.WriteAllText(Path.Combine(outputRoot, "pilot.stdout.log"), stdout.GetAwaiter().GetResult());
        File.WriteAllText(Path.Combine(outputRoot, "pilot.stderr.log"), stderr.GetAwaiter().GetResult());

        Check(completed && !cancelled && process.ExitCode == 0, "Empirical pilot failed, timed out or was cancelled; inspect pilot stdout/stderr and stage logs.");
        return Evaluate(protocolPath, Path.Combine(outputRoot, "predictions.csv"), Path.Combine(outputRoot, "fit_artifact.json"), outputRoot);
    }
}

