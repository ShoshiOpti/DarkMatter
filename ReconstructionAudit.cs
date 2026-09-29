using System.Globalization;
using Result = System.Collections.Generic.Dictionary<string, object?>;
using Row = System.Collections.Generic.Dictionary<string, string>;

namespace DarkUniverse;

/// <summary>
/// Audits the distinction between a norm-preserving field reconstruction and
/// an actual projection that discards density. This is an algebraic audit of
/// saved fits, not a new field solution, fit, or observational validation.
/// </summary>
public static class ReconstructionAudit
{
    sealed record Point(string Galaxy, int Index, double Radius, double CatalogueRadius,
        bool Outer, bool Active, double Density, double CoreFraction, double H,
        double OriginalCore, double TargetMass, double ChargedMass, double TargetSpeed,
        double ChargedSpeed);

    public static Result Run(string dataRoot, string outputRoot)
    {
        var safe = ExecutionSafety.ValidateRoots(dataRoot, outputRoot);
        dataRoot = safe.DataRoot;
        string output = Path.Combine(safe.OutputRoot, "reconstruction_audit");
        Directory.CreateDirectory(output);
        string pointsPath = Path.Combine(dataRoot, "publication", "diagnostics", "point_requirements.csv");
        string galaxiesPath = Path.Combine(dataRoot, "publication", "diagnostics", "galaxy_requirements.csv");
        string predictionsPath = Path.Combine(dataRoot, "publication", "statistics", "charged", "predicted_rows.csv");
        string frozenCorePath = Path.Combine(dataRoot, "publication", "diagnostics", "reconstruction", "frozen_core.csv");
        string envelopeFitsPath = Path.Combine(dataRoot, "publication", "diagnostics", "reconstruction", "envelope_fits.csv");
        var hashes = new[] { pointsPath, galaxiesPath, predictionsPath, frozenCorePath, envelopeFitsPath }.ToDictionary(p => p, Data.FileSha256);
        var inputGalaxies = Csv.Read(galaxiesPath);
        var active = new Dictionary<string, bool>(StringComparer.Ordinal);
        foreach (var row in inputGalaxies)
            Require(active.TryAdd(row["galaxy"], Flag(row, "envelope_active")), "Duplicate galaxy requirements row.");

        var observations = Sparc.Load(dataRoot);
        var raw = Csv.Read(pointsPath);
        var points = new Dictionary<(string Galaxy, int Index), Point>();
        int checks = 0;
        double maxFractionError = 0, maxNormError = 0, maxDensityClosure = 0, maxOverlapError = 0;
        var rows = new List<Result>();
        foreach (var row in raw)
        {
            string galaxy = row["galaxy"];
            int index = int.Parse(row["index"], CultureInfo.InvariantCulture);
            Require(active.ContainsKey(galaxy), "Diagnostic galaxy lacks an active-envelope flag.");
            var p = new Point(galaxy, index, N(row, "R_kpc"), N(row, "r_catalogue_kpc"),
                Flag(row, "is_outer"), active[galaxy], N(row, "target_density_msun_kpc3"),
                N(row, "target_core_fraction"), N(row, "projected_orthogonal_density_fraction"),
                N(row, "original_core_density_msun_kpc3"), N(row, "target_enclosed_mass_msun"),
                N(row, "charged_enclosed_mass_msun"), N(row, "target_velocity_physical_kms"),
                N(row, "charged_velocity_physical_kms"));
            Require(points.TryAdd((galaxy, index), p), "Duplicate diagnostic observation.");
            Require(p.Radius > 0 && p.CatalogueRadius > 0 && p.Density > 0 && p.OriginalCore >= 0,
                "Invalid diagnostic radius or density.");
            Require(p.H >= 0 && p.H <= 1 && p.CoreFraction >= 0 && p.CoreFraction <= 1,
                "Projection fractions must lie in [0,1].");
            Require(observations.Points.TryGetValue(galaxy, out var original) && index >= 0 && index < original.Length,
                "Diagnostic observation is absent from the original radial data.");
            var observed = original![index];
            Close(p.CatalogueRadius, observed.RadiusKpc, "Catalogue radius changed.");
            int innerCount = original.Length - Math.Max(2, (int)Math.Ceiling(.2 * original.Length));
            Require(p.Outer == (index >= innerCount), "Frozen radial split changed.");

            // h is a fraction in the fitted core/envelope projection basis.
            // It is not the microscopic parent-field equilibrium composition f*.
            double retainedFraction = 1 - p.H;
            double retainedDensity = p.Density * retainedFraction;
            double orthogonalDensity = p.Density * p.H;
            double archivedCoreDensity = p.Density * p.CoreFraction;
            double retainedAmplitude = Math.Sqrt(retainedFraction);
            double orthogonalAmplitude = Math.Sqrt(p.H);
            double norm = retainedAmplitude * retainedAmplitude + orthogonalAmplitude * orthogonalAmplitude;
            double reconstructedDensity = p.Density * norm;
            double fractionError = Math.Abs(retainedFraction - p.CoreFraction);
            double densityError = Math.Abs(retainedDensity + orthogonalDensity - p.Density) / p.Density;
            double overlap = N(row, "projected_scalar_overlap");
            double overlapError = Math.Abs(overlap * overlap - retainedFraction);
            Require(fractionError <= 2e-14 && Math.Abs(norm - 1) <= 2e-14 && densityError <= 2e-14 && overlapError <= 2e-14,
                "Projection identities fail beyond floating-point tolerance.");
            checks += 4;
            maxFractionError = Math.Max(maxFractionError, fractionError);
            maxNormError = Math.Max(maxNormError, Math.Abs(norm - 1));
            maxDensityClosure = Math.Max(maxDensityClosure, densityError);
            maxOverlapError = Math.Max(maxOverlapError, overlapError);
            rows.Add(new()
            {
                ["galaxy"] = galaxy, ["index"] = index, ["radius_physical_kpc"] = p.Radius,
                ["radius_catalogue_kpc"] = p.CatalogueRadius, ["is_outer"] = p.Outer,
                ["active_envelope"] = p.Active, ["observed_catalogue_kms"] = observed.ObservedSpeedKms,
                ["sigma_catalogue_kms"] = observed.ErrorSpeedKms,
                ["target_total_density_msun_kpc3"] = p.Density,
                ["archived_target_core_fraction"] = p.CoreFraction,
                ["discarded_projection_density_fraction_h"] = p.H,
                ["retained_projection_density_fraction_one_minus_h"] = retainedFraction,
                ["retained_projection_density_msun_kpc3"] = retainedDensity,
                ["orthogonal_projection_density_msun_kpc3"] = orthogonalDensity,
                ["archived_target_core_density_from_fraction_msun_kpc3"] = archivedCoreDensity,
                ["density_sum_relative_error"] = densityError, ["fraction_closure_absolute_error"] = fractionError,
                ["normalized_two_component_norm_squared"] = norm,
                ["norm_reconstructed_total_density_msun_kpc3"] = reconstructedDensity,
                ["norm_reconstruction_relative_roundoff"] = (reconstructedDensity - p.Density) / p.Density,
                ["same_total_density_added_leading_mass_msun"] = 0d,
                ["same_total_density_added_leading_speed_squared_kms2"] = 0d,
                ["original_core_density_msun_kpc3"] = p.OriginalCore,
                ["original_core_minus_retained_projection_density_msun_kpc3"] = p.OriginalCore - retainedDensity,
                ["original_core_minus_retained_projection_over_target_density"] = (p.OriginalCore - retainedDensity) / p.Density,
                ["retained_projection_to_original_core_density"] = p.OriginalCore > 0 ? retainedDensity / p.OriginalCore : null,
                ["saved_target_enclosed_mass_msun"] = p.TargetMass,
                ["saved_charged_enclosed_mass_msun"] = p.ChargedMass,
                ["saved_target_velocity_physical_kms"] = p.TargetSpeed,
                ["saved_charged_velocity_physical_kms"] = p.ChargedSpeed
            });
        }
        Require(points.Count == 3034 && points.Values.Count(p => p.Outer) == 659 && active.Count == 131,
            "The audit requires the frozen 131-galaxy, 3034-observation, 659-outer sample.");
        Require(points.Values.Select(p => p.Galaxy).Distinct().Count() == active.Count,
            "Point and galaxy samples differ.");

        var galaxyRows = new List<Result>();
        foreach (var group in rows.GroupBy(r => (string)r["galaxy"]!).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            var outer = group.Where(r => (bool)r["is_outer"]!).ToArray();
            Require(outer.Length > 0, "A selected galaxy has no outer observations.");
            galaxyRows.Add(new()
            {
                ["galaxy"] = group.Key, ["points"] = group.Count(), ["outer_points"] = outer.Length,
                ["active_envelope"] = active[group.Key],
                ["outer_median_discarded_projection_density_fraction"] = Data.Median(outer.Select(r => V(r, "discarded_projection_density_fraction_h"))),
                ["outer_median_retained_projection_density_fraction"] = Data.Median(outer.Select(r => V(r, "retained_projection_density_fraction_one_minus_h"))),
                ["outer_median_original_core_minus_retained_over_target_density"] = Data.Median(outer.Select(r => V(r, "original_core_minus_retained_projection_over_target_density"))),
                ["outer_median_retained_to_original_core_density"] = MedianNullable(outer.Select(r => r["retained_projection_to_original_core_density"])),
                ["outer_zero_original_core_rows"] = outer.Count(r => V(r, "original_core_density_msun_kpc3") == 0)
            });
        }

        // This identity propagates the unchanged leading Poisson source. It is
        // not a recalculated velocity for a newly excited/stressed field state.
        var predictions = Csv.Read(predictionsPath);
        var predictionKeys = new HashSet<(string, int, string)>();
        var invariantRows = new List<Result>();
        var modelRows = new List<Result>();
        foreach (var row in predictions)
        {
            string galaxy = row["galaxy"], model = row["model"];
            int index = int.Parse(row["index"], CultureInfo.InvariantCulture);
            Require(predictionKeys.Add((galaxy, index, model)), "Duplicate model prediction.");
            Require(points.TryGetValue((galaxy, index), out var p), "A prediction has no diagnostic observation.");
            double radius = N(row, "radius"), observed = N(row, "observed"), sigma = N(row, "sigma");
            var original = observations.Points[galaxy][index];
            Close(radius, original.RadiusKpc, "Prediction radius changed.");
            Close(observed, original.ObservedSpeedKms, "Prediction observation changed.");
            Close(sigma, original.ErrorSpeedKms, "Prediction error changed.");
            bool outer = Flag(row, "is_outer"), failed = Flag(row, "fit_failed");
            Require(outer == p!.Outer, "Prediction and diagnostic splits differ.");
            double before = Csv.Number(row, "predicted"), after = before;
            bool exact = BitConverter.DoubleToInt64Bits(before) == BitConverter.DoubleToInt64Bits(after);
            Require(exact, "Algebraic coordinate identity changed a prediction.");
            invariantRows.Add(new()
            {
                ["galaxy"] = galaxy, ["index"] = index, ["model"] = model, ["radius_catalogue_kpc"] = radius,
                ["is_outer"] = outer, ["fit_failed_preserved"] = failed,
                ["observed_catalogue_kms"] = observed, ["sigma_catalogue_kms"] = sigma,
                ["saved_predicted_catalogue_kms"] = before,
                ["predicted_after_same_source_identity_catalogue_kms"] = after,
                ["identity_speed_change_kms"] = double.IsFinite(before) ? 0d : null,
                ["finite_nonnegative_prediction"] = double.IsFinite(before) && before >= 0,
                ["bitwise_prediction_preserved"] = exact
            });
        }
        foreach (var model in invariantRows.GroupBy(r => (string)r["model"]!).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            Require(model.Count() == points.Count, "A saved model lacks observations.");
            modelRows.Add(new()
            {
                ["model"] = model.Key, ["rows"] = model.Count(), ["outer_rows"] = model.Count(r => (bool)r["is_outer"]!),
                ["failed_rows_preserved"] = model.Count(r => (bool)r["fit_failed_preserved"]!),
                ["nonfinite_or_negative_predictions_preserved"] = model.Count(r => !(bool)r["finite_nonnegative_prediction"]!),
                ["maximum_finite_identity_speed_change_kms"] = 0d, ["all_predictions_bitwise_preserved"] = true
            });
        }
        Require(modelRows.Count == 7 && predictions.Count == 21238, "Saved prediction model/sample coverage changed.");
        Csv.Write(Path.Combine(output, "density_decomposition.csv"), rows);
        Csv.Write(Path.Combine(output, "galaxy_projection_summary.csv"), galaxyRows);
        Csv.Write(Path.Combine(output, "saved_prediction_identity.csv"), invariantRows);
        Csv.Write(Path.Combine(output, "model_identity_summary.csv"), modelRows);
        Plot(output, galaxyRows);
        var forceBridge = ForceBridge(points, predictions, frozenCorePath, envelopeFitsPath, output);
        foreach (var (path, hash) in hashes) Require(Data.FileSha256(path) == hash, "An immutable audit input changed.");
        var all = galaxyRows.Select(r => V(r, "outer_median_discarded_projection_density_fraction")).ToArray();
        var activeOnly = galaxyRows.Where(r => (bool)r["active_envelope"]!).Select(r => V(r, "outer_median_discarded_projection_density_fraction")).ToArray();
        var report = new Result
        {
            ["schema_version"] = 1, ["status"] = "passed", ["algebraic_identity_checks"] = checks,
            ["galaxies"] = galaxyRows.Count, ["active_envelope_galaxies"] = activeOnly.Length,
            ["all_rows"] = rows.Count, ["outer_rows"] = points.Values.Count(p => p.Outer),
            ["prediction_rows"] = predictions.Count, ["prediction_models"] = modelRows,
            ["force_bridge"] = forceBridge,
            ["equal_galaxy_outer_median_discarded_density_fraction_all"] = Data.Median(all),
            ["equal_galaxy_outer_median_discarded_density_fraction_active"] = Data.Median(activeOnly),
            ["equal_galaxy_mean_of_outer_median_discarded_density_fraction_all"] = all.Average(),
            ["equal_galaxy_mean_of_outer_median_discarded_density_fraction_active"] = activeOnly.Average(),
            ["maximum_fraction_closure_absolute_error"] = maxFractionError,
            ["maximum_normalized_amplitude_norm_error"] = maxNormError,
            ["maximum_density_sum_relative_error"] = maxDensityClosure,
            ["maximum_scalar_overlap_squared_error"] = maxOverlapError,
            ["maximum_absolute_original_core_minus_retained_over_target_density"] = rows.Max(r => Math.Abs(V(r, "original_core_minus_retained_projection_over_target_density"))),
            ["source_sha256"] = hashes.ToDictionary(p => Path.GetRelativePath(dataRoot, p.Key).Replace('\\', '/'), p => p.Value),
            ["original_catalogue_sha256"] = observations.CatalogueSha256,
            ["original_mass_models_sha256"] = observations.MassModelsSha256,
            ["fraction_definition"] = "h = rho_orthogonal/rho_target = 1 - target_core_fraction in the archived fitted core/envelope projection basis; h is not the microscopic equilibrium component fraction f*.",
            ["summary_definition"] = "First take the median h among every galaxy's outer observations, without derivative masks; then give each galaxy equal weight. The median of galaxy medians is neither an enclosed-mass fraction nor a pooled-point average.",
            ["original_core_caveat"] = "original_core_density is the separate saved core quantity. The force bridge independently verifies rho_target*target_core_fraction = core_multiplier*original_core_density. This rescaled density is not generally the raw original core.",
            ["same_source_identity"] = "For a norm-preserving reconstruction at every radius, delta rho_rest = 0 and hence delta M_rest(<R) = 0, so delta V_halo^2 = G delta M/R = 0 at leading Newtonian order with the same boundary potential and baryons. Propagating this identity leaves saved model predictions and their failure flags unchanged by construction.",
            ["stress_caveat"] = "An exact invertible field-coordinate change preserves physical stress. A physically different normalized internal configuration can carry additional gradient/contact stress while preserving leading rest density. This audit does not solve that configuration, its stationarity, backreaction or higher-order gravitational source.",
            ["projection_caveat"] = "An actual one-component projection discards rho_target*h locally. Finite sampled radii alone do not determine the projected enclosed mass or a new rotation curve; neither is fabricated here.",
            ["scope"] = "Algebraic retrospective audit of fitted density components and saved predictions. No fit, dynamical reconstruction, halo occupation law, new statistical test or independent observational validation.",
            ["files"] = new[] { "density_decomposition.csv", "galaxy_projection_summary.csv", "saved_prediction_identity.csv", "model_identity_summary.csv", "projection_audit.svg", "projection_audit.series.csv", "force_decomposition.csv", "galaxy_force_summary.csv", "force_bridge.json", "force_projection_bridge.svg", "force_projection_bridge.series.csv" }
        };
        Data.SaveJson(Path.Combine(output, "report.json"), report);
        File.WriteAllText(Path.Combine(output, "README.md"),
            "# Reconstruction and projection audit\n\n" + report["scope"] + "\n\n" + report["fraction_definition"] + "\n\n" +
            report["summary_definition"] + "\n\n" + report["same_source_identity"] + "\n\n" + report["stress_caveat"] + "\n\n" +
            report["projection_caveat"] + "\n\n" + report["original_core_caveat"] + "\n\nThe force bridge uses the existing selected envelope fit and independently reproduces its saved curve. Orthogonal restoration is measured relative to the rescaled projected core. Envelope minus original core also contains the fitted core-amplitude change. Local density fractions differ from enclosed halo-force fractions; no new model curve or fit is produced." + "\n\nThe SVG uses only finite galaxy summaries and explicit axis limits. All original prediction failures remain in the identity table.\n");
        Console.WriteLine($"Reconstruction audit: {rows.Count} density rows; outer median h across galaxies {Data.Median(all):G6}, active galaxies {Data.Median(activeOnly):G6}; {predictions.Count} predictions unchanged by identity.");
        return report;
    }

    static Result ForceBridge(Dictionary<(string Galaxy, int Index), Point> points, List<Row> predictions,
        string frozenCorePath, string envelopeFitsPath, string output)
    {
        var fits = new Dictionary<string, Row>(StringComparer.Ordinal);
        foreach (var row in Csv.Read(envelopeFitsPath).Where(r => r["model"] == "envelope_a1_q20"))
            Require(fits.TryAdd(row["galaxy"], row), "Duplicate selected envelope fit.");
        Require(fits.Count == 131, "Selected envelope fit coverage changed.");
        var saved = predictions.Where(r => r["model"] is "envelope" or "original_core")
            .ToDictionary(r => (r["galaxy"], int.Parse(r["index"], CultureInfo.InvariantCulture), r["model"]));
        var forceRows = new List<Result>();
        var summaries = new List<Result>();
        double maxSpeedError = 0, maxOriginalSpeedError = 0, maxForceIdentityError = 0,
            maxDensityRescalingError = 0, maxPhysicalSpeedError = 0;
        const double signThreshold = 1e-8;
        foreach (var group in Csv.Read(frozenCorePath).GroupBy(r => r["galaxy"]).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            Require(fits.TryGetValue(group.Key, out var fit), "Frozen core has no selected envelope fit.");
            double alpha = N(fit!, "core_multiplier"), amplitude = N(fit, "B_catalogue"),
                a = N(fit, "a_catalogue"), b = N(fit, "b_catalogue");
            Require(alpha >= 0 && amplitude >= 0 && a > 0 && b > a, "Invalid selected envelope parameters.");
            bool active = amplitude > 1e-6;
            var local = new List<Result>();
            var originalRows = group.OrderBy(r => N(r, "r")).ToArray();
            for (int index = 0; index < originalRows.Length; index++)
            {
                var core = originalRows[index];
                Require(points.TryGetValue((group.Key, index), out var p), "Frozen core observation has no diagnostic row.");
                double radius = N(core, "r"), baryon = N(core, "baryon_v2"), coreV2 = N(core, "halo_v2"),
                    tilt = N(core, "inclination_factor"), distance = N(core, "distance_ratio");
                Close(radius, p!.CatalogueRadius, "Force bridge catalogue radius differs.");
                Close(distance * radius, p.Radius, "Force bridge physical radius differs.");
                Close(tilt, N(fit, "inclination_factor"), "Force bridge inclination calibration differs.");
                Close(distance, N(fit, "distance_ratio"), "Force bridge distance calibration differs.");
                Require(tilt > 0 && distance > 0 && coreV2 >= 0 && active == p.Active && Flag(core, "is_outer") == p.Outer,
                    "Force bridge physical calibration, activity or radial split differs.");
                double envelopeV2 = amplitude * (b * Math.Atan(radius / b) - a * Math.Atan(radius / a)) / radius;
                double retainedCoreV2 = alpha * coreV2;
                double projectedV2 = baryon + retainedCoreV2;
                double totalV2 = projectedV2 + envelopeV2;
                double originalV2 = baryon + coreV2;
                Require(envelopeV2 >= 0 && totalV2 > 0 && originalV2 > 0 && retainedCoreV2 + envelopeV2 > 0,
                    "Invalid squared force decomposition.");
                double envelopeSpeed = N(saved[(group.Key, index, "envelope")], "predicted");
                double originalSpeed = N(saved[(group.Key, index, "original_core")], "predicted");
                double reconstructedSpeed = Math.Sqrt(totalV2), reconstructedOriginal = Math.Sqrt(originalV2);
                double speedError = Math.Abs(reconstructedSpeed - envelopeSpeed);
                double originalError = Math.Abs(reconstructedOriginal - originalSpeed);
                Require(speedError <= 2e-9 && originalError <= 2e-9, "Exact force bridge does not reproduce saved curves.");
                double coreChange = (alpha - 1) * coreV2;
                double netChange = envelopeSpeed * envelopeSpeed - originalSpeed * originalSpeed;
                double identityError = Math.Abs(netChange - (coreChange + envelopeV2));
                Require(identityError <= 2e-10 * Math.Max(1, totalV2), "Envelope minus original-core force identity fails.");
                double rescalingError = Math.Abs(p.Density * p.CoreFraction - alpha * p.OriginalCore) / p.Density;
                Require(rescalingError <= 2e-13, "Projected target core is not the selected rescaled original core.");
                double physicalError = Math.Abs(reconstructedSpeed / tilt - p.TargetSpeed);
                Require(physicalError <= 2e-8, "Physical and catalogue target speeds disagree.");
                maxSpeedError = Math.Max(maxSpeedError, speedError);
                maxOriginalSpeedError = Math.Max(maxOriginalSpeedError, originalError);
                maxForceIdentityError = Math.Max(maxForceIdentityError, identityError);
                maxDensityRescalingError = Math.Max(maxDensityRescalingError, rescalingError);
                maxPhysicalSpeedError = Math.Max(maxPhysicalSpeedError, physicalError);
                var row = new Result
                {
                    ["galaxy"] = group.Key, ["index"] = index, ["radius_catalogue_kpc"] = radius,
                    ["radius_physical_kpc"] = p.Radius, ["is_outer"] = p.Outer, ["active_envelope"] = active,
                    ["core_multiplier_alpha"] = alpha, ["baryon_v2_catalogue_kms2"] = baryon,
                    ["original_core_v2_catalogue_kms2"] = coreV2,
                    ["rescaled_projected_core_v2_catalogue_kms2"] = retainedCoreV2,
                    ["orthogonal_envelope_v2_catalogue_kms2"] = envelopeV2,
                    ["projected_core_plus_baryon_v2_catalogue_kms2"] = projectedV2,
                    ["saved_envelope_speed_catalogue_kms"] = envelopeSpeed,
                    ["reproduced_envelope_speed_catalogue_kms"] = reconstructedSpeed,
                    ["saved_original_speed_catalogue_kms"] = originalSpeed,
                    ["reproduced_original_speed_catalogue_kms"] = reconstructedOriginal,
                    ["core_rescaling_force_change_catalogue_kms2"] = coreChange,
                    ["full_saved_envelope_minus_original_force_catalogue_kms2"] = netChange,
                    ["force_decomposition_absolute_error_kms2"] = identityError,
                    ["negative_full_force_change"] = netChange < -signThreshold,
                    ["positive_full_force_change"] = netChange > signThreshold,
                    ["orthogonal_force_fraction_of_total"] = envelopeV2 / totalV2,
                    ["orthogonal_force_fraction_of_halo"] = envelopeV2 / (retainedCoreV2 + envelopeV2),
                    ["local_discarded_density_fraction_h"] = p.H,
                    ["core_density_rescaling_relative_error"] = rescalingError,
                    ["projected_total_force_admissible"] = projectedV2 >= 0
                };
                local.Add(row); forceRows.Add(row);
            }
            var outer = local.Where(r => (bool)r["is_outer"]!).ToArray();
            int negative = outer.Count(r => (bool)r["negative_full_force_change"]!),
                positive = outer.Count(r => (bool)r["positive_full_force_change"]!);
            summaries.Add(new()
            {
                ["galaxy"] = group.Key, ["active_envelope"] = active, ["core_multiplier_alpha"] = alpha,
                ["outer_points"] = outer.Length, ["negative_outer_force_change_points"] = negative,
                ["positive_outer_force_change_points"] = positive,
                ["neutral_outer_force_change_points"] = outer.Length - negative - positive,
                ["all_outer_force_changes_negative"] = negative == outer.Length,
                ["all_outer_force_changes_positive"] = positive == outer.Length,
                ["both_outer_force_change_signs"] = negative > 0 && positive > 0,
                ["outer_median_local_discarded_density_fraction"] = Data.Median(outer.Select(r => V(r, "local_discarded_density_fraction_h"))),
                ["outer_median_orthogonal_force_fraction_of_total"] = Data.Median(outer.Select(r => V(r, "orthogonal_force_fraction_of_total"))),
                ["outer_median_orthogonal_force_fraction_of_halo"] = Data.Median(outer.Select(r => V(r, "orthogonal_force_fraction_of_halo")))
            });
        }
        Require(forceRows.Count == 3034 && summaries.Count == 131, "Force bridge observation coverage changed.");
        var activeGalaxies = summaries.Where(r => (bool)r["active_envelope"]!).ToArray();
        var activeOuter = forceRows.Where(r => (bool)r["active_envelope"]! && (bool)r["is_outer"]!).ToArray();
        Result Aggregate(IEnumerable<Result> values)
        {
            var sample = values.ToArray();
            return new()
            {
                ["galaxies"] = sample.Length,
                ["median_outer_median_local_discarded_density_fraction"] = Data.Median(sample.Select(r => V(r, "outer_median_local_discarded_density_fraction"))),
                ["median_outer_median_orthogonal_force_fraction_of_total"] = Data.Median(sample.Select(r => V(r, "outer_median_orthogonal_force_fraction_of_total"))),
                ["median_outer_median_orthogonal_force_fraction_of_halo"] = Data.Median(sample.Select(r => V(r, "outer_median_orthogonal_force_fraction_of_halo"))),
                ["median_core_multiplier"] = Data.Median(sample.Select(r => V(r, "core_multiplier_alpha")))
            };
        }
        var result = new Result
        {
            ["status"] = "passed", ["selected_model"] = "envelope_a1_q20", ["activity_definition"] = "B_catalogue > 1e-6",
            ["all_galaxies"] = Aggregate(summaries), ["active_galaxies"] = Aggregate(activeGalaxies),
            ["active_outer_points"] = activeOuter.Length,
            ["active_outer_negative_full_force_changes"] = activeOuter.Count(r => (bool)r["negative_full_force_change"]!),
            ["active_outer_positive_full_force_changes"] = activeOuter.Count(r => (bool)r["positive_full_force_change"]!),
            ["active_galaxies_with_any_negative_outer_change"] = activeGalaxies.Count(r => (int)r["negative_outer_force_change_points"]! > 0),
            ["active_galaxies_all_outer_negative"] = activeGalaxies.Count(r => (bool)r["all_outer_force_changes_negative"]!),
            ["active_galaxies_all_outer_positive"] = activeGalaxies.Count(r => (bool)r["all_outer_force_changes_positive"]!),
            ["active_galaxies_mixed_outer_signs"] = activeGalaxies.Count(r => (bool)r["both_outer_force_change_signs"]!),
            ["sign_threshold_catalogue_kms2"] = signThreshold,
            ["maximum_saved_envelope_speed_reproduction_error_kms"] = maxSpeedError,
            ["maximum_saved_original_speed_reproduction_error_kms"] = maxOriginalSpeedError,
            ["maximum_physical_speed_reproduction_error_kms"] = maxPhysicalSpeedError,
            ["maximum_force_identity_absolute_error_kms2"] = maxForceIdentityError,
            ["maximum_density_rescaling_relative_error"] = maxDensityRescalingError,
            ["projected_core_plus_baryons_negative_v2_rows"] = forceRows.Count(r => !(bool)r["projected_total_force_admissible"]!),
            ["negative_projected_total_force_rows"] = forceRows.Where(r => !(bool)r["projected_total_force_admissible"]!).Select(r => new Result
            {
                ["galaxy"] = r["galaxy"], ["index"] = r["index"], ["radius_catalogue_kpc"] = r["radius_catalogue_kpc"],
                ["is_outer"] = r["is_outer"], ["projected_core_plus_baryon_v2_catalogue_kms2"] = r["projected_core_plus_baryon_v2_catalogue_kms2"]
            }).ToArray(),
            ["projected_force_caveat"] = "The projected core plus signed baryonic force need not admit a real circular speed at every radius. Negative squared-force rows remain signed and are not square-rooted, clipped or omitted.",
            ["force_identity"] = "V_env^2 = V_b^2 + alpha_c*V_core^2 + V_e^2; V_env^2 - V_original^2 = (alpha_c-1)*V_core^2 + V_e^2. Orthogonal restoration from the rescaled projected core contributes V_e^2; it is not the full envelope-minus-original change.",
            ["envelope_force_formula"] = "V_e,cat^2 = B_cat*[b_cat*atan(r_cat/b_cat)-a_cat*atan(r_cat/a_cat)]/r_cat; physical V^2=V_cat^2/inclination_factor^2 and R=distance_ratio*r_cat.",
            ["fraction_interpretation"] = "h is a local density fraction. V_e^2/(alpha_c V_core^2+V_e^2) is the envelope fraction of the spherical halo's enclosed mass/force. V_e^2/V_env^2 includes baryonic force in its denominator. Their per-galaxy outer medians are different summaries.",
            ["double_counting_warning"] = "Adding V_e^2 to the saved envelope or to a scalar with the same total density counts the fitted envelope twice. No such counterfactual is treated as a fitted model or assigned significance.",
            ["scope"] = "Exact bookkeeping of an existing fitted analytic envelope, not a derived occupation, equilibrium or new observational fit. Negative net changes arise from fitted core rescaling, not negative envelope density."
        };
        Csv.Write(Path.Combine(output, "force_decomposition.csv"), forceRows);
        Csv.Write(Path.Combine(output, "galaxy_force_summary.csv"), summaries);
        Data.SaveJson(Path.Combine(output, "force_bridge.json"), result);
        var plot = new PlotDocument("Restoring the projected component is not the full curve difference", 2, 1)
        { Scope = "103 active fitted envelopes; same outer observations. Different fractions have different denominators." };
        var fractions = new PlotPanel("Local density and integrated force differ", "Galaxy median outer fraction", "Fraction of active galaxies")
        { XMin = -.02, XMax = 1.02, YMin = -.02, YMax = 1.02,
            XTicks = FractionTicks(), YTicks = FractionTicks(), Note = "Each curve gives equal weight to the same 103 galaxies." };
        fractions.Ecdf("Local discarded density h", activeGalaxies.Select(r => V(r, "outer_median_local_discarded_density_fraction")), "#24658c", true);
        fractions.Ecdf("Envelope / enclosed halo force", activeGalaxies.Select(r => V(r, "outer_median_orthogonal_force_fraction_of_halo")), "#b44e34", true);
        fractions.Ecdf("Envelope / total force", activeGalaxies.Select(r => V(r, "outer_median_orthogonal_force_fraction_of_total")), "#4c865c", true);
        plot.Panels.Add(fractions);
        double[] counts = [activeGalaxies.Count(r => (bool)r["all_outer_force_changes_positive"]!),
            activeGalaxies.Count(r => (bool)r["both_outer_force_change_signs"]!),
            activeGalaxies.Count(r => (bool)r["all_outer_force_changes_negative"]!)];
        var signs = new PlotPanel("Envelope minus original core can lower force", "Signs across a galaxy's outer observations", "Active galaxies")
        {
            XMin = .4, XMax = 3.6, YMin = 0, YMax = Math.Ceiling((counts.Max() + 5) / 10) * 10,
            YTicks = Enumerable.Range(0, (int)Math.Ceiling((counts.Max() + 5) / 10) + 1).ToDictionary(i => 10d * i, i => (10 * i).ToString(CultureInfo.InvariantCulture)),
            XTicks = new() { [1] = "All positive", [2] = "Mixed", [3] = "All negative" },
            Note = "Includes fitted core rescaling; threshold 1e-8 (km/s)^2."
        };
        signs.Add("Full saved force change", [1, 2, 3], counts, "#24658c").Bars = true;
        plot.Panels.Add(signs);
        plot.Save(Path.Combine(output, "force_projection_bridge.svg"));
        return result;
    }

    static void Plot(string output, List<Result> galaxies)
    {
        var document = new PlotDocument("Projection loss in the existing fitted density", 2, 1)
        {
            Scope = "Saved fits; h is local projected density, not an observed mode population. No new fit."
        };
        var ecdf = new PlotPanel("Each galaxy contributes one outer median", "Outer median discarded density fraction h", "Fraction of galaxies")
        { XMin = -.02, XMax = 1.02, YMin = -.02, YMax = 1.02,
            XTicks = FractionTicks(), YTicks = FractionTicks(), Note = "Outer rows retain their original split; no derivative mask." };
        ecdf.Ecdf("All 131 galaxies", galaxies.Select(r => V(r, "outer_median_discarded_projection_density_fraction")), "#24658c", true);
        ecdf.Ecdf("103 active envelopes", galaxies.Where(r => (bool)r["active_envelope"]!).Select(r => V(r, "outer_median_discarded_projection_density_fraction")), "#b44e34", true);
        document.Panels.Add(ecdf);
        string xKey = "outer_median_discarded_projection_density_fraction", yKey = "outer_median_original_core_minus_retained_over_target_density";
        double[] values = galaxies.Select(r => V(r, yKey)).ToArray();
        Require(values.All(double.IsFinite), "Plot summary contains a nonfinite value.");
        double low = Math.Min(0, values.Min()), high = Math.Max(0, values.Max());
        double padding = Math.Max(.05, .1 * (high - low));
        var comparison = new PlotPanel("Saved original core is a separate quantity", "Outer median discarded density fraction h", "Median (original core - retained core) / total")
        { XMin = -.02, XMax = 1.02, XTicks = FractionTicks(), YMin = low - padding, YMax = high + padding, Note = "Retained core = target density times (1-h); scaling can differ." };
        comparison.Add("Equal core densities", [0, 1], [0, 0], "#888888").Dash = "4 4";
        foreach (bool isActive in new[] { false, true })
        {
            var sample = galaxies.Where(r => (bool)r["active_envelope"]! == isActive).ToArray();
            comparison.Add(isActive ? "Active envelope" : "Inactive envelope", sample.Select(r => V(r, xKey)).ToArray(),
                sample.Select(r => V(r, yKey)).ToArray(), isActive ? "#b44e34" : "#24658c", true);
        }
        document.Panels.Add(comparison);
        document.Save(Path.Combine(output, "projection_audit.svg"));
    }

    static Dictionary<double, string> FractionTicks() => Enumerable.Range(0, 6).ToDictionary(i => i / 5d, i => (i / 5d).ToString("0.0", CultureInfo.InvariantCulture));
    static void Require([System.Diagnostics.CodeAnalysis.DoesNotReturnIf(false)] bool condition, string message)
    { if (!condition) throw new InvalidDataException("Reconstruction audit: " + message); }
    static void Close(double a, double b, string message) => Require(Math.Abs(a - b) <= 1e-10 * Math.Max(1, Math.Abs(b)), message);
    static double N(Row row, string key)
    { double value = Csv.Number(row, key); Require(double.IsFinite(value), "Nonfinite required input: " + key); return value; }
    static double V(Result row, string key) => (double)row[key]!;
    static double? MedianNullable(IEnumerable<object?> values)
    { double[] finite = values.OfType<double>().Where(double.IsFinite).ToArray(); return finite.Length == 0 ? null : Data.Median(finite); }
    static bool Flag(Row row, string key) => Csv.Text(row, key).Trim().ToLowerInvariant() switch
    {
        "true" or "1" or "1.0" => true, "false" or "0" or "0.0" => false,
        _ => throw new InvalidDataException("Invalid boolean in reconstruction audit: " + key)
    };
}
