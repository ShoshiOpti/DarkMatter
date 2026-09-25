using MathNet.Numerics.Statistics;
using System.Globalization;
using System.IO.Compression;
using System.Text;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

public static partial class PointwiseStatistics
{
    const string Analysis = "pointwise_comparison_2026_09_25";
    const string InputHash = "a057d37495824894874e3bcd0127a0752df145851cafd273fb9f378af38ed9a6";
    const int Draws = 99999;
    const ulong Seed = 2026092501;
    static readonly string[] Models = ["extended", "scalar_core", "compact_plummer", "nfw", "baryon_only", "profiled_baryons", "adaptive_primary"];
    static readonly string[] Comparators = ["scalar_core", "compact_plummer", "nfw", "baryon_only"];
    sealed record Point(string Galaxy, int Index, double Radius, double Observed, double Sigma, double[] Z, bool[] Failed);
    sealed record Descriptor(string Score, string Comparator, double? Rho = null);

    public static Result Run(string dataRoot, string outputRoot)
    {
        string numerics = Directory.Exists(Path.Combine(dataRoot, "numerics")) ? Path.Combine(dataRoot, "numerics") : dataRoot;
        string source = Path.Combine(numerics, Analysis), destination = Path.Combine(outputRoot, "pointwise_comparison", "results");
        Directory.CreateDirectory(destination);
        double[] rhos = ReadCorrelationSensitivities(source);
        string input = Path.Combine(source, "data", "pointwise_outer.csv");
        string hash = Data.FileSha256(input);
        Require(hash == InputHash, "Frozen pointwise input SHA-256 mismatch.");
        var points = ReadPoints(input);
        var galaxyGroups = points.GroupBy(p => p.Galaxy).ToArray();
        string[] galaxies = galaxyGroups.Select(g => g.Key).ToArray();
        int galaxyCount = galaxyGroups.Length, comparatorCount = Comparators.Length, pointCount = points.Length;
        int[] counts = galaxyGroups.Select(g => g.Count()).ToArray();
        Require(galaxyCount == 131 && pointCount == 659 && points.Sum(p => p.Failed[5] ? 1 : 0) == 5, "Incorrect frozen sample or failed profile count.");
        var (modelLosses, meanResiduals, residualVariances) = SummarizeGalaxyResiduals(galaxyGroups);
        Csv.Write(Path.Combine(destination, "model_galaxy_squared_losses.csv"), Enumerable.Range(0, galaxyCount).Select(g =>
        {
            var row = new Result { ["galaxy"] = galaxies[g] };
            for (int m = 0; m < Models.Length; m++)
                row[Models[m]] = modelLosses[g][m];
            return row;
        }));
        WritePairedDifferences(destination, galaxyGroups);

        var descriptors = CreateScoreDescriptors(rhos);
        int scoreCount = descriptors.Count;
        var contrasts = CalculateGalaxyContrasts(galaxyGroups, counts, descriptors, modelLosses, meanResiduals, residualVariances);
        double[] scoreMeans = new double[scoreCount], scoreErrors = new double[scoreCount];
        for (int j = 0; j < scoreCount; j++)
            (scoreMeans[j], scoreErrors[j]) = MeanSe(contrasts.Select(row => row[j]).ToArray());
        Console.WriteLine($"Pointwise: {pointCount} outer measurements, {galaxyCount} galaxy clusters; {Draws:N0} pairs and wild bootstrap draws.");
        var (bootstrapMeans, bootstrapTStatistics, pooledBootstrapMeans) = BootstrapPairedGalaxies(contrasts, counts, scoreMeans);

        int[] exceedances = WildBootstrap(contrasts, scoreMeans, scoreErrors);

        // Plus-one probabilities retain the Monte Carlo floor; zero tails remain flagged.
        double[] rawPValues = exceedances.Select(tailCount => (tailCount + 1.0) / (Draws + 1)).ToArray(), holmPValues = Statistics.Holm(rawPValues);
        var results = new List<Result>();
        var influence = new List<Result>();
        for (int j = 0; j < scoreCount; j++)
        {
            var descriptor = descriptors[j];
            double[] galaxyDifferences = contrasts.Select(row => row[j]).ToArray();
            double sumSquares = galaxyDifferences.Sum(x => Square(x - scoreMeans[j])), sum = galaxyDifferences.Sum();
            double[] varianceShares = galaxyDifferences.Select(x => Square(x - scoreMeans[j]) / sumSquares).ToArray(), leaveOneOutMeans = galaxyDifferences.Select(x => (sum - x) / (galaxyCount - 1)).ToArray();
            double[] tQuantiles = QuantilePair(bootstrapTStatistics[j]), meanQuantiles = QuantilePair(bootstrapMeans[j]);

            // This pairs-bootstrap interval does not invert the separate wild-bootstrap test.
            double low = scoreMeans[j] - tQuantiles[1] * scoreErrors[j], high = scoreMeans[j] - tQuantiles[0] * scoreErrors[j];
            var row = new Result
            {
                ["score"] = descriptor.Score,
                ["comparator"] = descriptor.Comparator,
                ["n_galaxies"] = galaxyCount,
                ["n_points"] = pointCount,
                ["effect"] = scoreMeans[j],
                ["cluster_se"] = scoreErrors[j],
                ["bootstrap_t_ci_low"] = low,
                ["bootstrap_t_ci_high"] = high,
                ["percentile_ci_low"] = meanQuantiles[0],
                ["percentile_ci_high"] = meanQuantiles[1],
                ["ci_crosses_zero"] = low <= 0 && high >= 0,
                ["variance_effective_clusters"] = 1 / varianceShares.Sum(Square),
                ["max_variance_share"] = varianceShares.Max(),
                ["max_variance_galaxy"] = galaxies[Array.IndexOf(varianceShares, varianceShares.Max())],
                ["leave_one_out_min"] = leaveOneOutMeans.Min(),
                ["leave_one_out_max"] = leaveOneOutMeans.Max(),
                ["leave_one_out_reversals"] = leaveOneOutMeans.Count(x => Math.Sign(x) != Math.Sign(scoreMeans[j]))
            };
            string[] primaryFields = ["t_statistic", "cluster_t_p_sensitivity", "wild_exceedances", "wild_draws", "p_wild", "p_holm4", "Z_holm4", "raw_MC_95low", "raw_MC_95high", "monte_carlo_floor", "inference"];
            foreach (string key in primaryFields)
                row[key] = null;
            row["rho"] = descriptor.Rho;
            if (j < comparatorCount)
            {
                double t = scoreMeans[j] / scoreErrors[j];
                int tailCount = exceedances[j];
                row["t_statistic"] = t;
                row["cluster_t_p_sensitivity"] = Statistics.RegularizedBeta((galaxyCount - 1.0) / (galaxyCount - 1 + t * t), (galaxyCount - 1.0) / 2, .5);
                row["wild_exceedances"] = tailCount;
                row["wild_draws"] = Draws;
                row["p_wild"] = rawPValues[j];
                row["p_holm4"] = holmPValues[j];
                row["Z_holm4"] = Statistics.NormalMagnitude(holmPValues[j]);
                row["raw_MC_95low"] = tailCount == 0 ? 0 : Statistics.BetaInverse(.025, tailCount, Draws - tailCount + 1);
                row["raw_MC_95high"] = tailCount == Draws ? 1 : Statistics.BetaInverse(.975, tailCount + 1, Draws - tailCount);
                row["monte_carlo_floor"] = tailCount == 0;
                row["inference"] = "nominal_retrospective_fixed_prediction_mean_loss_test";
                for (int g = 0; g < galaxyCount; g++)
                    influence.Add(new Result
                    {
                        ["comparator"] = descriptor.Comparator,
                        ["galaxy"] = galaxies[g],
                        ["delta_galaxy_loss"] = galaxyDifferences[g],
                        ["mean_contribution"] = galaxyDifferences[g] / galaxyCount,
                        ["variance_share"] = varianceShares[g],
                        ["mean_without"] = leaveOneOutMeans[g],
                        ["reverses_direction"] = Math.Sign(leaveOneOutMeans[g]) != Math.Sign(scoreMeans[j])
                    });
            }
            results.Add(row);
        }

        var primary = results.Take(comparatorCount).ToList();
        Csv.Write(Path.Combine(destination, "score_sensitivity.csv"), results);
        Csv.Write(Path.Combine(destination, "primary_cluster_results.csv"), primary);
        Csv.Write(Path.Combine(destination, "galaxy_influence.csv"), influence);
        Csv.Write(Path.Combine(destination, "galaxy_contrasts.csv"), Enumerable.Range(0, galaxyCount).Select(g =>
        {
            var row = new Result { [""] = galaxies[g] };
            for (int j = 0; j < scoreCount; j++)
            {
                var descriptor = descriptors[j];
                row[descriptor.Score + "__" + descriptor.Comparator + (descriptor.Rho is null ? "" : "__rho" + descriptor.Rho.Value.ToString(CultureInfo.InvariantCulture))] = contrasts[g][j];
            }
            return row;
        }));
        var pooledResults = new List<Result>();
        for (int j = 0; j < comparatorCount; j++)
        {
            double[] q = QuantilePair(pooledBootstrapMeans[j]);
            pooledResults.Add(new Result
            {
                ["comparator"] = Comparators[j],
                ["n_galaxies"] = galaxyCount,
                ["n_points"] = pointCount,
                ["effect"] = Enumerable.Range(0, galaxyCount).Sum(g => counts[g] * contrasts[g][j]) / pointCount,
                ["percentile_ci_low"] = q[0],
                ["percentile_ci_high"] = q[1],
                ["estimand"] = "pooled_point_weighted_ratio; clusters resampled and denominator recomputed",
                ["interpretation"] = "secondary_weighting_sensitivity_no_primary_p"
            });
        }
        Csv.Write(Path.Combine(destination, "pooled_point_sensitivity.csv"), pooledResults);

        Result profile = ProfileSensitivity(galaxies, counts, modelLosses);
        Data.SaveJson(Path.Combine(destination, "profiled_baseline_failure_sensitivity.json"), profile);
        var residuals = ResidualSummary(galaxyGroups);
        Csv.Write(Path.Combine(destination, "residual_distribution_summary.csv"), residuals);
        SaveBootstrap(Path.Combine(destination, "bootstrap_summaries.npz"), bootstrapMeans, bootstrapTStatistics, pooledBootstrapMeans);

        var checks = new Result
        {
            ["status"] = "pass",
            ["galaxies"] = galaxyCount,
            ["points"] = pointCount,
            ["draws"] = Draws,
            ["seed"] = Seed,
            ["bootstrap_zero_se"] = new int[scoreCount],
            ["wild_zero_se"] = new int[comparatorCount],
            ["input_sha256"] = hash,
            ["methods"] = "PROTOCOL.json",
            ["primary_tests"] = comparatorCount,
            ["finite_profiled_subset"] = "130galaxies/654points",
            ["no_refit"] = true,
            ["scope"] = "conditional retrospective inference assuming independent galaxy clusters; no history correction"
        };
        Data.SaveJson(Path.Combine(destination, "checks.json"), checks);
        Result validation = Validate(source, destination, profile, checks);
        validation["self_checks"] = SelfChecks();
        Data.SaveJson(Path.Combine(destination, "validation.json"), validation);
        var report = new Result
        {
            ["primary"] = primary,
            ["score_sensitivity"] = results,
            ["residual_distribution"] = residuals,
            ["pooled_point_sensitivity"] = pooledResults,
            ["profiled_baseline_failure"] = profile,
            ["validation"] = validation,
            ["monte_carlo_floor_note"] = "Three zero-exceedance comparisons reach the Monte Carlo resolution limit; p=0.00001 raw and p=0.00004 Holm4 (Z=4.10748) are finite-resolution reporting values, not precise extreme significance.",
            ["scope"] = checks["scope"]
        };
        Data.SaveJson(Path.Combine(destination, "results.json"), report);
        var core = results.Single(r => (string)r["score"]! == "absolute_kms" && (string)r["comparator"]! == "scalar_core");
        Console.WriteLine($"Pointwise: core MAE {residuals[1]["mean_abs_kms_equal_galaxy"]:F5} -> envelope {residuals[0]["mean_abs_kms_equal_galaxy"]:F5} km/s; reduction {-(double)core["effect"]!:F5}, 95% CI [{-(double)core["bootstrap_t_ci_high"]!:F5}, {-(double)core["bootstrap_t_ci_low"]!:F5}].");
        Console.WriteLine($"Pointwise: NFW Holm4 p={holmPValues[2]:F5}, Z={Statistics.NormalMagnitude(holmPValues[2]):F5}; {exceedances.Count(tailCount => tailCount == 0)} comparisons at Monte Carlo resolution limit.");
        return report;
    }

    static double[] ReadCorrelationSensitivities(string source)
    {
        var protocol = Data.Json(Path.Combine(source, "PROTOCOL.json"));
        Require(protocol.GetProperty("bootstrap_draws").GetInt32() == Draws && protocol.GetProperty("seed").GetUInt64() == Seed, "Changed frozen bootstrap protocol.");
        Require(protocol.GetProperty("primary_comparators").EnumerateArray().Select(x => x.GetString()).SequenceEqual(Comparators), "Changed primary comparison family.");
        double[] rhos = protocol.GetProperty("sensitivities").GetProperty("correlation_rho").EnumerateArray().Select(x => x.GetDouble()).ToArray();
        Require(rhos.SequenceEqual(new double[] { 0, .25, .5, .75 }), "Changed correlation sensitivities.");
        return rhos;
    }

    static (double[][] Loss, double[][] MeanResiduals, double[][] ResidualVariances) SummarizeGalaxyResiduals(IGrouping<string, Point>[] galaxyGroups)
    {
        int galaxyCount = galaxyGroups.Length;
        double[][] modelLosses = new double[galaxyCount][], meanResiduals = new double[galaxyCount][], residualVariances = new double[galaxyCount][];
        for (int g = 0; g < galaxyCount; g++)
        {
            modelLosses[g] = new double[Models.Length];
            meanResiduals[g] = new double[Models.Length];
            residualVariances[g] = new double[Models.Length];
            for (int m = 0; m < Models.Length; m++)
            {
                double[] z = galaxyGroups[g].Select(p => p.Z[m]).ToArray();
                meanResiduals[g][m] = z.Average();
                residualVariances[g][m] = z.Sum(v => Square(v - meanResiduals[g][m])) / z.Length;
                modelLosses[g][m] = galaxyGroups[g].Any(p => p.Failed[m]) ? double.PositiveInfinity : z.Sum(Square) / z.Length;
            }
        }
        return (modelLosses, meanResiduals, residualVariances);
    }

    static void WritePairedDifferences(string destination, IGrouping<string, Point>[] galaxyGroups)
    {
        var paired = new List<Result>();
        foreach (string comparator in Comparators.Append("profiled_baryons"))
        {
            int m = Array.IndexOf(Models, comparator);
            foreach (var group in galaxyGroups)
                foreach (Point point in group)
                {
                    double a = point.Z[0], b = point.Z[m];
                    paired.Add(new Result
                    {
                        ["galaxy"] = point.Galaxy,
                        ["point_index"] = point.Index,
                        ["radius_kpc"] = point.Radius,
                        ["radius_fraction_of_last"] = point.Radius / group.Max(p => p.Radius),
                        ["comparator"] = comparator,
                        ["z_envelope"] = a,
                        ["z_comparator"] = b,
                        ["delta_squared"] = a * a - b * b,
                        ["delta_absolute"] = Math.Abs(a) - Math.Abs(b),
                        ["delta_absolute_kms"] = (Math.Abs(a) - Math.Abs(b)) * point.Sigma,
                        ["comparison_valid"] = double.IsFinite(a) && double.IsFinite(b),
                        ["comparator_fit_failed"] = point.Failed[m]
                    });
                }
        }
        Csv.Write(Path.Combine(destination, "pointwise_paired_differences.csv"), paired);
    }

    static List<Descriptor> CreateScoreDescriptors(double[] rhos)
    {
        var descriptors = new List<Descriptor>();
        foreach (string score in new[] { "squared_standardized", "absolute_standardized", "absolute_kms" })
            descriptors.AddRange(Comparators.Select(c => new Descriptor(score, c)));
        foreach (double rho in rhos.Skip(1))
            descriptors.AddRange(Comparators.Select(c => new Descriptor("correlated_quadratic", c, rho)));
        return descriptors;
    }

    static double[][] CalculateGalaxyContrasts(IGrouping<string, Point>[] galaxyGroups, int[] counts,
        List<Descriptor> descriptors, double[][] modelLosses, double[][] meanResiduals, double[][] residualVariances)
    {
        int galaxyCount = galaxyGroups.Length, scoreCount = descriptors.Count;
        var contrasts = new double[galaxyCount][];
        for (int g = 0; g < galaxyCount; g++)
        {
            contrasts[g] = new double[scoreCount];
            for (int j = 0; j < scoreCount; j++)
            {
                var desc = descriptors[j];
                int m = Array.IndexOf(Models, desc.Comparator);
                contrasts[g][j] = desc.Score switch
                {
                    "squared_standardized" => modelLosses[g][0] - modelLosses[g][m],
                    "absolute_standardized" => galaxyGroups[g].Average(p => Math.Abs(p.Z[0]) - Math.Abs(p.Z[m])),
                    "absolute_kms" => galaxyGroups[g].Average(p => (Math.Abs(p.Z[0]) - Math.Abs(p.Z[m])) * p.Sigma),
                    _ => Correlated(residualVariances[g][0], meanResiduals[g][0], counts[g], desc.Rho!.Value) - Correlated(residualVariances[g][m], meanResiduals[g][m], counts[g], desc.Rho.Value)
                };
            }
        }
        Require(contrasts.SelectMany(x => x).All(double.IsFinite), "Non-finite predeclared continuous contrast.");
        return contrasts;
    }

    static (double[][] Means, double[][] TStatistics, double[][] PooledMeans) BootstrapPairedGalaxies(
        double[][] contrasts, int[] counts, double[] means)
    {
        int galaxyCount = contrasts.Length, scoreCount = means.Length, comparatorCount = Comparators.Length;
        // Every sampled galaxy contributes its full radial block, with equal primary weight.
        double[][] bootstrapMeans = Matrix(scoreCount, Draws), bootstrapTStatistics = Matrix(scoreCount, Draws), pooledBootstrapMeans = Matrix(comparatorCount, Draws);
        var random = new PointwiseRandom(Seed);
        int[] sampledGalaxies = new int[galaxyCount];
        double[] sampleMeans = new double[scoreCount], sampleSquaredDeviations = new double[scoreCount], sampleWeightedSums = new double[comparatorCount];
        for (int b = 0; b < Draws; b++)
        {
            Array.Clear(sampleMeans);
            Array.Clear(sampleSquaredDeviations);
            Array.Clear(sampleWeightedSums);
            int denominator = 0;
            for (int i = 0; i < galaxyCount; i++)
            {
                int g = sampledGalaxies[i] = random.NextInt(galaxyCount);
                double[] row = contrasts[g];
                denominator += counts[g];
                for (int j = 0; j < scoreCount; j++)
                    sampleMeans[j] += row[j];
                for (int j = 0; j < comparatorCount; j++)
                    sampleWeightedSums[j] += row[j] * counts[g];
            }
            for (int j = 0; j < scoreCount; j++)
                sampleMeans[j] /= galaxyCount;
            foreach (int g in sampledGalaxies)
                for (int j = 0; j < scoreCount; j++)
                    sampleSquaredDeviations[j] += Square(contrasts[g][j] - sampleMeans[j]);
            for (int j = 0; j < scoreCount; j++)
            {
                double se = Math.Sqrt(sampleSquaredDeviations[j] / (galaxyCount - 1)) / Math.Sqrt(galaxyCount);
                Require(se > 0 && double.IsFinite(se), "Undefined pairs bootstrap standard error.");
                bootstrapMeans[j][b] = sampleMeans[j];
                bootstrapTStatistics[j][b] = (sampleMeans[j] - means[j]) / se;
            }
            // Pooled weighting recomputes its denominator after each cluster resample.
            for (int j = 0; j < comparatorCount; j++)
                pooledBootstrapMeans[j][b] = sampleWeightedSums[j] / denominator;
        }
        return (bootstrapMeans, bootstrapTStatistics, pooledBootstrapMeans);
    }
    static Point[] ReadPoints(string path)
    {
        var rows = Csv.Read(path);
        Require(rows.Count == 659 * Models.Length, "Wrong frozen prediction row count.");
        return rows.GroupBy(r => (Galaxy: r["galaxy"], Index: int.Parse(r["point_index"], CultureInfo.InvariantCulture)))
            .OrderBy(g => g.Key.Galaxy, StringComparer.Ordinal).ThenBy(g => g.Key.Index).Select(group =>
            {
                Require(group.Count() == Models.Length && group.Select(r => r["model"]).ToHashSet().SetEquals(Models), "Missing or duplicate paired model prediction.");
                var first = group.First();
                double radius = Csv.Number(first, "radius_kpc"), observed = Csv.Number(first, "observed_kms"), sigma = Csv.Number(first, "sigma_kms");
                Require(double.IsFinite(radius) && radius > 0 && double.IsFinite(observed) && double.IsFinite(sigma) && sigma > 0, "Invalid observation or uncertainty.");
                var z = new double[Models.Length];
                var failed = new bool[Models.Length];
                foreach (var row in group)
                {
                    Require(Csv.Number(row, "radius_kpc") == radius && Csv.Number(row, "observed_kms") == observed && Csv.Number(row, "sigma_kms") == sigma, "Unmatched paired observations.");
                    Require(Csv.Flag(row, "is_outer") && !Csv.Flag(row, "is_training"), "Training row in outer analysis.");
                    int m = Array.IndexOf(Models, row["model"]);
                    z[m] = (Csv.Number(row, "predicted_kms") - observed) / sigma;
                    failed[m] = Csv.Flag(row, "fit_failed");
                    Require(m == 5 || double.IsFinite(z[m]) && !failed[m], "Failed primary or provenance model.");
                    Require(Csv.Flag(row, "valid_prediction") == double.IsFinite(z[m]), "Invalid prediction flag.");
                }
                return new Point(group.Key.Galaxy, group.Key.Index, radius, observed, sigma, z, failed);
            }).ToArray();
    }

    static int[] WildBootstrap(double[][] contrasts, double[] means, double[] errors)
    {
        int galaxyCount = contrasts.Length, comparatorCount = Comparators.Length;
        var random = new PointwiseRandom(Seed + 1);
        var exceedances = new int[comparatorCount];
        var signedSums = new double[comparatorCount];
        double[] squaredSums = Enumerable.Range(0, comparatorCount).Select(j => contrasts.Sum(row => row[j] * row[j])).ToArray();
        for (int b = 0; b < Draws; b++)
        {
            // NumPy starts a fresh byte buffer for each 5,000-draw call.
            if (b % 5000 == 0)
                random.StartInt8Batch();
            Array.Clear(signedSums);
            for (int g = 0; g < galaxyCount; g++)
            {
                int sign = 2 * random.NextInt8(2) - 1;
                for (int j = 0; j < comparatorCount; j++)
                    signedSums[j] += sign * contrasts[g][j];
            }
            for (int j = 0; j < comparatorCount; j++)
            {
                double sampleMean = signedSums[j] / galaxyCount, variance = (squaredSums[j] - galaxyCount * sampleMean * sampleMean) / (galaxyCount - 1), sampleError = Math.Sqrt(Math.Max(variance, 0) / galaxyCount);
                Require(sampleError > 0 && double.IsFinite(sampleError), "Undefined wild bootstrap standard error.");
                if (Math.Abs(sampleMean / sampleError) >= Math.Abs(means[j] / errors[j]))
                    exceedances[j]++;
            }
        }
        return exceedances;
    }

    // The failed profiled fit is retained; only its explicitly conditional subset is bootstrapped.
    static Result ProfileSensitivity(string[] galaxies, int[] counts, double[][] modelLosses)
    {
        int[] feasibleGalaxies = Enumerable.Range(0, galaxies.Length).Where(g => double.IsFinite(modelLosses[g][5])).ToArray();
        Require(feasibleGalaxies.Length == 130 && feasibleGalaxies.Sum(g => counts[g]) == 654, "Incorrect feasible profiled subset.");
        double[] finiteDifferences = feasibleGalaxies.Select(g => modelLosses[g][0] - modelLosses[g][5]).ToArray();
        var (mean, se) = MeanSe(finiteDifferences);
        var random = new PointwiseRandom(Seed + 2);
        var profileMeans = new double[Draws];
        var profileTStatistics = new double[Draws];
        var sample = new double[finiteDifferences.Length];
        for (int b = 0; b < Draws; b++)
        {
            for (int i = 0; i < finiteDifferences.Length; i++)
                sample[i] = finiteDifferences[random.NextInt(finiteDifferences.Length)];
            var (sampleMean, sampleError) = MeanSe(sample);
            Require(sampleError > 0 && double.IsFinite(sampleError), "Undefined finite-profile bootstrap standard error.");
            profileMeans[b] = sampleMean;
            profileTStatistics[b] = (sampleMean - mean) / sampleError;
        }
        var tQuantiles = QuantilePair(profileTStatistics);
        var meanQuantiles = QuantilePair(profileMeans);
        return new Result
        {
            ["model"] = "profiled_baryons",
            ["full_sample_galaxies"] = galaxies.Length,
            ["full_sample_points"] = counts.Sum(),
            ["failed_galaxies"] = Enumerable.Range(0, galaxies.Length).Where(g => !double.IsFinite(modelLosses[g][5])).Select(g => galaxies[g]).ToArray(),
            ["full_sample_effect"] = "minus_infinity_due_to_comparator_failure; finite mean inference undefined",
            ["full_sample_p"] = null,
            ["finite_subset_galaxies"] = feasibleGalaxies.Length,
            ["finite_subset_points"] = feasibleGalaxies.Sum(g => counts[g]),
            ["conditional_finite_effect"] = mean,
            ["conditional_finite_cluster_se"] = se,
            ["conditional_finite_ci_low"] = mean - tQuantiles[1] * se,
            ["conditional_finite_ci_high"] = mean - tQuantiles[0] * se,
            ["conditional_finite_percentile_low"] = meanQuantiles[0],
            ["conditional_finite_percentile_high"] = meanQuantiles[1],
            ["primary_p_assigned"] = false,
            ["scope"] = "conditional-on-feasible130 sensitivity, not fullsample claim"
        };
    }

    static List<Result> ResidualSummary(IGrouping<string, Point>[] galaxyGroups)
    {
        var result = new List<Result>();
        for (int m = 0; m < Models.Length; m++)
        {
            double total = 0, squared = 0, abs = 0, kms = 0, gt3 = 0, gt5 = 0;
            int finite = 0, failed = 0;
            foreach (var group in galaxyGroups)
                foreach (Point point in group)
                {
                    double z = point.Z[m], weight = 1.0 / group.Count();
                    if (!double.IsFinite(z))
                    {
                        failed++;
                        continue;
                    }
                    finite++;
                    total += weight;
                    squared += weight * z * z;
                    abs += weight * Math.Abs(z);
                    kms += weight * Math.Abs(z) * point.Sigma;
                    if (Math.Abs(z) > 3)
                        gt3 += weight;
                    if (Math.Abs(z) > 5)
                        gt5 += weight;
                }
            result.Add(new Result
            {
                ["model"] = Models[m],
                ["finite_points"] = finite,
                ["failed_points"] = failed,
                ["mean_squared_equal_galaxy"] = squared / total,
                ["mean_abs_standardized_equal_galaxy"] = abs / total,
                ["mean_abs_kms_equal_galaxy"] = kms / total,
                ["fraction_abs_z_gt3"] = gt3 / total,
                ["fraction_abs_z_gt5"] = gt5 / total,
                ["finite_conditioning"] = failed > 0
            });
        }
        return result;
    }

    static void SaveBootstrap(string path, double[][] means, double[][] tstats, double[][] pooled)
    {
        using var file = File.Create(path);
        using var zip = new ZipArchive(file, ZipArchiveMode.Create);
        foreach (var (name, matrix) in new[] { ("means", means), ("tstats", tstats), ("pooled_means", pooled) })
        {
            using var stream = zip.CreateEntry(name + ".npy", CompressionLevel.Optimal).Open();
            using var writer = new BinaryWriter(stream, Encoding.ASCII, true);
            string header = $"{{'descr': '<f8', 'fortran_order': False, 'shape': ({Draws}, {matrix.Length}), }}";
            header += new string(' ', (64 - (10 + header.Length + 1) % 64) % 64) + "\n";
            writer.Write(new byte[] { 0x93, (byte)'N', (byte)'U', (byte)'M', (byte)'P', (byte)'Y', 1, 0 });
            writer.Write((ushort)header.Length);
            writer.Write(Encoding.ASCII.GetBytes(header));
            for (int b = 0; b < Draws; b++)
                for (int j = 0; j < matrix.Length; j++)
                    writer.Write(matrix[j][b]);
        }
    }

    static double[][] Matrix(int columns, int rows) => Enumerable.Range(0, columns).Select(_ => new double[rows]).ToArray();
    static double Square(double x) => x * x;
    static double Correlated(double variance, double mean, int n, double rho) => variance / (1 - rho) + mean * mean / (1 + (n - 1) * rho);
    static (double Mean, double Se) MeanSe(double[] d)
    {
        double mean = d.Average();
        return (mean, Math.Sqrt(d.Sum(x => Square(x - mean)) / (d.Length - 1)) / Math.Sqrt(d.Length));
    }
    static double[] QuantilePair(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        // R7 matches the protocol's NumPy linear quantiles; these draws are finite.
        return [SortedArrayStatistics.QuantileCustom(sorted, .025, QuantileDefinition.R7),
            SortedArrayStatistics.QuantileCustom(sorted, .975, QuantileDefinition.R7)];
    }
    static void Require(bool condition, string message)
    {
        if (!condition)
            throw new InvalidDataException("Pointwise analysis: " + message);
    }
}
