using System.Diagnostics;
using Result = System.Collections.Generic.Dictionary<string, object?>;

namespace DarkUniverse;

/// <summary>
/// Synthetic finite-sample calibration, separate from manuscript reproduction.
/// Shared-systematic scenarios intentionally violate the engine's independence assumption.
/// </summary>
public static class CalibrationStudy
{
    sealed record Scenario(string Name, string Description, bool Independent);
    static readonly Scenario[] Scenarios =
    [
        new("normal", "Independent standard normal galaxy effects.", true),
        new("skewed", "Independent centered exponential galaxy effects, mean zero and variance one.", true),
        new("heteroskedastic", "Independent normal effects with a smooth 0.2-to-3 standard-deviation gradient, normalized to average variance one.", true),
        new("influential", "Independent normal effects with two high-variance galaxies (scale 8) and remaining scale 0.3, normalized to average variance one.", true),
        new("shared_systematic", "A common N(0,0.25) offset plus independent N(0,0.75) terms; pair correlation 0.25 intentionally violates galaxy independence.", false)
    ];

    public static Result Run(string dataRoot, string outputRoot, AnalysisStudyProtocol? protocol = null)
    {
        _ = dataRoot; // Synthetic experiments do not read or refit observed galaxies.
        protocol ??= new(); protocol.Validate();
        string output = Path.Combine(outputRoot, "research", "calibration");
        Directory.CreateDirectory(output);
        var stopwatch = Stopwatch.StartNew();
        int R = protocol.SimulationReplicates, G = protocol.GalaxyCount;
        var results = new Result[checked(Scenarios.Length * R * 2)];
        Console.WriteLine($"Calibration: {Scenarios.Length} synthetic scenarios, null and known effect, {R} replicates each, {protocol.BootstrapDraws} draws.");
        for (int s = 0; s < Scenarios.Length; s++)
        {
            int scenarioIndex = s;
            var scenario = Scenarios[s];
            Parallel.For(0, R, new ParallelOptions { MaxDegreeOfParallelism = protocol.MaxParallelism }, replicate =>
            {
                ulong stream = checked((ulong)(scenarioIndex + 1) * 1000000000UL + (ulong)replicate * 3);
                ulong dataSeed = StatisticalEngine.DeriveSeed(protocol.Seed, stream);
                double[] errors = Generate(scenarioIndex, G, dataSeed);
                var options = new PairedBootstrapOptions
                {
                    Draws = protocol.BootstrapDraws, ConfidenceLevel = protocol.ConfidenceLevel,
                    PairsSeed = StatisticalEngine.DeriveSeed(protocol.Seed, stream + 1),
                    WildSeed = StatisticalEngine.DeriveSeed(protocol.Seed, stream + 2)
                };
                for (int effectIndex = 0; effectIndex < 2; effectIndex++)
                {
                    double truth = effectIndex == 0 ? 0 : protocol.StandardizedEffect;
                    var input = errors.Select((e, g) => new PairedGalaxyObservation($"g{g:D6}", 2 + g % 8, e + truth)).ToArray();
                    var run = StatisticalEngine.Run(input, options);
                    bool? covered = run.BootstrapTInterval is { } interval ? interval.Lower <= truth && truth <= interval.Upper : null;
                    bool? rejected = run.WildPValue is double p ? p <= protocol.TestAlpha : null;
                    results[(scenarioIndex * R + replicate) * 2 + effectIndex] = new()
                    {
                        ["scenario"] = scenario.Name, ["case"] = effectIndex == 0 ? "null" : "known_effect",
                        ["replicate"] = replicate, ["true_effect"] = truth, ["estimated_effect"] = run.Effect,
                        ["error"] = run.Effect - truth, ["cluster_se"] = run.ClusterStandardError,
                        ["ci_lower"] = run.BootstrapTInterval?.Lower, ["ci_upper"] = run.BootstrapTInterval?.Upper,
                        ["covers_true_effect"] = covered, ["wild_p"] = run.WildPValue, ["rejects_zero"] = rejected,
                        ["wild_exceedances"] = run.WildExceedances, ["mc_floor"] = run.MonteCarloFloor,
                        ["maximum_variance_share"] = run.MaximumVarianceShare, ["variance_concentration_count"] = run.VarianceConcentrationCount,
                        ["status"] = run.Status, ["degenerate_pairs_draws"] = run.DegeneratePairsDraws,
                        ["degenerate_wild_draws"] = run.DegenerateWildDraws,
                        ["data_seed"] = dataSeed, ["pairs_seed"] = options.PairsSeed, ["wild_seed"] = options.WildSeed
                    };
                }
            });
            Console.WriteLine($"Calibration {scenario.Name}: completed {R} null and {R} known-effect replicates.");
        }
        var summary = new List<Result>();
        foreach (var scenario in Scenarios)
        foreach (string effectCase in new[] { "null", "known_effect" })
        {
            var rows = results.Where(r => (string)r["scenario"]! == scenario.Name && (string)r["case"]! == effectCase).ToArray();
            int validTests = rows.Count(r => r["rejects_zero"] is bool), rejects = rows.Count(r => r["rejects_zero"] is true);
            int validIntervals = rows.Count(r => r["covers_true_effect"] is bool), covers = rows.Count(r => r["covers_true_effect"] is true);
            var rejectionMc = validTests > 0 ? StatisticalEngine.BinomialInterval(rejects, validTests) : null;
            var coverageMc = validIntervals > 0 ? StatisticalEngine.BinomialInterval(covers, validIntervals) : null;
            summary.Add(new()
            {
                ["scenario"] = scenario.Name, ["case"] = effectCase, ["description"] = scenario.Description,
                ["independence_assumption_holds"] = scenario.Independent, ["attempted_replicates"] = R,
                ["valid_tests"] = validTests, ["invalid_tests"] = R - validTests, ["rejects_zero"] = rejects,
                ["rate_interpretation"] = effectCase == "null" ? "type_I_rejection" : "power_for_configured_nonzero_mean",
                ["rejection_rate_valid_only"] = validTests > 0 ? rejects / (double)validTests : null,
                ["rejection_mc95_lower"] = rejectionMc?.Lower, ["rejection_mc95_upper"] = rejectionMc?.Upper,
                ["nominal_test_alpha"] = protocol.TestAlpha,
                ["null_rejection_mc_interval_contains_nominal"] = effectCase == "null" && rejectionMc is not null ? rejectionMc.Lower <= protocol.TestAlpha && protocol.TestAlpha <= rejectionMc.Upper : null,
                ["valid_intervals"] = validIntervals, ["invalid_intervals"] = R - validIntervals, ["covers_true_effect"] = covers,
                ["coverage_valid_only"] = validIntervals > 0 ? covers / (double)validIntervals : null,
                ["coverage_mc95_lower"] = coverageMc?.Lower, ["coverage_mc95_upper"] = coverageMc?.Upper,
                ["nominal_interval_coverage"] = protocol.ConfidenceLevel,
                ["coverage_mc_interval_contains_nominal"] = coverageMc is not null ? coverageMc.Lower <= protocol.ConfidenceLevel && protocol.ConfidenceLevel <= coverageMc.Upper : null,
                ["mean_bias"] = rows.Average(r => (double)r["error"]!),
                ["rmse"] = Math.Sqrt(rows.Average(r => StatisticalEngine.Square((double)r["error"]!))),
                ["mean_reported_cluster_se"] = rows.Average(r => (double)r["cluster_se"]!),
                ["mean_ci_width_valid_only"] = validIntervals > 0 ? rows.Where(r => r["covers_true_effect"] is bool).Average(r => (double)r["ci_upper"]! - (double)r["ci_lower"]!) : null,
                ["mean_maximum_variance_share"] = rows.Average(r => (double)r["maximum_variance_share"]!),
                ["mc_floor_replicates"] = rows.Count(r => r["mc_floor"] is true)
            });
        }
        Csv.Write(Path.Combine(output, "replicates.csv"), results);
        Csv.Write(Path.Combine(output, "scenario_summary.csv"), summary);
        var report = new Result
        {
            ["schema_version"] = 1, ["status"] = "completed", ["protocol"] = protocol,
            ["random_algorithm"] = Pcg64Random.ResearchAlgorithm, ["elapsed_seconds"] = stopwatch.Elapsed.TotalSeconds,
            ["scenario_count"] = Scenarios.Length, ["total_simulated_datasets"] = results.Length,
            ["purpose"] = "Measure finite-sample null rejection, known-effect power and interval coverage; successful execution does not mean a method is calibrated.",
            ["scope"] = "Synthetic stress tests, not validation of actual SPARC galaxies or a fitted physical model. Shared-systematic scenarios intentionally violate independent-galaxy assumptions.",
            ["paired_scenarios"] = "Null and known-effect datasets reuse the same error realization and bootstrap streams; comparisons across these cases are correlated by design.",
            ["monte_carlo_intervals"] = "Marginal 95% Clopper-Pearson intervals for simulation rejection and coverage proportions; not simultaneous across scenarios.",
            ["invalid_policy"] = "Invalid tests and intervals are counted separately. Reported rates condition explicitly on valid outcomes; degenerate draws are never silently discarded within a run.",
            ["bootstrap_resolution_floor"] = 1.0 / (protocol.BootstrapDraws + 1),
            ["alpha_below_resolution_floor"] = protocol.TestAlpha < 1.0 / (protocol.BootstrapDraws + 1),
            ["summary"] = summary
        };
        Data.SaveJson(Path.Combine(output, "report.json"), report);
        return report;
    }

    static double[] Generate(int scenario, int galaxies, ulong seed)
    {
        var rng = Pcg64Random.FromSeed(seed);
        var result = new double[galaxies];
        double[] scales = Enumerable.Range(0, galaxies).Select(i => scenario switch
        {
            2 => .2 + 2.8 * i / (galaxies - 1),
            3 => i < 2 ? 8 : .3,
            _ => 1.0
        }).ToArray();
        double rms = Math.Sqrt(scales.Average(StatisticalEngine.Square));
        double shared = scenario == 4 ? .5 * rng.NextNormal() : 0;
        for (int g = 0; g < galaxies; g++)
            result[g] = scenario switch
            {
                1 => -Math.Log(rng.NextOpenDouble()) - 1,
                2 or 3 => rng.NextNormal() * scales[g] / rms,
                4 => shared + Math.Sqrt(.75) * rng.NextNormal(),
                _ => rng.NextNormal()
            };
        return result;
    }
}

