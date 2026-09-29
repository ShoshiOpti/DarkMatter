# Statistical constraints for research studies

These studies extend the analysis without replacing the frozen manuscript protocols, seeds, or results. The example JSON files deserialize to AnalysisStudyProtocol; omitted fields retain its defaults. This type is separate from the locked fresh-fitting ResearchWorkflow protocol.

From the repository root:

~~~powershell
dotnet run --project DarkMatter -c Release -- calibrate --protocol DarkMatter/research/calibration.example.json
dotnet run --project DarkMatter -c Release -- resolution --protocol DarkMatter/research/resolution.example.json
~~~

The application records each run separately. Within its run directory, calibration writes research/calibration/report.json, scenario_summary.csv and replicates.csv. Resolution writes research/mechanism_resolution/report.json, galaxy_contrasts.csv and pointwise_differences.csv.

## Reusable paired-galaxy API

~~~csharp
var options = new PairedBootstrapOptions
{
    Draws = 999,
    PairsSeed = 123,
    WildSeed = 456,
    ConfidenceLevel = 0.95
};
PairedBootstrapResult result = StatisticalEngine.Run(
    new PairedGalaxyObservation[]
    {
        // Supply one finite, paired mean-loss difference per galaxy.
        // This abbreviated example describes the type, not an adequate study.
        new("galaxy_A", PointCount: 5, Difference: -0.4),
        new("galaxy_B", PointCount: 8, Difference: 0.2)
    }, options);
~~~

Each observation contains a unique galaxy identifier, its number of evaluation points, and the mean paired loss difference within that galaxy. The primary effect gives galaxies equal weight. The separate pooled-point effect weights by point count and recomputes the denominator after every whole-galaxy resample. Sorting identifiers makes the seeded result invariant to input row order.

The engine returns one **unadjusted** comparison. A caller must retain its declared comparison family and apply its chosen multiplicity adjustment. It must also define the loss, evaluation split, eligible population, direction of the contrast and failure policy before evaluation. A minimum galaxy-count rule is an operational gate; no fixed count establishes adequate calibration.

Research streams use versioned PCG64 with SplitMix64 seed expansion. Their seed mapping is intentionally distinct from NumPy SeedSequence. The frozen manuscript adapters retain their original NumPy states, buffering rules and full-stream checksums. Research options cannot silently change a published run.

## Failures and degenerate draws

The reusable engine rejects duplicate identifiers, nonpositive point counts and nonfinite differences. The caller must handle a failed fit at the galaxy/model level; it must not hide a training failure behind finite outer predictions.

For the manuscript convention, finite models beat failed models in signs, joint failures have no sign, and an undefined full-sample continuous contrast receives no primary p value. A feasible-subset calculation is a separately labeled conditional sensitivity. A p=1 placeholder may preserve a fixed Holm family internally, but it is not the failed comparison's published p value.

The engine handles degeneracy explicitly:

- All differences exactly zero produce an exact zero interval and p=1.
- A nonzero constant effect has undefined studentized inference.
- Any zero-variance pairs draw invalidates the bootstrap-t interval.
- Any zero-variance wild draw invalidates the wild p value.
- No degenerate draws are silently dropped. Counts and status remain in the result. Descriptive percentile intervals, or a valid wild test when only pairs draws failed, may still be available.

Callers must inspect the nullable inferential fields and status; successful program execution does not guarantee an interpretable statistical test.

## Intervals, tests and Monte Carlo precision

The engine computes a marginal whole-galaxy pairs-bootstrap-t interval and a separate null-imposed Rademacher wild test. These procedures are not inverses. An interval can exclude zero while the wild p value exceeds the selected threshold. Neither the unadjusted interval nor the test becomes simultaneous across a family merely because other p values receive Holm adjustment.

The wild test is approximate under a mean-only null. Exact sign-randomization reasoning requires stronger symmetry assumptions. Independent galaxy clusters and adequate contribution from multiple clusters matter; the reported variance-concentration count is a diagnostic, not a literal effective sample size or a replacement degrees of freedom.

The wild p value is (exceedances + 1)/(draws + 1). Zero exceedances are flagged at the simulation resolution floor. A Clopper-Pearson interval describes uncertainty in the Monte Carlo tail count, not uncertainty from fitting, galaxy selection or shared systematics. More draws refine simulation precision; they do not correct an invalid sampling assumption.

## What calibration establishes

The calibration example uses 500 independent synthetic realizations per case and 999 bootstrap draws per realization. Every scenario has a null case and a known nonzero mean. The scenarios cover normal, skewed, smoothly heteroskedastic, strongly influential and shared-systematic effects.

Reports contain null rejection rates, known-effect power and coverage of the true mean. Each rate has a marginal 95% binomial Monte Carlo interval. Invalid tests and intervals are counted, and reported rates explicitly condition on valid outcomes. Comparisons across scenarios are not simultaneous tests.

Null and shifted cases share their simulated error realization and random streams. Their paired results are correlated by design; they must not be pooled as independent simulations. The shared-systematic scenario deliberately violates independent-galaxy assumptions. Poor performance there demonstrates a limitation of the method under that constructed dependence; it does not estimate the dependence present in SPARC.

The simulations do not certify the real fitted models. Frozen-prediction inference does not propagate fitting uncertainty, common microscopic parameters, inherited calibration choices or model-development history. Those require a separately defined refitting or shared-uncertainty analysis.

## Choosing mechanism margins and numerical budgets

The resolution example leaves all three margin/budget settings null. Its default run measures differences and uncertainty without asserting equivalence. A nonsignificant p value is never interpreted as equivalence.

Before setting these fields, record the scientific justification in the study protocol:

- **mechanismSpeedToleranceKms** is an absolute tolerance for the largest charged/control speed difference among the retained outer observations. It is not a bound over unobserved radii.
- **mechanismLossMargin** is a symmetric practical margin for the equal-galaxy mean squared-standardized-loss difference. Its units differ from the speed tolerance.
- **solverSpeedErrorBoundKms** is a supplied uniform absolute error bound on each model's predicted speed at every evaluated point. Establish it from suitable numerical-convergence/error evidence for both calculations. The application neither estimates nor certifies this bound. Zero asserts no numerical error and requires justification.

The budget is propagated adversarially: the speed-difference bound adds twice the supplied speed error, and the loss bound includes both residual-linear and squared-error terms before equal-galaxy averaging. No cancellation or independence of solver errors is assumed.

A speed assessment requires both its margin and the solver bound. A loss assessment likewise requires its margin and the solver bound, and expands a pairs-bootstrap-t interval by the propagated numerical error. The inclusion interval uses confidence 1 - 2*testAlpha (90% at alpha=0.05), following two-one-sided-test-style interval inclusion. It is an approximate, conditional diagnostic rather than an inversion of the wild test. It remains subject to bootstrap calibration and galaxy-independence assumptions.

A configured margin is not a detected effect, a physical equivalence proof, or evidence for a formation mechanism. Preserve an unconfigured run and the scientific rationale for any subsequent margin; do not choose a margin just to obtain a preferred conclusion.

## Reconstruction update: unchanged inference and new diagnostics

The 27 September reconstruction update is not a new fitted model family. At fixed
total leading density, baryons, source geometry, boundary data and observation
mapping, the Poisson rotation curve is unchanged. Thus the same predictions,
residuals, galaxy weights, failure flags and original resampling protocol must
produce the same inferential results. The publication result contract checks this
numerically before showing the 25 baseline scenes. `statistical_impact.json` records
the current model summaries and that contract; it does not treat an identity as
an independent test or assign it a new p-value.

The density and force audits are descriptive and condition on the fitted target.
Local density fractions, enclosed-force fractions, and fractions of galaxies have
different denominators. Reported galaxy medians first take each galaxy's median
across its declared outer rows, then the median across the stated all/active sample.
The active subset uses the archived envelope threshold. Counting force decreases
relative to the original core diagnoses its amplitude change; it is not a rejection
of a newly fitted physical mechanism. A hypothetical double-counted envelope is
only a bookkeeping diagnostic, not an additional inferential competitor.

The periodic field fixture tests equations at a single instant and supplies no
sampling distribution or galaxy calibration. Its internal-stress changes cannot
be inserted into SPARC velocities without a physical normalization and evolved
state. Future forward predictions must use the locked research protocol, preserve
failure policy, and separate training from independent observational validation.
