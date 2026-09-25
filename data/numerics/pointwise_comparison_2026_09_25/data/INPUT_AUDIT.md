# Point-level prediction input audit

The new table reuses saved, frozen predictions. No model is refitted, no observed speed or error is transformed, and no missing prediction is replaced by a finite penalty.

## Audited sample and schema

- 131 galaxies; 3,034 distinct measured radii; 659 original outer radii.
- The outer block is the last max(2, ceil(0.2*n_g)) sorted radii in each galaxy; it contains 2–23 points per galaxy.
- Seven models give 4,613 outer model–point rows. All seven have exactly the same galaxy/radius/observation/error/split indexing.
- The canonical coordinate metadata come from the saved matched scalar-core rows. Independent profiled-baryon data match the canonical observations, errors and radii to CSV roundoff, with maximum absolute discrepancy below 1.5e-14. The integer point_index is the original zero-based radius-row index within a galaxy.
- pointwise_outer.csv is the principal long table. outer_wide.csv has one row per outer measurement and each model's predicted speed as its model-named column.
- pointwise_all.csv additionally retains every original training row for audit, not for pooling with held-out evaluation.
- residual_kms = predicted_kms - observed_kms; standardized_residual = residual_kms/sigma_kms; squared_standardized_error is its square.
- valid_prediction, fit_failed, and failure_reason must be checked before finite-loss calculations. Missing values remain missing.
- galaxy_availability.csv, model_inventory.csv, reproduced_galaxy_losses.csv, and checks.json report all availability and agreement checks.

## Model meanings and conditioning

| ID | Saved model | Conditioning |
|---|---|---|
| extended | Fixed finite envelope, a/Rd=1 and cutoff ratio q=20, added to scalar core | Inner-fitted amplitudes; inherited matched-core nuisance parameters and core length held fixed |
| scalar_core | Matched scalar core | Original inner-only matched nuisance fit |
| compact_plummer | Compact Plummer control | Original inner-fitted control, sharing the specified conditional calibration |
| nfw | Matched NFW | Its own original inner-only matched nuisance fit |
| baryon_only | GR baryons with the scalar dark component removed | Conditional ablation of the inner-fitted core calibration; not a separate baryon-only optimization |
| profiled_baryons | Independently profiled GR baryons | Saved inner-only fit with matched stellar and geometry priors/bounds; no halo warm start |
| adaptive_primary | Conditional inner-selected envelope/core model | Its selector uses an inner validation subset, but inherited nuisance parameters used all original inner rows |

All speeds and quoted errors are in fixed catalogue coordinates. Distance and inclination transformations have already been applied to each model prediction; their fitted factors are preserved as metadata. Do not rescale the observations or errors again.

The independently profiled stellar priors have width 0.10 dex, with catalogue distance and inclination pulls, declared six-pull limits, and physical geometry bounds. This is a saved finite-multistart penalized fit, not a marginal posterior or a proven global optimum.

The original outer velocities were withheld from each model's direct fitting objective and the adaptive selector. However, SPARC and these outer outcomes have been used repeatedly in development. Frozen predictions do not turn this into a newly independent confirmatory sample. The adaptive comparison is additionally conditional on reused inner nuisance estimation. These limitations remain even when a more informative loss statistic replaces a win-rate statistic.

## Missing predictions and baseline feasibility

The independently profiled baryon pipeline fails for UGC01281 before nonlinear fitting: its signed baryon force is infeasible on the training rows within the declared stellar-prior bounds. Its five outer rows remain present with missing predictions, valid_prediction=False, fit_failed=True, and failure_reason=training_infeasible_within_prior_bounds.

There is no mathematically defined finite continuous loss difference for that galaxy against this comparator. A finite-loss sensitivity must be explicitly labeled 130 feasible galaxies / 654 outer points and report the retained feasibility failure separately. It must not be presented as the all-131-galaxy result, and the missing prediction must not be given an arbitrary large finite loss to manufacture significance.

All six other models have finite predictions at all 659 outer radii. The core-calibrated baryon ablation has invalid predictions at some inner radii; these are preserved in the all-row audit and do not affect the finite outer table.

## Covariance and uncertainty limits

The supplied SPARC fitting inputs contain individual tabulated speed errors, not a full within-galaxy radial covariance matrix. This is explicitly documented by the original uncertainty_sensitivity.py metadata (line 238 in the archived source), and the fitting residual uses elementwise division by those errors. The frozen source snapshot is inputs/nuisance_fitting_source.py.

Consequently, squaring and summing individual standardized residuals defines a useful measurement-error-weighted loss, but does not establish a chi-square null distribution with 659 independent degrees of freedom. Shared distance, inclination, stellar calibration, beam/radial systematics, and inherited fitted parameters can correlate points. New inference should preserve the galaxy as the resampling or clustering unit while allowing every radial point to contribute to a galaxy's continuous score.

The per-point sigma_kms values are measurement errors. They are not posterior predictive standard deviations and do not include uncertainty of the trained parameters. A frozen-prediction population bootstrap does not propagate a repeated full fitting-and-selection pipeline or unknown common systematic errors.

## Provenance and reproduction

input_manifest.json gives the original workspace path, SHA256 and byte count for every immutable source snapshot in inputs/. assemble_pointwise.py reads only those snapshots after the initial collection. It does not call a fitting routine. Run:

    python assemble_pointwise.py

All seven saved per-galaxy outer losses reproduce from the point table. The maximum absolute difference is 9.1e-13; model-specific differences are in checks.json. Every observation and split cross-match is checked, every output CSV is hashed in output_manifest.json, and source snapshots are rehashed after assembly. The outer-conditioned inverse curve is deliberately absent because it uses the very outer observations it reconstructs and is not a comparable frozen prediction.
