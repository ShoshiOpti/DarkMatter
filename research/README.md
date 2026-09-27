# Research workflows separate from the publication baseline

The native research commands create and lock study inputs, enforce prediction and fitting contracts, evaluate a predeclared comparison family, and run a small empirical refitting pilot. They do not revise the paper, the frozen publication statistics, or the retained original Python sources. [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md) documents the reusable inference engine, synthetic calibration and mechanism resolution commands; their example JSON protocols are a separate type from the fitting protocol described here.

## Run a development pilot

Run from the repository root. Use a Python environment with NumPy, SciPy, pandas and matplotlib (the retained scientific requirements are in `reference/requirements.txt`). `--packages` is optional and supplies a directory containing installed packages; no machine-specific package path is stored in the study protocol.

```powershell
dotnet run --project DarkMatter -c Release -- research-create --out output/research-drafts
```

The command prints its new run directory and records it in `output/research-drafts/latest.json`. Copy **only** `protocol.json`, `predictors.csv`, `training.csv`, `evaluation.csv`, `metadata.json` and `split.json` from that run into a new editable folder, for example `output/studies/empirical-001`. This keeps the original run artifact intact. Review and edit the draft and inputs there before locking. The default sample contains 131 eligible SPARC galaxies; the pilot selects DDO154 and NGC0300 in advance and runs a baseline and two perturbations.

```powershell
dotnet run --project DarkMatter -c Release -- research-lock --protocol output/studies/empirical-001/protocol.json --out output/research-runs
dotnet run --project DarkMatter -c Release -- research-validate --protocol output/studies/empirical-001/protocol.json --out output/research-runs
dotnet run --project DarkMatter -c Release -- research-pilot --protocol output/studies/empirical-001/protocol.json --reference reference --python python --out output/research-runs
```

Every command creates a separate run. Pilot output contains `predictions.csv`, `fit_artifact.json`, `fit_records.json`, `perturbations.json`, `adapter_provenance.json`, stdout/stderr logs, and copied original fit modules plus each replicate's envelope inputs/results. The native evaluator then writes `research_galaxy_scores.csv` and `research_evaluation.json`. The protocol's `pilot.timeout_seconds` is 3600 by default; timeout or Ctrl+C terminates the child process and retains logs. An incomplete run fails rather than publishing partial results.

To evaluate an external fitter's result, supply its actual prediction and artifact paths:

```powershell
dotnet run --project DarkMatter -c Release -- research-evaluate --protocol output/studies/empirical-001/protocol.json --predictions PATH_TO_PREDICTIONS.csv --fit-artifact PATH_TO_FIT_ARTIFACT.json --out output/research-runs
```

The numerical pilot was exercised with Python 3.12, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.1 and matplotlib 3.11.2. Baseline predictions for its two galaxies reproduced the corresponding archived predictions to within 3.1e-6 km/s; all 24 named fits across the baseline and two perturbations succeeded. That is an adapter regression check, not evidence about statistical power or independent predictive accuracy.

## What the empirical pilot actually refits

The adapter copies seven retained modules into an isolated output tree, records original and adapted hashes, and runs their original optimizers. It solves the dimensionless core profile, refits core and NFW halo parameters and galaxy nuisance parameters, invokes the depth-guard refit when required, solves the constrained profiled-baryon problem, and reruns the conditional-envelope fit on the newly fitted core. The envelope's `a1,q20` case is the declared comparison; the original optimization body additionally runs its configuration grid. Grids execute separately per galaxy; any grid failure explicitly fails that galaxy's envelope slot while preserving other galaxies and replicates. Model failures remain explicit rows with no fabricated finite loss.

Each perturbation changes inner observed velocities by independent Gaussian errors and changes distance/inclination **prior centres** with their stated errors, rejecting draws outside the physical domain. Prediction coordinates stay at the original catalogue distance/inclination. Disk/bulge mass-to-light parameters are reoptimized under the fixed original priors; their prior hyperparameters are not sampled. The original catalogue disk scale remains the envelope's catalogue-coordinate scale; the newly fitted distance/inclination factors enter its fitted-core inputs.

The same draw is passed to every comparator. Distance/inclination/stellar prior centres and bounds are common where applicable. Each model retains its original declared optimizer and multistart budget; this is not a claim that all models use identical algorithms or computational effort. The envelope is conditional on a freshly fitted core rather than a simultaneous joint refit of every core and envelope parameter.

Change `fitting_rules.training_weight` before locking to compare:

- `pooled_points`: sum of inner squared standardized residuals plus squared prior pulls, as in the original empirical fits.
- `equal_galaxy`: mean inner squared standardized residuals plus the **same** squared prior pulls. Implemented by multiplying inner observational errors by the square root of the inner point count for optimization only. Outer errors are unchanged. This changes the likelihood/prior tradeoff; it is not multiplication of the entire objective by a harmless constant.

The default is a **sensitivity ensemble** around already observed data. It is not posterior sampling, a parametric bootstrap around a fitted latent truth, a calibrated uncertainty interval, population resampling, shared stellar-zero-point uncertainty, or a shared physical-mass recalibration. No unmeasured covariance model is invented. Those choices require a scientifically justified and predeclared generative model. Pressure/charged full physical fit generation remains an external backend, not an implemented native or Python pilot.

## Locked input format and independent data

Replace all five input files and update the draft source/acquisition statement to use new observations. The native protocol accepts arbitrary galaxy identifiers and complete split manifests. The empirical pilot itself remains restricted to `retrospective_development` and the original outer-20% suffix convention; genuinely new validation requires an appropriately reviewed external fitter.

| File | Required contents |
| --- | --- |
| `predictors.csv` | `galaxy,index,radius_kpc,sigma_kms,gas_kms,disk_kms,bulge_kms,partition,population_role`; no observed outcome column |
| `training.csv`, `evaluation.csv` | `galaxy,index,observed_kms,sigma_kms`, disjoint and jointly covering every predictor row |
| `metadata.json` | Per-galaxy `D,e_D,inc,e_inc,Rdisk,catalogue_D,catalogue_inc`; empirical original pipeline also uses `distance_method` |
| `split.json` | Array of `galaxy,index,partition,population_role`, exactly matching every predictor key |

Indices begin at zero and are contiguous. Radii are positive and nondecreasing; observational errors are finite and positive; signed component velocities are finite. Inner rows precede outer rows, with at least three inner rows and one outer row per galaxy. Distance/disk scale are positive, inclination is in (0,90], and metadata errors are nonnegative. A galaxy belongs entirely to population `train` or `transfer`; this population split is distinct from radial `inner`/`outer`.

`evaluation_galaxies` defines the entire scoring sample. `models` and `primary_model` lock the family. `shared_parameter_models` and `shared_training_galaxies` declare any shared calibration explicitly; shared training galaxies must have population role `train`. Local state can use a transfer galaxy's inner rows; shared parameters cannot. The default empirical models have no population-shared fitted parameters. For a strict transfer-only study, select only transfer galaxies for evaluation and separately declare the training population.

Locking hashes the protocol and every input, validates complete unique splits, and rejects linked/aliased files and path escapes. Every later command checks exact lock coverage and all hashes. Locks are local integrity records, not trusted timestamps or an enforceable preregistration service; someone able to replace protocol and lock together can create another apparent history.

`prospective_candidate` requires `previously_inspected:false`, a dataset source, and an acquisition statement. This is a user declaration. `independence_verified` remains false even for that label. These tools cannot establish acquisition provenance, absence of prior inspection, or observational independence. Hashing the existing SPARC observations does not turn them into prospective evidence.

## Prediction artifact contract

`predictions.csv` has `replicate,galaxy,index,model,predicted_kms,fit_failed`. Replicate zero is the baseline. Supply exactly every combination of replicate, declared model and selected galaxy's outer row, including failures; duplicates, omitted rows, undeclared keys and nonfinite successful predictions are rejected. Any failed row makes that model/galaxy failed for scoring.

The accompanying JSON must include:

- `protocol_sha256`, `predictions_sha256`, and `fitting_rules_json` (copy the latter exact string from the lock).
- `read_evaluation_outcomes:false`, `training_keys` containing every selected evaluation galaxy's inner key exactly once, and `replicates`.
- `ensemble_kind` matching `ensemble_contract.kind`: `baseline_only`, `nuisance_sensitivity`, or `full_fitted_procedure`.
- `has_shared_parameters` consistent with the locked model declarations, and `shared_training_keys` exactly covering the declared shared-training galaxies' inner rows, with no transfer or outer rows.
- For shared fits: `shared_parameter_fits`, one record per replicate/shared model, each containing `replicate,model,training_keys,parameters_file,parameters_sha256,transfer_parameters_sha256`. Parameter files are local to the artifact, and training/transfer hashes must agree.

A `full_fitted_procedure` contract additionally requires at least two replicates, a declared shared training population, and `replicate_provenance` for every replicate. Each record contains:

```json
{
  "replicate": 0,
  "population_multiplicities": {"training_galaxy_A": 1, "training_galaxy_B": 1},
  "refitted_stages": ["every stage declared in ensemble_contract.required_refit_stages"],
  "paired_models": ["every model declared in the protocol"],
  "measurement_draws_file": "replicate_0/draws.json",
  "measurement_draws_sha256": "SHA256",
  "fit_records_file": "replicate_0/fits.json",
  "fit_records_sha256": "SHA256"
}
```

Population multiplicities name every declared training galaxy, including zero counts, and preserve the population size; baseline counts are all one. These common replicate records establish the claimed pairing across comparators. The external backend must actually rerun affected local and shared fits under that replicate and freeze each shared state before transfer. The native contract verifies completeness/hashes and declarations; it cannot prove the content of an optimizer's evidence files or that an external process avoided reading outcomes. External source versions, solver convergence, constraint checks, weighting implementation and generative assumptions still need scientific review. No calibrated full-procedure interval is manufactured from arbitrary submitted replicate files.

## Scoring and statistical limits

The report gives each galaxy's mean outer squared standardized error and MAE. Full-sample model means are undefined if any declared galaxy failed; a separately labeled finite subset mean is descriptive only. Baseline comparisons use comparator loss minus primary loss, so positive effects favor the primary. The whole predeclared family remains represented.

Conditional pairs-bootstrap-t/percentile intervals, null-imposed wild tests and Holm correction use the shared native `StatisticalEngine` only for complete baseline comparisons when enabled and above `inference.minimum_galaxies` (default 20). This is a configurable operational gate, not a universal validity threshold. If any declared shared training galaxy also belongs to the evaluation sample, inference is suppressed for the entire family; disjoint transfer results remain conditional on the frozen shared fit. It does not account for shared-training uncertainty or repair correlations between galaxies. Failed or undefined family members contribute one internally to Holm but receive no invented reported p-value. Sensitivity replicates are not treated as independent samples; the two-galaxy pilot reports descriptive effects only. Read the concentration, degeneracy and Monte Carlo diagnostics before interpreting an eligible baseline result.

`ResearchWorkflow.SelfCheck()` exercises synthetic lock tampering, missing hash coverage, noncanonical replicate IDs, derived-score overflow, shared-training inference gates, nonfinite inputs, overlapping splits, metadata/radius errors, false prospective declarations, shared-model flags, shared training leakage, changed transfer parameters, valid shared transfer, and descriptive failure retention. The main console test harness invokes it alongside the inference and run-safety checks.

