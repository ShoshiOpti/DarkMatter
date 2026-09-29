# Dark Universe: consolidated .NET analysis

Native .NET 10 application with Math.NET Numerics 5.0.0 and System.CommandLine
2.0.12. Runtime reproduction needs no Python. Package versions are locked.

## Commands

Run from this directory with `dotnet run -c Release -- COMMAND` or from the root
with `dotnet run --project DarkMatter -c Release -- COMMAND`.

| Command | Work performed |
|---|---|
| `sparc-full` | Native full-profile refits for all 175 SPARC galaxies, updated Figure 3, all-galaxy atlas, grouped outer statistics and full-fit conservative p/Z bounds |
| `publication` (default) | Three statistical families, inverse/parent diagnostics, reconstruction checks, 25 guarded publication figures and current-result panels |
| `verify` / `all` | Publication workflow, historical checks/plots, and the latest full-SPARC native refits, plots and statistics |
| `statistics` | Current and historical saved-prediction statistics |
| `diagnostics` | Inverse scalar closure, parent compatibility algebra, persistence audit and reconstructed-field checks |
| `reconstruction` | Parent/polar equation and momentum checks, internal-stress response and frozen-profile density/force audit |
| `formation --experiment PATH` | Hash-checked conditional formation experiments, prescribed-design descriptive statistics and new SVG plots |
| `plots` | Current statistics/diagnostics, 25 publication scenes and 47 older SVG variants |
| `pointwise` | Original 25 September empirical comparison and its three native figures |
| `import` | SPARC import and validation |
| `calibrate` | Synthetic null/power and interval-coverage experiments with Monte Carlo uncertainty |
| `resolution` | Charged/control differences, explicit numerical-error and practical-margin assessment |
| `research-create` | Editable development protocol with separate predictors, training and evaluation files |
| `research-lock` / `research-validate` | Validate/hash-lock the declared study inputs and fit policy |
| `research-pilot` | Actual empirical Python refits under paired perturbations in a staged working copy |
| `research-evaluate` | Evaluate newly supplied predictions and fitter provenance against the locked study |

Every command requires a complete input manifest and checks it again at the end.
`--data PATH` selects inputs; `--out PATH` selects an output root. Canonical paths
must not overlap, including through Windows junctions. Unsafe linked input/output
entries are rejected. `--help` lists command-specific options. Defaults are this
project's `data/` and `output/`, or those beside a published executable.

Each invocation gets a fresh `output/runs/<UTC-time>-<id>/`. Its report records
status, stage, timings, executable hash, inputs and errors. Successful runs alone
replace `output/latest.json` and the `output/index.html` redirect; failed attempts
remain visible through `output/latest-attempt.json`. Existing results are kept.
These pointer files are replaced atomically per file with rollback on failure;
they are not a transactional multi-file database. Older direct `output/` files
are historical and should not be mistaken for the latest run.

## Results and scientific scope

### Full SPARC update (29 September 2026)

```powershell
dotnet run --project DarkMatter -c Release -- sparc-full
```

From inside `DarkMatter/`, omit `--project DarkMatter`. Open `DarkMatter/output/index.html`
(or the `index.html` under a supplied `--out` directory). `publication` retains the
historical workflow. `verify` and `all` also run the new full-SPARC analysis and link
both sets of results from the gallery. Existing scientific results are preserved.

The `sparc_full/` run directory contains 32 freshly rendered SVG figures, including
the same six Figure 3 examples and 30 atlas sheets covering every galaxy. Every SVG
has a `.series.csv` with the plotted scientific coordinates. `predictions.csv`,
`fit_states.csv`, `galaxy_scores.csv`, `RESULTS.json`, `full_fit_significance.json`,
`nuisance_sensitivity.json` and `verification.json` contain the numerical results.

The console recomputes positive density interpolation, exact volume integrals,
infinite exterior tails and full central-potential constraints. It performs 1,400
native mass/length fits (four halo templates, two radial protocols, 175 galaxies).
No Python, plotting service, network or repository-relative reference tree is
needed at runtime. The existing SVG renderer creates all figures from this run's
predictions, rather than replaying frozen pictures.

The numerical **field profiles are fixed theory inputs**, supplied with their
original solver, equations and refinement diagnostics in `data/sparc_full/`.
This command does not solve baryon-forced field states or fit universal microscopic
parameters. The matched distance/inclination/stellar sensitivity uses archived
nuisance-fitted predictions: .NET recomputes their scores and resampling, but does
not rerun that nuisance optimizer. The original Python regeneration sources and
requirements are included as provenance, separate from native runtime needs.

Primary group-bootstrap draws are generated natively using PCG64 and sequential
conditional-binomial inversion. Complete count-stream SHA-256 contracts match
NumPy's 99,999 draws for both sample definitions; outcome statistics are calculated
afresh. Every run verifies native profile values, every velocity, inference counts,
adjusted p/Z and green-fit bounds against independent Python reference results.
Disagreement fails the run before the latest successful gallery is promoted.

The green all-data fit's significance is reported as **conservative p upper bounds
and equivalent two-sided Z lower bounds**, distinct from outer-prediction bootstrap
and sign tests. Gaussian marginal errors give p<=2.60295e-121 and Z>=23.4211;
assuming only correct marginal variances gives p<=0.00177061 and Z>=3.12623.
Both allow arbitrary error correlations but condition on the fixed grey means
and quoted error scales. They are not probabilities that the action is correct.
Two physically undefined central grey predictions are flagged; the physical bound
uses 3,389 rows across all 175 galaxies, and the 3,391-row zero-speed convention is
reported separately. Exterior halo mass is never truncated.

`tools/stage_sparc_full.py` is a developer-only staging script for deliberate input
updates; it preserves and verifies unrelated entries in the complete input manifest.
Its random-count contract uses the standard conditional-binomial algorithm (see
[NumPy distributions source](https://github.com/numpy/numpy/blob/v2.3.5/numpy/random/src/distributions/distributions.c)).

`verification/sparc_full_upgrade_validation.json` records the 43 passing regression
tests, standalone published execution, all 175 atlas panels, 16,955 checked model
coordinates, and 73 byte-identical result files across source/published runs.
`tools/check_sparc_full_delivery.py` is the optional standard-library-only developer
audit that checks those output identities; it is not required by the console.

### Historical workflow outputs

Paths below are relative to a successful run directory.

| Output | Content |
|---|---|
| `index.html` | Current gallery and links to reports/coordinate tables |
| `run_report.json` | Overall verification, runtime versions, scope and consumed-input hashes |
| `publication_statistics/{envelope,pressure,charged}/` | Recomputed exact signs, continuous effects, paired intervals, influence and covariance diagnostics |
| `pointwise_comparison/` | Original empirical pointwise analysis, including its primary continuous bootstrap |
| `publication_diagnostics/` | Recomputed inverse and parent algebra, accepted persistence records and checks |
| `publication/figures/` | 25 SVGs, `.series.csv` scientific artist coordinates and `.axes.json` metadata |
| `publication/figure_replay_report.json` | Scene hashes, producer/input provenance, geometry and element checks |
| `publication/result_compatibility.json` | Original input hashes and 15 recomputed tables checked against the figure baseline |
| `current_results/` | Three freshly rendered effect/interval panels, plotted series and readable comparison tables |
| `reconstruction/` | Combined theory/plot gallery, audit report and (publication runs) statistical-impact record |
| `field_reconstruction/` | Native parent/polar equations, same-density different-stress fixture, CSVs and SVGs |
| `reconstruction_audit/` | Exact density/projection and force bookkeeping for saved profiles; no new fit |
| `formation/` | Outcome-source checks, design-weighted descriptive statistics and conditional formation SVG/CSV/HTML reports |
| `research/` | Calibration, resolution, protocol or refit results for the selected research command |
| `historical_figures.html` | Additional earlier workflows after `all`, `verify` or `plots` |

The envelope, fixed-projector pressure and charged statistical families remain
separate. Continuous inference uses 99,999 galaxy-level pairs/wild draws, the
original seeds and fixed Holm families. The profiled-baryon failure is retained in both physical comparisons: 130 galaxies / 654 points form the explicitly restricted diagnostic;
no full-sample continuous p-value is invented. Failed predictions and exact/tight
sign conventions are audited. See `../CONTEXT.md` for interpretation.

Inverse diagnostics recalculate scalar potential/support, mass feasibility,
minimax errors, contact reduction, nodes and the parent necessary equation from
retained density derivatives or amplitude-Laplacian terms. Saved derivative
uncertainty/domain certificates are inputs. The persistence audit summarizes 14
accepted runs in seven galaxies and preserves six rejected numerical attempts.
Neither calculation is fresh differentiation, fitting or evolution.

The 25 publication figures use compact frozen vector scenes exported from the
original Matplotlib producers, including axis metadata and 71,902 scientific
coordinate rows. Native code validates and renders those scenes. Before publishing the gallery,
it checks the original input identities and numerical agreement of 15 newly
computed tables with the explicit `data/publication/figures/result_contract.json`
baseline. Incompatible results fail the run instead of showing stale figures as
current. Deliberate scientific changes require regenerated scenes and a reviewed
updated contract. The separate `current_results/` panels are built directly from
this run's statistics, with feasible-only comparisons and undefined tests marked. `tools/extract_publication_scenes.py` and retained
original producers in `../reference/` support deliberate regeneration. Original
PDF figure baselines remain at `../figures/` for manuscript builds.

SPARC gas includes helium, signed negative components represent signed force,
undefined circular speeds remain undefined, and observations/errors are never
rescaled twice. See the root context document for distinctions between the
ordinary-GR models, empirical envelope and inverse diagnostics.

## Build, publication and verification

From the project root:

```powershell
dotnet restore DarkUniverse.slnx --locked-mode
dotnet build DarkUniverse.slnx -c Release --no-restore
dotnet run --project DarkMatter.Tests -c Release --no-build
dotnet run --project DarkMatter -c Release --no-build -- verify
dotnet publish DarkMatter -c Release --no-restore -o output/native
```

The publish output includes all required `data/` inputs. With a .NET 10 runtime,
`dotnet output/native/DarkUniverse.Console.dll publication` works independently
of the source/archive directory. Python science regeneration is a separate
workflow documented in `../reference/README.md`.

`DarkMatter.Tests` is a dependency-light regression executable: use `dotnet run`,
not `dotnet test`. It exits nonzero if any safety/scientific regression fails.

`verification/research_upgrade_validation.json` records the current upgrade checks.
`verification/consolidated_validation.json` records the earlier migration checks.
Other pre-existing files in `verification/` are historical evidence and may
state older figure/package counts. They are not current completion certificates.

## Research workflow

```powershell
dotnet run --project DarkMatter -c Release --no-build -- calibrate --protocol DarkMatter/research/calibration.example.json
dotnet run --project DarkMatter -c Release --no-build -- resolution --protocol DarkMatter/research/resolution.example.json
dotnet run --project DarkMatter -c Release --no-build -- research-create
```

Run these examples from the solution root. Calibration and resolution protocols
are typed JSON settings; unknown fields are rejected. Defaults use 500 simulated
datasets per scenario/case and 999 resampling draws, independently of the frozen
99,999-draw publication protocol. Resolution margins default to null and therefore
produce no equivalence claim.

For the study/fitter contract, locking and Python pilot commands, see
[research/README.md](research/README.md). For statistical assumptions, simulation
interpretation and numerical budgets, see
[research/STATISTICAL_METHODS.md](research/STATISTICAL_METHODS.md). The native
publication workflow remains Python-free; the optional pilot requires Python
and the scientific dependencies described there. `tools/` and `research/` are
included when publishing, but the original `reference/` tree must be supplied
separately to a published pilot via `--reference`.

## Reconstructed dynamics and current plots (27 September 2026)

`reconstruction` checks the original two-field quartic GPP equations against
an independent polar-variable implementation on 64/128 Fourier grids. A
periodic fixture has the same initial density and zero bulk current for three
relative-phase amplitudes, but different internal energy and density acceleration.
This is an instantaneous equation check, not time integration or a halo fit.
The numerical tolerances are explicitly separated for ordinary equations and
roundoff-sensitive fourth spatial derivatives.

A second audit distinguishes normalized total-density reconstruction from a
projection that discards orthogonal density. The inverse fitted target projects
to the rescaled core, alpha_c times the original core. Its empirical-minus-original
force also contains the core-amplitude change. Two unchanged archived fit tables
are now native inputs, bringing the complete data manifest from 211 to 213 files;
all prior input hashes and the 25-figure baseline contract are unchanged.

The publication run recomputes all three original statistical families and writes
`reconstruction/statistical_impact.json` only after the result contract passes.
There is no new significance test of an algebraic change of coordinates. A new
predictive rotation curve requires internal initial/boundary data or a justified
stationary-state selection and a forward solution.

## Conditional formation experiments (27 September 2026)

```powershell
dotnet run --project DarkMatter -c Release -- formation --experiment research/population_law_2026_09_27/formation_bundle.json --out output/formation
```

This additive command consumes a strict schema-1 JSON bundle and the hashed JSON
sources it names. Each plotted outcome must equal its declared source value at
an explicit JSON pointer. Unknown fields, duplicate properties, nonfinite data,
invalid weights and mismatched hashes are rejected. Sources and the bundle are
checked again after report generation. The ordinary 213-file native data manifest
is also verified; no frozen publication inputs or figure contracts are changed.

The command recomputes the weighted mean, population standard deviation and
range **over the supplied design cases**, separately for every study and group.
Design weights are prescribed inputs. These descriptors and their cumulative
plots are not galaxy population frequencies, uncertainty intervals or statistical
significance tests. A group with one case correctly has zero design dispersion,
which carries no uncertainty claim. Numerical diagnostics are carried into the
summary as their maximum absolute values. The formation integrators and their
numerical convergence checks are independent scientific producers; this native
command verifies provenance and renders their results without requiring Python.

The bundle contract is documented in
`../research/population_law_2026_09_27/native/INPUT_SCHEMA.md`.
