# Dark Universe: consolidated .NET analysis

Native .NET 10 application with Math.NET Numerics 5.0.0 and System.CommandLine
2.0.12. Runtime reproduction needs no Python. Package versions are locked.

## Commands

Run from this directory with `dotnet run -c Release -- COMMAND` or from the root
with `dotnet run --project DarkMatter -c Release -- COMMAND`.

| Command | Work performed |
|---|---|
| `publication` (default) | Three statistical families, inverse/parent diagnostics, persistence audit, 25 guarded publication figures and three current-result panels |
| `verify` / `all` | Current workflow plus SPARC checks, historical statistics and 47 older plot variants |
| `statistics` | Current and historical saved-prediction statistics |
| `diagnostics` | Inverse scalar closure, parent compatibility algebra and persistence-record audit |
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
