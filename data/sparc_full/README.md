# Full SPARC complete-mass profile test

Completed 29 September 2026. Coverage: all 175 SPARC galaxies and 3,391 measurements; no poor-fitting galaxy removed. The historical quality subset (Q < 3, inclination >= 30 degrees, at least eight rows) contains 131 galaxies. The report revisits the same six Figure 3 examples; the companion atlas includes every galaxy.

## Scope

This is a retrospective, conditional test of isolated spherical halo templates from the supplied two-field/Friedmann extension. Each halo fits its total mass and length, with shape fixed independently of these new fits. It is not a common-particle-mass fit or a self-consistent baryon-forced two-field equilibrium for each galaxy. The cosmological solution alone does not specify galaxy occupations. All models use weak-field GR; the null is GR sourced by catalogued baryons alone.

The halo extends to infinity, with enclosed mass used for radial force and the exterior potential integral retained for the central-potential bound. SPARC's signed baryonic force contributions are preserved. No truncation is introduced at the last observed radius. Some enormous extrapolated outer masses are weakly identified; they are not measured galaxy inventories.

## Statistical results

Primary evaluation fits the inner N-max(2,ceil(0.2N)) radii and scores the remaining outer rows (761 total). The loss averages squared velocity residual/error first within each galaxy and then equally across galaxies. All-row fits are also provided, but their fitted residual improvements are not assigned predictive p-values.

For all 175 galaxies at fixed catalogue distance, inclination and disk/bulge mass-to-light ratios 0.5/0.7, baryon loss is 461.199204 and the new analytic-envelope loss is 44.963256. The four-comparison Holm-adjusted Monte Carlo p is 0.00004, equivalent to two-sided Z=4.10748. This value is a resampling-resolution floor (zero exceedances in 99,999 draws), not an exactly measured extreme tail probability. The quality-subset values are p=0.00008 and Z=3.94440.

Allowing the same distance, inclination and common stellar multiplier to vary with the same priors in both models gives all-sample losses 81.940538 and 25.558056, exploratory unadjusted p=0.01083 and Z=2.54814. Quality-subset p=0.01453 and Z=2.44389. There is no single calibration-independent significance. The response-versus-scalar improvement is small; the new analytic approximation and old exact dPIE formula give practically identical fits under the same calibration.

The bootstrap uses calibration groups, including a single Ursa Major distance-method-4 group, rather than treating all radial measurements as independent. Missing radial covariance, historical hypothesis search and unknown halo occupations remain limitations. Exact group-sign diagnostics answer a different question and must not be advertised as a detection of two-field physics.

UGC01281 has two inner rows with negative net catalogue baryonic force. The flagged zero-speed continuation is a descriptive/optimization convention, not a physical circular solution. No fixed-calibration outer score uses that convention. The galaxy is retained, and a separately labelled admissible-subset sensitivity is supplied.

## Reproduction

Use Python 3 with NumPy, SciPy and Matplotlib; pypdf is used only for delivery checks. The report renderer additionally requires pdfLaTeX on PATH. From this directory:

```text
python -B build_profiles.py
python -B test_sparc.py
python -B robustness.py
python -B audit.py
python -B check_mesh.py
python -B render_report.py
```

All inputs needed for these commands are in `inputs/`. Local project runtime paths in the scripts are optional import fallbacks; ordinary installed packages work without them. Outside the original project, the renderer writes a sibling `report/` directory. The existing CSV/JSON outputs allow figures and the report to be regenerated without refitting. Repeating the 99,999-draw inference uses the recorded seed; small floating-point differences can depend on numerical-library versions.

## Files

- `PROTOCOL.json`: declared scope, fitting and inference rules.
- `RESULTS.json`: complete primary summaries, intervals, p/Z and sign diagnostics.
- `nuisance_sensitivity.json`, `guard_sensitivity.json`, `mesh_sensitivity.json`: calibration, potential-bound and resolution checks.
- `fit_states.csv`, `predictions.csv`, `galaxy_scores.csv`: all fitted states and row/galaxy results.
- `full_nuisance_*.csv`: matched distance/inclination/stellar fits and predictions.
- `field_profiles.csv`, `profile_diagnostics.json`: independently reconstructed field profiles and numerical diagnostics.
- `independent_audit.json`, `verification.json`: mass integrals, prediction reconstruction and protocol checks.
- `input_manifest.json`: portable input paths, SHA-256 hashes and source locations.
- `figure3_updated.png`, `full_profile_summary.png`: report figures.

SPARC source: Lelli, McGaugh & Schombert (2016), AJ 152, 157, https://arxiv.org/abs/1606.09251 and https://astroweb.case.edu/SPARC/ . User-supplied equations and diagnostics and the previously archived theory-only compression are included with their hashes.
