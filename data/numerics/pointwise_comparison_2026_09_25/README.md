# Pointwise residual comparison

This analysis retains error magnitudes from 659 outer rotation observations, grouped into 131 galaxies. It does not refit, select or alter a galaxy curve. The main paper Section 6.2 and Reference Section R19 report its results. Older sign tests remain a distinct historical comparison.

## Reproduce

From the investigation directory, with NumPy, pandas, SciPy and Matplotlib installed:

    python numerics/pointwise_comparison_2026_09_25/analyze_pointwise.py
    python numerics/pointwise_comparison_2026_09_25/independent_check/verify_pointwise.py
    python numerics/pointwise_comparison_2026_09_25/plot_pointwise.py

These commands use the bundled frozen data and overwrite the corresponding results. The optional lookup of tmp/python_packages in ancestor directories is an author-environment convenience. A normal Python installation suffices. Original input-source extraction is audited separately in data/INPUT_AUDIT.md and may require the parent workspace; it is not needed for the commands above.

## Protocol and random streams

PROTOCOL.json and PROTOCOL.md were specified before the new bootstrap outcomes, after the data and model family had already been explored. The comparison is retrospective. Base seed 2026092501 controls the galaxy-pairs draws; the production implementation uses the deterministic offsets +1 for wild draws and +2 for feasible-profiled-baseline draws. Each uses 99,999 replicates. These stream allocations are documented here for exact reproduction; no protocol, fit, sample or inferential choice was changed in response to an outcome.

Each point contributes its envelope-minus-comparator squared standardized residual difference. Those differences are averaged within galaxy, then galaxies receive equal weight. Negative effects favor the envelope. Whole-galaxy resampling retains within-galaxy dependence, while assuming independent galaxies. No full radial covariance or nuisance-refit uncertainty is estimated.

Primary comparisons: original core, compact Plummer control, matched NFW, and core-calibrated baryon ablation. Marginal 95% studentized pairs-bootstrap intervals and two-sided null-imposed studentized wild-bootstrap tests are different procedures; the intervals are not test inversions. The four mean-loss tests receive Holm adjustment.

## Principal results

| Comparator | Mean squared-score difference | Marginal 95% bootstrap-t interval | Wild p, Holm4 |
|---|---:|---:|---:|
| Core | -5.174 | [-8.451,-3.280] | 0.00004, MC floor |
| Compact control | -5.017 | [-8.323,-3.052] | 0.00004, MC floor |
| NFW | +10.533 | [+1.069,+27.669] | 0.07146 |
| Core-calibrated baryons | -391.808 | [-583.886,-288.322] | 0.00004, MC floor |

The three Monte Carlo floors mean zero exceedances in 99,999 draws, using (1+k)/(B+1) before adjustment. The displayed Z=4.11 is the Gaussian probability equivalent of that numerical convention, not a measured extreme tail or a bound on physical significance. NFW's interval/test disagreement is reported, together with weighting and covariance sensitivities, rather than selecting whichever method crosses a preferred threshold.

Mean absolute velocity error is 12.879 km/s for core, 10.533 for envelope and 8.109 for NFW. Envelope-minus-core MAE is -2.347 km/s with 95% bootstrap-t interval [-3.316,-1.705]. The squared-loss directions survive every one-galaxy omission, but variance is concentrated in a few galaxies.

Separately profiled baryons fail training for UGC01281. The full continuous comparison is not finite and has no reported p-value. Its labeled feasible-only sensitivity retains 130 galaxies/654 points, alongside the original failure.

## Files and checks

- data/: pointwise tidy tables, input snapshots, original-loss reproduction and hashes.
- results/: pointwise differences, galaxy scores, bootstrap summaries, intervals, test counts, covariance/weighting/influence sensitivities and failure handling.
- independent_check/: separate bootstrap implementation and numerical/manuscript review.
- plots/: all plotted values, figures, counts, input/output hashes and visual review.
- METHOD_REVIEW.md and RESULT_INTERPRETATION_REVIEW.md: independent method and interpretation checks.
- PLOT_NOTES.md: exact figure scope and captions.

The independent numerical reproduction passes all 21 comparison groups, maximum saved-output discrepancy 2.73e-12. The historical fitted curves remain unchanged. A stronger prediction test still requires unexploited data or a complete nested refitting and selection procedure with justified observation covariance.
