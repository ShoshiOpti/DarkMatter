# Plot construction and caption notes

Run from the investigation directory:

    python numerics/pointwise_comparison_2026_09_25/plot_pointwise.py

The script reads only saved point observations/predictions and saved comparison results. It performs no fitting, bootstrap resampling, p-value calculation, or model selection. Numerical plotting tables, source hashes, figure hashes, and checks are saved in plots/. It writes PDF, PNG and SVG for each figure and publishes byte-identical copies in the investigation's figures/ directory.

## pointwise_cluster_comparison

Suggested caption: Point-level accuracy and galaxy-cluster mean contrasts on the original outer block. (a) Empirical cumulative distributions of absolute standardized speed residuals for all 659 observations in 131 galaxies, weighted by 1/(131 n_g) per point so every galaxy contributes equally. The horizontal axis is linear up to 3 and logarithmic thereafter; the full tails remain visible. (b) Mean paired squared-standardized loss differences, fixed envelope minus comparator, with marginal 95% studentized galaxy-bootstrap intervals. Negative values favor the envelope. The baryon ablation uses a separate horizontal scale and the core's nuisance calibration. All predictions are frozen; intervals condition on those predictions and do not include full refitting, shared systematic errors, or the prior model-development history.

Important interpretation: These are marginal bootstrap-t intervals, not simultaneous intervals or inversions of the separate wild-bootstrap hypothesis test. In particular, the NFW interval and wild-test p-value need not give identical threshold decisions. The figure deliberately has no significance stars or independent-point p-values. The separately profiled baryon baseline, with its failed UGC01281 fit, belongs to the labeled finite-subset sensitivity rather than this four-comparison primary figure.

## pointwise_radial_differences

Suggested caption: Every original outer point's paired squared-standardized loss difference, envelope minus scalar core (left) or NFW (right), against radius divided by that galaxy's last measured radius. Each panel contains all 659 points from 131 galaxies; there is no radial binning, fitting, smoothing, or point deletion. Negative differences favor the envelope. The vertical axes are linear between -1 and 1 and logarithmic outside; the panels retain their full individual ranges. Points within each galaxy are correlated. These descriptive point patterns do not supply an independently tested radial law.

## pointwise_covariance_sensitivity

Suggested caption: Sensitivity of equal-galaxy mean paired correlated quadratic losses to an assumed common within-galaxy residual correlation rho=0, 0.25, 0.5 and 0.75. All four primary comparators retain all 131 galaxies and 659 points. Error bars are the saved marginal 95% galaxy-bootstrap-t intervals. Each panel uses its own vertical scale; the zero line separates directions favoring the envelope (negative) or comparator (positive). Rho=0 reproduces the diagonal-error primary score. The common correlation is specified rather than estimated from the data. These are sensitivities, not additional significance tests.

The covariance assumption is R_rho=(1-rho)I+rho 11^T applied to standardized residuals within each galaxy. It is not claimed to recover the unknown actual radial covariance, its distance dependence, or common cross-galaxy systematics.

## Verification

- Each ECDF includes 659 finite point residuals, carries total weight one, assigns each galaxy total weight 1/131, and reproduces the saved equal-galaxy mean absolute residual.
- All four forest means and interval endpoints are read directly from primary_cluster_results.csv. The figure distinguishes the baryon ablation scale.
- Each radial panel contains all 659 paired differences; displayed axis bounds explicitly include every point with transformed-scale padding.
- All 16 covariance effects and interval endpoints are read directly from primary_cluster_results.csv and score_sensitivity.csv.
- plot_checks.json records counts, numerical ranges, intervals, and input hashes.
- The input hashes are checked again after rendering. Published figures have the same SHA256 as the plotted originals.
