# Significance of the green all-data envelope versus the grey fixed baryonic curve

This addendum directly evaluates the existing green curves, fitted to all radial observations, against the existing grey catalogue-baryon predictions. It supplies **finite-sample conservative p-values**, rather than an exact likelihood-ratio tail or a simulation-based estimate. It does not reuse the outer-prediction sign statistic.

## Numerical result

All 175 galaxies are retained. The primary statistical calculation uses 3,389 physically defined grey predictions. At two central UGC01281 observations the catalogue's net baryonic force is negative, so the baryon-only circular speed is undefined; those two observations cannot supply a physical Gaussian null mean. They remain in the data and the original green fit is unchanged. No galaxy, outer observation or exterior halo mass is removed.

| Quantity | Value |
|---|---:|
| Baryon squared standardized residual, Q0 | 1,960,465.195407 |
| Green-envelope squared standardized residual, Q1 | 46,434.777028 |
| Observed improvement, T = Q0 - Q1 | 1,914,030.418379 |
| Fitted halo coordinates | 350 (mass and scale per galaxy) |
| Conservative p-value with Gaussian marginal errors | 2.6029476401e-121 |
| Equivalent two-sided Gaussian Z | 23.42114030 |
| Conservative p-value using marginal variances only | 0.001770609269 |
| Equivalent two-sided Gaussian Z | 3.12623346 |

These valid conservative p-values are **upper bounds on the actual null tail probability of the improvement statistic**, so their Z conversions are lower bounds on the corresponding tail significance. They are not measured tail probabilities of exactly those magnitudes.

Keeping all 3,391 observations and explicitly using the previous nonphysical zero-speed continuation at the two flagged rows gives Q0=1,960,466.562888, Q1=46,442.019108 and T=1,914,024.543780. Its Gaussian bound is p<=3.0800567444e-121 (Z>=23.41396627), and its variance-only bound is p<=0.001771659622 (Z>=3.12605907). This is a convention-dependent diagnostic, not a circular-orbit solution at those two radii.

## Null, statistic and calibration

The null is that the **fixed plotted grey curve** gives the conditional mean velocity at each catalogue radius, and each tabulated error is the true marginal standard deviation. The baryonic force templates, distances, inclinations and stellar mass-to-light ratios are conditioned on as fixed, exact inputs. This excludes unmodelled baryonic or calibration uncertainty. It tests that specified baryon-only prediction, not all possible baryon-only GR models.

Define

\[
Q_0=\sum_{i=1}^N\frac{(V_i-\mu_{0i})^2}{\sigma_i^2},\qquad
Q_1=\sum_{i=1}^N\frac{(V_i-\widehat\mu_{1i})^2}{\sigma_i^2},\qquad
T=\max(0,Q_0-Q_1).
\]

Under independent Gaussian errors and exact optimization, Q0-Q1 is twice the log-likelihood improvement. With unknown correlation, we treat it as a standardized squared-residual improvement statistic, not the full joint likelihood ratio.

For **any** fitted or selected alternative, including a fit selected on these same data,

\[
0\leq T\leq Q_0,
\]

because Q1 is nonnegative. Consequently a null-tail bound for Q0 also bounds the null tail of T. This remains true for all 350 fitted coordinates, a boundary at zero halo mass, an unidentified halo scale under the null, and even a hypothetical perfect fit with Q1=0. It is deliberately more conservative than calibrating the actual fit family. It does not use Wilks' theorem or sqrt(delta chi-square).

### Gaussian marginal errors, arbitrary dependence

If every standardized null residual Zi=(Vi-mu0i)/sigmai has a standard-normal marginal distribution, then

\[
\{Q_0\geq t\}\subseteq\bigcup_{i=1}^N
\{|Z_i|\geq\sqrt{t/N}\}.
\]

The union bound therefore gives

\[
\Pr_0(T\geq t)\leq\Pr_0(Q_0\geq t)
\leq\min\{1,\,2N\overline\Phi(\sqrt{t/N})\}.
\]

At the observed T, this yields the Gaussian p-value in the table. **No independence assumption is needed**, either between radial observations or between galaxies. The statement nevertheless requires Gaussian marginal tails and correct error scales; it does not establish those assumptions empirically. The standard union inequality is described by [NIST](https://itl.nist.gov/div898/handbook/prc/section4/prc473.htm); its application to this fitted-improvement statistic is the derivation above.

### Only second moments, arbitrary distribution and dependence

If E0[Zi]=0 and E0[Zi^2]=1, then E0[Q0]=N without independence. Markov's inequality gives

\[
\Pr_0(T\geq t)\leq\Pr_0(Q_0\geq t)\leq\min\{1,N/t\}.
\]

This yields p<=0.001770609269. It permits non-Gaussian tails and arbitrary correlations, but still assumes the quoted variances and the fixed null means are correct. It is a much weaker bound, not a contradictory result.

For either bound the reported convention is

\[
Z_{\rm eq}=\Phi^{-1}(1-p_{\rm conservative}/2).
\]

Because both bounds are monotone in T and dominate its actual null survival function, evaluating them at observed T gives conservative valid p-values. No Monte Carlo floor is involved. The two calculations are nested assumption scenarios for the same comparison, not separate discoveries.

## Interpretation

The green fit improves substantially over the specified grey prediction. Under the stated null-error assumptions, the improvement remains significant even when calibration is made conservative enough to allow a perfectly flexible alternative.

This rejects a specified baryon-only mean/error model; it does **not** assign a probability that the modified action is correct. The bound is intentionally insensitive to which alternative produced the improvement. Comparing the two-field profile against conventional dark halos, enforcing shared microscopic scales, and propagating baryonic/calibration uncertainty are needed for more specific physical conclusions. The former matched-nuisance outer-prediction result is a different experiment and cannot be replaced by the numbers here.

The green fit is also not a perfect description at the quoted error scale: Q1/N=13.7016. Its much smaller residual relative to baryons should not be confused with residuals consistent with unit variance or with a globally verified physical solution.

## Reproduction and checks

Run `python -B full_fit_significance.py` in this directory. It reads the delivered row-level predictions without refitting or altering any galaxy profile and writes `full_fit_significance.json` and `full_fit_significance_rows.csv`. It independently reconstructs all sums, checks them against the original full-data results, verifies the two flagged rows are inner observations, and checks both numerical Gaussian-tail and Z conversions. The input prediction SHA-256 is recorded in the JSON. Full exterior mass treatment is inherited unchanged from those predictions.
