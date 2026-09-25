# Independent statistical-method review: pointwise rotation-curve comparisons

Date: 25 September 2026. Scope: design/review only; no galaxy fit or paper edited by this review.

## Recommendation

Use the individual outer-radius residuals, retain their magnitudes, and keep galaxies as the sampling units. The main comparison should be the mean paired difference in squared standardized residual over galaxies, with whole-galaxy bootstrap intervals and a studentized null-imposed wild bootstrap. This is a different scientific question from the existing fair-sign test: it measures average size of the improvement, including large gains and large failures. It does not turn 659 correlated radii into 659 independent galaxies.

Use a compact primary family of specified model comparisons, declare the loss, weighting, failure policy, number of replicates and sensitivities before calculating new outcomes, and report all its members. This is a frozen computational specification, not prospective registration: the observations have already informed model development.

A one-line user explanation is: "We now compare how far each predicted point misses the observation, rather than counting only which model wins a galaxy; uncertainty still resamples entire galaxies so neighboring radii are not treated as independent evidence."

## Data and estimand

For fixed saved predictions define

    z_mgj = (V_mgj - V_observed,gj) / sigma_gj
    d_gj = z_Agj^2 - z_Bgj^2
    d_g = sum_j d_gj / n_g
    Delta = mean_g d_g.

Negative Delta favors model A. All finite differences, including exact/practical ties, enter; the practical sign-test tie threshold is irrelevant to continuous loss. There are G=131 galaxies and N=659 outer observations, with n_g=2--23 and median 4.

This is an equal-galaxy estimand: every galaxy is weighted equally, with its observations averaged internally. It is exactly computable from the full pointwise residual matrix; aggregation to d_g loses no information needed by this additive statistic. The old code already computed the pointwise losses before reducing them to a sign. The improvement is to retain loss magnitude in inference, not to claim that pointwise residuals had never been used.

The different pooled-point estimand is

    Delta_point = sum_g n_g d_g / sum_g n_g.

It asks about a randomly chosen observed radius in this catalogue and gives more weight to more densely sampled galaxies. Report it as a sensitivity, not a replacement chosen for its p value. A whole-galaxy bootstrap of this ratio must resample both numerator and denominator together.

The empirical Delta for these fixed 131 galaxies is known from the saved predictions. A population confidence interval additionally assumes independent, sufficiently representative galaxy draws from a relevant population. SPARC is an observationally selected catalogue; bootstrap intervals do not automatically generalize to every galaxy. It is safest to describe them as conditional catalogue-resampling uncertainty and keep the population assumption explicit.

## Primary inference and exact null

The primary null is equality of expected, equally galaxy-weighted loss:

    H0: E[d_g] = 0,

under an independent-galaxy population/sampling interpretation. It is not equality of win probabilities, equality of models, or absence of a physical dark sector. The procedure conditions on the saved fitted predictions and the chosen model family.

1. Whole-galaxy pairs bootstrap: draw G galaxy indices with replacement and carry every observation, all model predictions, and all requested scores together. Compute the mean d in each draw. Report a marginal 95% percentile interval, with an explicit number of replicates and RNG seed. Intervals for multiple rows are not simultaneous unless separately constructed that way.
2. Observed studentized statistic:

       T = mean(d) / sqrt[sum_g (d_g-mean(d))^2 / {G(G-1)}].

3. Null-imposed Rademacher wild bootstrap for the intercept-only galaxy regression d_g=Delta+u_g: under H0 the restricted fit is zero, so the restricted residual is d_g, not d_g-mean(d). Generate d_g* = w_g d_g with independent w_g in {-1,+1}, and recompute the mean and its standard error for each replicate. Use the same w_g across paired model comparisons, preserving their shared-galaxy structure.
4. Two-sided Monte Carlo tail count k=#{|T*| >= |T|}, with p_MC=(k+1)/(B+1). Keep the exact comparison convention in the protocol.

This is a standard asymptotic bootstrap use of cluster signs. It is not an exact randomization/permutation test for a mean-only null unless additional symmetry/exchangeability assumptions hold. Calling it "exact" would be unjustified. A centered multiplier bootstrap is a different, also asymptotic construction; do not mix its residual formula with a claim to implement the restricted wild bootstrap.

For equal galaxy weighting this is an intercept-only heteroskedastic mean problem after aggregation. For point weighting the robust standard error is instead

    se_point^2 = G/(G-1) * sum_g [n_g(d_g-Delta_point)]^2 / (sum_g n_g)^2.

Do not reuse the unweighted standard error for the weighted estimand. Generate whole-galaxy weights for the pointwise weighted wild test too. If that secondary test is not needed, descriptive paired point means plus whole-galaxy ratio-bootstrap intervals are adequate.

Cameron and Miller explain why resampling whole clusters and allowing within-cluster dependence matter. Their restricted wild bootstrap estimates the model under the tested null before generating cluster-weighted residuals. Their treatment also cautions that a few influential clusters can undermine otherwise large nominal cluster counts. [Primary article, Secs. II and VI](https://faculty.econ.ucdavis.edu/faculty/cameron/research/Cameron_Miller_JHR_2015_February.pdf). A more recent primary overview is [MacKinnon, Nielsen and Webb, 2023](https://arxiv.org/abs/2205.03285).

## Review of the written implementation protocol

The subsequently written PROTOCOL.json chooses a pairs bootstrap-t interval as primary and the percentile interval as a sensitivity. This is consistent with the mean-loss question. For each unrestricted pairs resample compute

    T_b = (Delta_b - Delta) / SE_b,
    CI_t = [Delta - q_.975(T_b) SE, Delta - q_.025(T_b) SE].

The quantile order and signs matter. The null-imposed wild bootstrap supplies the p value separately. Bootstrap-t interval inversion, percentile intervals and the wild-bootstrap test are different procedures and can disagree in finite samples; their agreement should not be assumed.

No draw with zero or nonfinite SE should be silently omitted. Count such draws and implement an explicit degeneracy policy. Inspect distributions if variance is concentrated in very few galaxies. The protocol's four estimable full-sample comparisons, feasible-only profiled-baryon sensitivity, fixed covariance grid, and no-new-policy-model selection are appropriate for the bounded request.

## Numerical resolution and reporting

B=99,999 gives a minimum reported p_MC of 10^-5, corresponding to a two-sided Gaussian equivalent Z about 4.42 before multiplicity adjustment. B=199,999 gives 5 x 10^-6 and Z about 4.56. A numerical bootstrap cannot support arbitrarily large "sigma" from zero exceedances. Report k and B; if k=0, explicitly mark Monte Carlo resolution reached. The true bootstrap tail is not known to equal that floor. A binomial Monte Carlo uncertainty interval for k/B is useful when a p value is close to a reporting threshold.

The plus-one convention avoids a zero estimated tail; it does not make this asymptotic bootstrap exact or repair model selection. Phipson and Smyth's primary result concerns random permutation tests; use its finite-resolution lesson without claiming it proves bootstrap calibration here. [Primary paper](https://arxiv.org/abs/1603.05766).

If equivalent Z is retained, define it solely as Phi^{-1}(1-p/2) and report which model the signed effect favors. It is not a measurement signal-to-noise ratio or physical detection significance. Prefer the effect, interval and p in the main table; move Z to the reproducibility table if space is tight.

For a declared primary family use Holm adjustment on the nominal primary bootstrap p values and label the scope. It neither compensates for trying several score functions and keeping the best nor adjusts unknown historical exploration. Robust and covariance sensitivities should all be displayed, with no best-p selection. Marginal confidence intervals need not invert the wild-bootstrap/Holm test and must not be described as doing so.

Inspect all generated statistics: no silently discarded zero-variance replicates, nonfinite scores or failed fits. If every d_g is exactly zero, report Delta=0, interval [0,0], p=1 and Z=0. Preserve explicit handling of degenerate studentization rather than deleting draws.

## Influence and robustness are essential here

Raw squared residuals can be dominated by one or a few galaxies. The prior audit already found individual-galaxy reversals of model-selection means. Alongside the primary mean, report:

- Every leave-one-galaxy-out mean, only as an influence diagnostic, with no galaxy deleted from the primary analysis.
- The largest absolute contributions d_g/G to the sample mean and the largest shares of variance, (d_g-Delta)^2 / sum(d-Delta)^2. An optional concentration count, 1/sum(variance_share^2), is a descriptive diagnostic, not calibrated degrees of freedom.
- A pointwise absolute standardized residual score, mean_g mean_j |z_mgj|, as an interpretable less-outlier-sensitive sensitivity.
- Optionally Huber residual loss with a declared fixed cutoff (for example 2 standardized-error units), or mean_g log(1+L_mg). These answer different questions; the cutoff/transformation must not be tuned to recover a preferred ranking.

A log(1+L) comparison downweights catastrophic fits and retains performance ratios approximately at high loss. Its possible disagreement with squared-loss inference is information about the error distribution, not a reason to select whichever is more significant. A median paired loss change is also descriptive and robust but can be zero because models share exact boundary solutions.

The square-loss difference has a useful pointwise identity:

    d_gj = (V_Agj-V_Bgj)(V_Agj+V_Bgj-2 V_observed,gj) / sigma_gj^2.

Thus an improvement depends on how far the models differ and which side of their midpoint the observation occupies. A tiny victory and a large correction are no longer treated equally.

## Unknown radial covariance: score sensitivity, not invented precision

Whole-galaxy inference already tolerates arbitrary radial dependence for the chosen empirical score, subject to independent galaxies and usual moment/non-dominance assumptions. It does not recover an unknown observational covariance or convert the diagonal loss into a calibrated chi-square likelihood.

One transparent covariance sensitivity uses a declared equicorrelation matrix on standardized errors,

    R_g(rho) = (1-rho) I + rho 11',
    L_mg(rho) = z_mg' R_g(rho)^(-1) z_mg / n_g,

for fixed rho values such as 0, .25, .5, .75. For 0<=rho<1 this is positive definite. Report all values and do not fit rho to residuals and then treat it as known. The identity

    L_mg(rho)
      = Var_population,j(z_mgj)/(1-rho)
        + mean_j(z_mgj)^2/[1+(n_g-1)rho]

separates radial shape mismatch from a coherent standardized offset. Assumed correlation reduces the weight of a coherent offset while increasing the weight of shape variation. That may change a model ranking; it is not paradoxical and is not a universal inflation of error bars. This simple model does not reproduce every distance, inclination, calibration or gas-covariance pattern.

A distance-based exponential covariance is a possible later sensitivity but needs a fixed correlation-length rule and more choices. With only two to four outer points in many galaxies, fitting an independent covariance length per galaxy would be poorly supported.

There is also an exact conditional measurement-error identity. For fixed predictions and Gaussian y with known covariance Sigma_g,

    Var(d_g | predictions) =
      4/n_g^2 * (V_A-V_B)' D_sigma^(-2) Sigma_g D_sigma^(-2) (V_A-V_B).

The paired loss difference is linear in y, even though each loss is quadratic. This could support a future covariance-based test. Here Sigma_g is unavailable and saved predictions inherit uncertain training calibration; treating this diagonal, frozen-prediction formula as a discovery Z would omit those uncertainties. Do not publish it as the primary test.

## Baseline feasibility and scope

The separately profiled baryon model has an infeasible training pipeline outcome for UGC01281. The established all-galaxy loss is infinite. No ordinary finite-mean bootstrap can include infinity and remain a finite continuous-loss test. Preserve that failure and state that any continuous comparison using the remaining 130 galaxies is conditional on baseline feasibility. Report excluded galaxy and outer-point count. Do not clip loss, replace infinity with a convenient number, or silently compare unmatched samples. A bounded alternative could be specified, but it would be a new estimand requiring explicit interpretation.

The core-calibrated baryon ablation has finite outer predictions for all 131 and can enter the full comparison; label it as an ablation, not the optimal baryon-only model.

NFW has its own inner-fitted nuisance calibration. Envelope and compact control share inherited calibration and only refit amplitudes. A paired score comparison can compare these saved algorithms but is not a likelihood-ratio test with a guessed parameter-count difference. AIC/BIC or Wilks chi-square tails are unsuitable replacements on these held-out, nonnested, conditionally calibrated scores.

No new method repairs reused outer data, inherited nuisance leakage into adaptive validation, or missing microscopic field predictions. Neither pairs nor wild bootstrap re-estimates nuisance parameters or the full model-selection history. If re-estimation uncertainty is later wanted, rerun the full inner-fit/selection pipeline inside a galaxy/measurement resampling design that respects the intended target and supplies a defensible error model. That is a different, larger experiment.

A hierarchical residual model could eventually describe galaxy offsets, radial correlation, and heterogeneity jointly. It adds assumptions and identifiable hyperparameters; it is not automatically superior to the transparent cluster analysis with the present sparse outer sampling. The current task can deliver meaningful new results using existing saved predictions without introducing a new physical fit.

## Suggested paper presentation

Lead with one table of continuous effect estimates, whole-galaxy 95% intervals and nominal adjusted bootstrap p values. Keep model mean/median losses and all feasibility information visible. Put the existing win-rate p/Z table in the Reference Note as a complementary fair-sign question. A figure should show per-galaxy continuous differences or their distribution, a confidence interval forest panel, and the covariance/robustness sensitivities. Preserve the existing rotation-curve gallery: predictions have not changed.

The main wording should say "smaller average retained loss" or "nominal comparison of saved predictive errors," not "theory established at N sigma." Report disagreements among win rate, average error, robust scores and covariance sensitivities; they reveal whether the claim concerns frequent small gains or reliable reduction in total discrepancy.
