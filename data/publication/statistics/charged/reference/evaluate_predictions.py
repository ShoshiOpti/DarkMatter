"""Evaluate frozen charged-halo predictions; never fit or select a model.

python evaluate_predictions.py --input ../predicted_rows.csv --output results
See PROTOCOL.json for the frozen retrospective scope and failure policy.
"""
from pathlib import Path
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import hashlib
import json
import math
import sys

HERE = Path(__file__).resolve().parent
for root in HERE.parents:
    if (root / "tmp/python_packages").is_dir():
        sys.path.insert(0, str(root / "tmp/python_packages"))
        break
import numpy as np
import pandas as pd
from scipy.stats import beta, norm

PROTOCOL = json.loads((HERE / "PROTOCOL.json").read_text(encoding="utf-8"))
PRIMARY = PROTOCOL["primary_model"]
COMPARATORS = PROTOCOL["primary_comparators"]
MODELS = [PRIMARY] + COMPARATORS
KEYS = ["galaxy", "index"]
META = ["radius", "observed", "sigma", "is_outer"]
ALIASES = {
    "baryon_only": "inherited_baryons", "baryon_ablation": "inherited_baryons",
    "scalar_core": "original_core", "extended": "envelope",
    "zero_contact": "free_control", "free_halo": "free_control",
    "charged_population": "charged", "persistent_charge": "charged",
}
COL_ALIASES = {"point_index": "index", "radius_kpc": "radius",
               "observed_kms": "observed", "sigma_kms": "sigma",
               "predicted_kms": "predicted"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def boolean(series, name):
    values = series.astype(str).str.strip().str.lower()
    allowed = {"true": True, "false": False, "1": True, "0": False,
               "1.0": True, "0.0": False}
    if not values.isin(allowed).all():
        raise ValueError(f"Invalid boolean in {name}: {sorted(set(values)-set(allowed))}")
    return values.map(allowed).astype(bool)


def holm(probabilities):
    p = np.asarray(probabilities, float)
    assert np.isfinite(p).all() and np.all((p >= 0) & (p <= 1))
    order = np.argsort(p, kind="stable")
    adjusted = np.empty_like(p)
    adjusted[order] = np.minimum(1., np.maximum.accumulate((len(p)-np.arange(len(p))) * p[order]))
    return adjusted


def normal_equivalent(p):
    return float(norm.isf(float(p) / 2))


def sign_counts(left, right):
    wins = losses = ties = joint_failed = 0
    for a, b in zip(left, right):
        af, bf = np.isfinite(a), np.isfinite(b)
        if not af and not bf:
            joint_failed += 1
        elif not af:
            losses += 1
        elif not bf:
            wins += 1
        elif abs(a-b) <= 1.e-4 * (1 + max(a, b)):
            ties += 1
        elif a < b:
            wins += 1
        else:
            losses += 1
    n = wins + losses
    # Exact integer sum avoids subtracting nearly equal tail probabilities.
    p = min(1., sum(math.comb(n, k) for k in range(min(wins, losses)+1)) / 2.**max(n-1, 0)) if n else 1.
    return dict(wins=wins, losses=losses, ties=ties, joint_failures=joint_failed,
                informative_signs=n, p_sign_raw=p)


def validate_input(path, reference):
    raw = pd.read_csv(path).rename(columns=COL_ALIASES)
    required = set(KEYS+META+["model", "predicted"])
    missing = required-set(raw.columns)
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")
    raw["model"] = raw.model.replace(ALIASES)
    ignored = sorted(set(raw.model)-set(MODELS))
    raw = raw.loc[raw.model.isin(MODELS)].copy()
    if set(raw.model) != set(MODELS):
        raise ValueError(f"Missing required models: {sorted(set(MODELS)-set(raw.model))}")
    if raw.duplicated(KEYS+["model"]).any():
        raise ValueError("Duplicate galaxy/index/model rows, possibly an alias collision")
    raw["is_outer"] = boolean(raw.is_outer, "is_outer")
    if "fit_failed" not in raw:
        raw["fit_failed"] = False
        explicit_failure_flags = False
    else:
        raw["fit_failed"] = boolean(raw.fit_failed, "fit_failed")
        explicit_failure_flags = True
    for col in ["index", "radius", "observed", "sigma"]:
        raw[col] = pd.to_numeric(raw[col], errors="raise")
        assert np.isfinite(raw[col]).all(), col
    assert np.equal(raw["index"], np.floor(raw["index"])).all(), "Noninteger point index"
    raw["index"] = raw["index"].astype(int)
    raw["predicted"] = pd.to_numeric(raw.predicted, errors="raise")
    assert (raw.sigma > 0).all() and (raw.radius > 0).all()
    observations = raw[KEYS+META].drop_duplicates().set_index(KEYS).sort_index()
    if observations.index.duplicated().any():
        raise ValueError("Models disagree on observation metadata")
    if len(observations) not in (659, 3034):
        raise ValueError(f"Expected 3034 all rows or 659 outer-only rows, got {len(observations)}")
    for model in MODELS:
        got = raw.loc[raw.model.eq(model)].set_index(KEYS).sort_index().index
        if not got.equals(observations.index):
            raise ValueError(f"Incomplete row coverage for {model}")
    saved = pd.read_csv(reference).rename(columns=COL_ALIASES)
    saved["is_outer"] = boolean(saved.is_outer, "reference.is_outer")
    saved = saved[KEYS+META].drop_duplicates().set_index(KEYS).sort_index()
    assert len(saved) == 3034 and int(saved.is_outer.sum()) == 659
    target = saved if len(observations) == 3034 else saved.loc[saved.is_outer]
    if not observations.index.equals(target.index):
        raise ValueError("Input observation keys differ from the original frozen sample")
    for col in ["radius", "observed", "sigma"]:
        if not np.allclose(observations[col], target[col], rtol=1e-12, atol=1e-12):
            raise ValueError(f"Changed original observation coordinate: {col}")
    if not np.array_equal(observations.is_outer, target.is_outer):
        raise ValueError("Outer/training split differs from original")
    outer_obs = observations.loc[observations.is_outer].copy()
    assert len(outer_obs) == 659 and outer_obs.index.get_level_values("galaxy").nunique() == 131
    outer = raw.loc[raw.is_outer].copy()
    pred = outer.pivot(index=KEYS, columns="model", values="predicted").loc[outer_obs.index, MODELS]
    bad_outer = (~np.isfinite(pred)) | (pred < 0)
    failed = raw.groupby(["galaxy", "model"]).fit_failed.any().unstack().loc[:, MODELS]
    failed = failed | bad_outer.groupby(level="galaxy").any()
    pred = pred.mask(bad_outer)
    return raw, outer_obs, pred, failed, dict(
        original_rows=len(saved), supplied_observation_rows=len(observations), outer_rows=len(outer_obs),
        galaxies=131, explicit_fit_failed_column=explicit_failure_flags, ignored_models=ignored,
        metadata_unchanged=True, split_unchanged=True, input_sha256=sha(path), reference_sha256=sha(reference))


def bootstrap_effects(D, n, draws, seed, wild_seed, wild=True):
    """Columns are contrasts, rows independent galaxy clusters; all finite."""
    D = np.asarray(D, float)
    if D.ndim == 1:
        D = D[:, None]
    G, K = D.shape
    assert G >= 2 and np.isfinite(D).all()
    n = np.asarray(n, float)
    mu = D.mean(axis=0)
    se = D.std(axis=0, ddof=1) / np.sqrt(G)
    exact_zero = np.all(D == 0., axis=0)
    if np.any((se == 0) & ~exact_zero):
        raise ValueError("Nonzero constant galaxy effects make studentized inference undefined")
    rng = np.random.default_rng(seed)
    means = np.empty((draws, K)); tstars = np.empty_like(means); pooled = np.empty_like(means)
    for start in range(0, draws, 2000):
        stop = min(start+2000, draws)
        ids = rng.integers(0, G, size=(stop-start, G))
        x = D[ids]
        bm = x.mean(axis=1)
        bs = x.std(axis=1, ddof=1)/np.sqrt(G)
        if np.any((bs <= 0) & ~exact_zero):
            raise ValueError("Degenerate nonconstant bootstrap draw: do not silently discard")
        means[start:stop] = bm
        tstars[start:stop] = np.divide(bm-mu, bs, out=np.zeros_like(bm), where=bs > 0)
        nn = n[ids]
        pooled[start:stop] = (x*nn[:, :, None]).sum(axis=1)/nn.sum(axis=1)[:, None]
    lo = mu - np.quantile(tstars, .975, axis=0)*se
    hi = mu - np.quantile(tstars, .025, axis=0)*se
    out = []
    for j in range(K):
        d = D[:, j]; centered = d-mu[j]
        denom = float(np.sum(centered**2))
        shares = centered**2/denom if denom else np.zeros(G)
        loo = (d.sum()-d)/(G-1)
        trim = int(np.floor(.1*G)); sorted_d = np.sort(d)
        out.append(dict(effect=float(mu[j]), cluster_se=float(se[j]),
            bootstrap_t_ci_low=float(lo[j]), bootstrap_t_ci_high=float(hi[j]),
            percentile_ci_low=float(np.quantile(means[:, j], .025)),
            percentile_ci_high=float(np.quantile(means[:, j], .975)),
            median_paired_difference=float(np.median(d)),
            trimmed10_paired_difference=float(np.mean(sorted_d[trim:G-trim])),
            pooled_point_effect=float(np.sum(n*d)/np.sum(n)),
            pooled_percentile_ci_low=float(np.quantile(pooled[:, j], .025)),
            pooled_percentile_ci_high=float(np.quantile(pooled[:, j], .975)),
            leave_one_out_min=float(loo.min()), leave_one_out_max=float(loo.max()),
            leave_one_out_reversals=int(np.sum(np.sign(loo) != np.sign(mu[j]))),
            max_variance_share=float(shares.max()),
            top3_variance_share=float(np.sort(shares)[-3:].sum()),
            variance_concentration_count=float(1/np.sum(shares**2)) if denom else None,
            exact_identical_scores=bool(exact_zero[j])))
    if wild:
        rng = np.random.default_rng(wild_seed)
        hits = np.zeros(K, int)
        observed = np.divide(mu, se, out=np.zeros_like(mu), where=se > 0)
        for start in range(0, draws, 5000):
            k = min(5000, draws-start)
            signs = 2*rng.integers(0, 2, size=(k, G), dtype=np.int8)-1
            meanstar = signs@D/G
            variance = ((D*D).sum(axis=0)-G*meanstar**2)/(G-1)
            sestar = np.sqrt(np.maximum(variance, 0)/G)
            if np.any((sestar <= 0) & ~exact_zero):
                raise ValueError("Nonzero degenerate wild-bootstrap draw")
            tstar = np.divide(meanstar, sestar, out=np.zeros_like(meanstar), where=sestar > 0)
            hits += np.sum(np.abs(tstar) >= np.abs(observed), axis=0)
        for j, r in enumerate(out):
            k = int(hits[j]); p = (k+1)/(draws+1)
            r.update(wild_exceedances=k, wild_draws=draws, p_wild_raw=p,
                t_statistic=float(observed[j]), monte_carlo_floor=bool(k == 0),
                raw_MC_95low=float(beta.ppf(.025, k, draws-k+1)) if k else 0.,
                raw_MC_95high=float(beta.ppf(.975, k+1, draws-k)) if k < draws else 1.)
    return out


def evaluate(path, output, reference, draws=None):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    draws = PROTOCOL["bootstrap_draws"] if draws is None else int(draws)
    assert draws >= 99
    raw, obs, pred, failed, audit = validate_input(path, reference)
    gals = sorted(obs.index.get_level_values("galaxy").unique())
    failed = failed.loc[gals, MODELS]
    n = obs.groupby(level="galaxy").size().loc[gals].to_numpy(float)
    z = pred.sub(obs.observed, axis=0).div(obs.sigma, axis=0)
    scores = (z*z).groupby(level="galaxy").mean().loc[gals, MODELS].mask(failed, np.inf)
    absz = z.abs().groupby(level="galaxy").mean().loc[gals, MODELS].mask(failed, np.inf)
    mae = pred.sub(obs.observed, axis=0).abs().groupby(level="galaxy").mean().loc[gals, MODELS].mask(failed, np.inf)
    scores.to_csv(output/"model_galaxy_squared_losses.csv")
    failed.to_csv(output/"model_galaxy_failures.csv")
    model_rows = []
    for model in MODELS:
        finite = ~failed[model].to_numpy()
        lv = scores[model].to_numpy()[finite]
        model_rows.append(dict(model=model, total_galaxies=131, total_outer_points=659,
            failed_galaxies=int((~finite).sum()), failed_outer_points=int(n[~finite].sum()),
            failed_galaxy_names=";".join(np.array(gals)[~finite]),
            finite_galaxies=int(finite.sum()), finite_outer_points=int(n[finite].sum()),
            full_sample_mean_squared_loss=float(lv.mean()) if finite.all() else None,
            finite_subset_mean_squared_loss=float(lv.mean()) if len(lv) else None,
            finite_subset_median_squared_loss=float(np.median(lv)) if len(lv) else None,
            finite_subset_mean_absolute_standardized=float(absz[model].to_numpy()[finite].mean()) if finite.any() else None,
            finite_subset_mean_absolute_kms=float(mae[model].to_numpy()[finite].mean()) if finite.any() else None,
            finite_subset_galaxies_rms_z_above3=int(np.sum(lv > 9)),
            conditioning="all131" if finite.all() else "explicitly_conditional_on_feasible_galaxies"))
    pd.DataFrame(model_rows).to_csv(output/"model_summary.csv", index=False)
    signs = []
    finite_comps = []
    for comp in COMPARATORS:
        signs.append(dict(primary=PRIMARY, comparator=comp, **sign_counts(scores[PRIMARY], scores[comp])))
        if np.isfinite(scores[[PRIMARY, comp]].to_numpy()).all():
            finite_comps.append(comp)
    ph = holm([r["p_sign_raw"] for r in signs])
    for r, p in zip(signs, ph):
        r.update(p_sign_holm6=float(p), Z_sign_holm6=normal_equivalent(p),
                 inference="nominal_retrospective_independent_fair_sign_diagnostic")
    pd.DataFrame(signs).to_csv(output/"sign_comparisons.csv", index=False)
    continuous = {}; contrast_rows = []; point_rows = []; influence_rows = []
    if finite_comps:
        D = np.column_stack([(scores[PRIMARY]-scores[c]).to_numpy() for c in finite_comps])
        results = bootstrap_effects(D, n, draws, PROTOCOL["pairs_seed"], PROTOCOL["wild_seed"])
        continuous.update(zip(finite_comps, results))
    sensitivities = []
    for j, comp in enumerate(COMPARATORS):
        finite = np.isfinite(scores[[PRIMARY, comp]].to_numpy()).all(axis=1)
        excluded = np.array(gals)[~finite].tolist()
        d = (scores.loc[finite, PRIMARY]-scores.loc[finite, comp]).to_numpy()
        if comp in continuous:
            r = continuous[comp]
            r.update(full_sample_continuous_defined=True, primary_p_assigned=True)
        elif len(d) >= 2:
            r = bootstrap_effects(d, n[finite], draws, PROTOCOL["finite_subset_seed"]+j, 0, wild=False)[0]
            r.update(full_sample_continuous_defined=False, primary_p_assigned=False,
                     p_wild_raw=None, wild_exceedances=None, wild_draws=None,
                     reason="Failed fit makes full-sample continuous mean inference undefined; finite-pair sensitivity only")
        else:
            r = dict(full_sample_continuous_defined=False, primary_p_assigned=False,
                     p_wild_raw=None, reason="Fewer than two jointly feasible galaxies")
        r.update(primary=PRIMARY, comparator=comp, n_galaxies=int(finite.sum()),
                 n_points=int(n[finite].sum()), excluded_galaxies=";".join(excluded),
                 conditioning="all131" if finite.all() else "jointly_feasible_subset_no_primary_p")
        r["mean_absolute_kms_difference"] = float((mae.loc[finite, PRIMARY]-mae.loc[finite, comp]).mean()) if finite.any() else None
        r["mean_absolute_standardized_difference"] = float((absz.loc[finite, PRIMARY]-absz.loc[finite, comp]).mean()) if finite.any() else None
        r["pairs_seed_used"] = (PROTOCOL["pairs_seed"] if comp in finite_comps else PROTOCOL["finite_subset_seed"]+j) if len(d) >= 2 else None
        r["wild_seed_used"] = PROTOCOL["wild_seed"] if comp in finite_comps else None
        continuous[comp] = r
        if len(d):
            mu = d.mean(); centered = d-mu; denom = np.sum(centered**2)
            shares = centered**2/denom if denom else np.zeros(len(d))
            for gal, delta, share in zip(np.array(gals)[finite], d, shares):
                influence_rows.append(dict(comparator=comp, galaxy=gal, delta_galaxy_loss=float(delta),
                    variance_share=float(share), mean_without=float((d.sum()-delta)/(len(d)-1)) if len(d)>1 else None,
                    conditioning=r["conditioning"]))
        for gal in gals:
            a, b = scores.loc[gal, [PRIMARY, comp]]
            contrast_rows.append(dict(galaxy=gal, comparator=comp, left_loss=a, right_loss=b,
                delta=a-b if np.isfinite(a) and np.isfinite(b) else None,
                left_failed=bool(failed.loc[gal, PRIMARY]), right_failed=bool(failed.loc[gal, comp])))
        for (gal, idx), row in obs.iterrows():
            a, b = z.loc[(gal, idx), [PRIMARY, comp]]
            valid = bool(not failed.loc[gal, PRIMARY] and not failed.loc[gal, comp])
            point_rows.append(dict(galaxy=gal, index=idx, radius=row.radius, comparator=comp,
                z_primary=a, z_comparator=b, delta_squared=a*a-b*b if valid else None,
                delta_absolute_kms=(abs(a)-abs(b))*row.sigma if valid else None,
                comparison_valid=valid))
        for rho in [0., .25, .5, .75]:
            cscore = []
            for model in [PRIMARY, comp]:
                vv = z[model].groupby(level="galaxy").var(ddof=0).loc[gals].to_numpy()
                mm = z[model].groupby(level="galaxy").mean().loc[gals].to_numpy()
                cscore.append(vv/(1-rho)+mm*mm/(1+(n-1)*rho))
            dd = (cscore[0]-cscore[1])[finite]
            sensitivities.append(dict(comparator=comp, assumed_rho=rho, n_galaxies=int(finite.sum()),
                mean_difference=float(dd.mean()) if len(dd) else None,
                median_difference=float(np.median(dd)) if len(dd) else None,
                conditioning=r["conditioning"], new_p_assigned=False))
    placeholders = [continuous[c]["p_wild_raw"] if continuous[c]["primary_p_assigned"] else 1. for c in COMPARATORS]
    adjusted = holm(placeholders)
    for comp, p in zip(COMPARATORS, adjusted):
        continuous[comp]["p_wild_holm6"] = float(p) if continuous[comp]["primary_p_assigned"] else None
        continuous[comp]["Z_wild_holm6"] = normal_equivalent(p) if continuous[comp]["primary_p_assigned"] else None
        continuous[comp]["inference"] = "nominal_retrospective_fixed_prediction_cluster_mean_loss"
    pd.DataFrame([continuous[c] for c in COMPARATORS]).to_csv(output/"continuous_comparisons.csv", index=False)
    pd.DataFrame(contrast_rows).to_csv(output/"galaxy_contrasts.csv", index=False)
    pd.DataFrame(point_rows).to_csv(output/"pointwise_paired_differences.csv", index=False)
    pd.DataFrame(influence_rows).to_csv(output/"galaxy_influence.csv", index=False)
    pd.DataFrame(sensitivities).to_csv(output/"covariance_sensitivity.csv", index=False)
    audit.update(status="pass", bootstrap_draws=draws, full_protocol_draws=draws == PROTOCOL["bootstrap_draws"],
        protocol_sha256=sha(HERE/"PROTOCOL.json"), evaluator_sha256=sha(__file__),
        primary=PRIMARY, family=COMPARATORS, finite_primary_continuous_comparators=finite_comps,
        pairs_seed=PROTOCOL["pairs_seed"], wild_seed=PROTOCOL["wild_seed"],
        finite_subset_base_seed=PROTOCOL["finite_subset_seed"],
        finite_subset_seed_rule="base plus zero-based comparator index",
        no_refitting=True, no_model_selection=True, no_new_validation=True,
        continuous_failed_slots_are_adjustment_placeholders_only=True)
    (output/"checks.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (output/"PROTOCOL_USED.json").write_text(json.dumps(PROTOCOL, indent=2), encoding="utf-8")
    print(pd.DataFrame(signs).to_string(index=False))
    columns = ["comparator", "n_galaxies", "effect", "bootstrap_t_ci_low", "bootstrap_t_ci_high", "p_wild_holm6", "Z_wild_holm6"]
    print(pd.DataFrame([continuous[c] for c in COMPARATORS]).reindex(columns=columns).to_string(index=False))
    return audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE/"results")
    parser.add_argument("--reference", type=Path, default=HERE.parents[1]/"pointwise_comparison_2026_09_25/data/pointwise_all.csv")
    parser.add_argument("--draws", type=int, help="Override for testing only; nonprotocol draw count explicitly flagged")
    args = parser.parse_args()
    evaluate(args.input, args.output, args.reference, args.draws)

