"""Finite-sample conservative p-values for the green all-row fit.

No Wilks degrees-of-freedom assertion, radial independence, or tail simulation.
These test the fixed grey conditional mean and stated marginal error scales.
"""
from pathlib import Path
import csv,json,math,hashlib
import test_sparc as t
import numpy as np
from scipy.stats import norm

HERE=Path(__file__).resolve().parent
rows=list(csv.DictReader((HERE/'predictions.csv').open()))
lookup={(r['galaxy'],int(r['index']),r['model']):r for r in rows if r['phase']=='all'}
records=[]
for (galaxy,index,model),b in lookup.items():
    if model!='baryons':continue
    g=lookup[galaxy,index,'new_analytic']
    y,sigma=float(b['observed']),float(b['sigma'])
    assert y==float(g['observed']) and sigma==float(g['sigma'])
    q0=((y-float(b['predicted']))/sigma)**2
    q1=((y-float(g['predicted']))/sigma)**2
    records.append({'galaxy':galaxy,'index':index,'radius_kpc':float(b['r_kpc']),
                    'quality':b['quality_sample']=='True',
                    'physical_null_defined':b['negative_catalogue_baryon_force']=='False',
                    'outer':b['outer']=='True','q_baryons':q0,'q_green':q1,'improvement':q0-q1})

def summary(selected):
    n=len(selected);q0=math.fsum(r['q_baryons'] for r in selected)
    q1=math.fsum(r['q_green'] for r in selected);delta=q0-q1
    # For ANY fitted curve, T=max(0,Q0-Q1)<=Q0.
    # Under H0 E[Q0]=n, no independence or Gaussian assumption required.
    p_moment=min(1.,n/max(delta,1e-300))
    # If each standardized residual has a N(0,1) marginal, arbitrary dependence:
    # {Q0>=delta} subset union_i {|z_i|>=sqrt(delta/n)}.
    logp_gauss=min(0.,math.log(2*n)+float(norm.logsf(math.sqrt(max(delta,0)/n))))
    p_gauss=math.exp(logp_gauss)
    return {'rows':n,'galaxies':len({r['galaxy'] for r in selected}),
            'chi2_baryons':q0,'chi2_green':q1,'delta_chi2':delta,
            'green_parameters':2*len({r['galaxy'] for r in selected}),
            'gaussian_marginal_conservative_p':p_gauss,
            'gaussian_marginal_log10_p':logp_gauss/math.log(10),
            'gaussian_marginal_equivalent_two_sided_Z':float(norm.isf(p_gauss/2)),
            'variance_only_conservative_p':p_moment,
            'variance_only_equivalent_two_sided_Z':float(norm.isf(p_moment/2)),
            'green_Q_per_row':q1/n,
            'interpretation':'Conservative valid p-values / upper bounds on the null improvement tail; not the exact LRT distribution.'}

out={'scope':'Existing green full-data fit versus existing fixed-calibration grey baryon curve.',
     'null':'Conditional mean equals grey curve; each reported sigma is the true marginal standard deviation. No baryonic-template or calibration uncertainty included.',
     'fitting':'The envelope can be selected and fitted to the same data: its nonnegative residual sum ensures T<=Q0 for every dataset.',
     'gaussian_bound':'p_conservative=min(1,2*N*NormalSF(sqrt(T/N))); Gaussian marginal errors, arbitrary dependence.',
     'variance_bound':'p_conservative=min(1,N/T); zero-mean errors with stated second moments, arbitrary dependence and error shapes.',
     'Z_convention':'Two-sided equivalent Z=NormalISF(p/2). These bounds give lower limits on equivalent tail significance under the stated null.',
     'physical_caveat':'The grey circular-speed prediction is undefined at two central UGC01281 rows with negative net baryonic force. Primary statistic uses the 3389 physically defined rows, retaining all 175 galaxies and all outer mass; the green fit is unchanged. The all3391 diagnostic uses the prior explicit zero-speed convention at these two rows.',
     'primary_physical_rows':summary([r for r in records if r['physical_null_defined']]),
     'all3391_zero_speed_convention':summary(records),
     'quality_physical_rows':summary([r for r in records if r['quality'] and r['physical_null_defined']]),
     'unphysical_null_rows':[r for r in records if not r['physical_null_defined']],
     'source_predictions_sha256':hashlib.sha256((HERE/'predictions.csv').read_bytes()).hexdigest()}
# Independent consistency against the delivered full-data summary.
prior=json.loads((HERE/'RESULTS.json').read_text())['summaries']['all175']
for m,key in [('baryons','chi2_baryons'),('new_analytic','chi2_green')]:
    expected=next(x['total_chi2'] for x in prior if x['phase']=='all' and x['model']==m)
    assert abs(out['all3391_zero_speed_convention'][key]-expected)<1e-7
assert len(records)==3391 and out['primary_physical_rows']['rows']==3389
assert len(out['unphysical_null_rows'])==2 and not any(r['outer'] for r in out['unphysical_null_rows'])
assert all(r['q_green']>=0 and r['q_baryons']>=0 for r in records)
for k in ['primary_physical_rows','all3391_zero_speed_convention','quality_physical_rows']:
    a=out[k];assert 0<a['delta_chi2']<=a['chi2_baryons']
    assert abs(2*norm.sf(a['gaussian_marginal_equivalent_two_sided_Z'])/a['gaussian_marginal_conservative_p']-1)<1e-10
    # Same Gaussian p from a direct union probability expression, independent code path.
    assert abs(2*a['rows']*norm.sf(np.sqrt(a['delta_chi2']/a['rows']))/a['gaussian_marginal_conservative_p']-1)<1e-10
(HERE/'full_fit_significance.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
t.csvout('full_fit_significance_rows.csv',records)
print(json.dumps(out,indent=2))
