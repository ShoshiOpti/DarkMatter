from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2] / 'student_plot_code')
ap.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'output')
a = ap.parse_args(); root = a.source / 'numerics'; out = a.output / 'figures'
frames = {}; inputs = {}; checked = {}; numerical_values = 0; max_relative_error = 0.

def read(path):
    path = root / path
    inputs[str(path.relative_to(a.source)).replace('\\','/')] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)

def select(frame, **kw):
    for key, value in kw.items(): frame = frame[frame[key] == value]
    return frame

def plot(name):
    if name not in frames: frames[name] = pd.read_csv(out / (name + '.series.csv')); checked[name] = set()
    return frames[name]

def check(name, panel, series, x, y, **intervals):
    global numerical_values, max_relative_error
    frame = plot(name); z = frame[(frame.panel == panel) & (frame.series == series)]
    if not len(z): raise AssertionError((name,panel,series,'missing series'))
    vals = dict(x=x,y=y,**intervals)
    for key, values in vals.items():
        expected = np.atleast_1d(np.asarray(values,float)); actual = z[key].to_numpy(float)
        if len(expected) != len(actual):
            raise AssertionError((name,panel,series,key,len(expected),len(actual)))
        np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=2e-10, equal_nan=True, err_msg=f'{name}/{panel}/{series}/{key}')
        mask = np.isfinite(expected) & np.isfinite(actual)
        if mask.any(): max_relative_error = max(max_relative_error, np.max(abs(actual[mask]-expected[mask])/np.maximum(1,abs(expected[mask]))))
        numerical_values += len(expected)
    checked[name].add((panel,series))

def ecdf(name,p,s,values,origin=False):
    values = np.asarray(values,float); x = np.sort(values[np.isfinite(values)]); y = np.arange(1,len(x)+1)/len(values)
    if origin: x = np.r_[x[0],x]; y = np.r_[0,y]
    check(name,p,s,x,y)

def observed(name,p,s,rows,r='radius_kpc',v='observed_kms',e='error_kms',train=None):
    masks = [np.ones(len(rows),bool)] if train is None else [np.asarray(train,bool),~np.asarray(train,bool)]
    for mask in masks:
        if not mask.any(): continue
        z = rows[mask]; check(name,p,s,z[r],z[v],y_lower=z[v]-z[e],y_upper=z[v]+z[e]); s += 1
    return s

def loss(rows):
    v = rows.holdout_chi2_per_row.to_numpy(float)
    return np.where(rows.holdout_admissible & np.isfinite(v) & (v>=0),v,np.inf)

chosen = ['DDO154','UGC11557','F568-3','NGC4559','UGC06614','NGC2955']
raw = {}; morphology = {}
for line in (root/'observations/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
    p = line.split()
    if len(p)==19:
        try: morphology[p[0]] = int(p[1])
        except ValueError: pass
for line in (root/'observations/data/MassModels_Lelli2016c.mrt').read_text().splitlines():
    p = line.split()
    if len(p)==10 and p[0] in morphology: raw.setdefault(p[0],[]).append(list(map(float,p[1:])))
raw = {k:np.asarray(v) for k,v in raw.items()}
assert len(raw)==175 and sum(map(len,raw.values()))==3391

for folder, prefix in [('results',''),('wide_bounds','legacy/wide_bounds/')]:
    fits=read(f'observations/{folder}/galaxy_fit_results.csv'); fits=fits[fits.ml_disk==.5]
    points=read(f'observations/{folder}/observed_points_and_predictions.csv'); curves=read(f'observations/{folder}/model_curves.csv')
    names=json.loads((root/f'observations/{folder}/fit_summary.json').read_text())['displayed']
    rank=points.groupby('galaxy').Vobs_kms.median().sort_values(kind='stable')
    assert names==rank.index[np.floor((np.arange(6)+.5)*len(rank)/6).astype(int)].tolist()
    name=prefix+'sparc_rotation_comparison'
    for p,g in enumerate(names):
        z=select(curves,galaxy=g); obs=select(points,galaxy=g)
        np.testing.assert_array_equal(obs[['r_kpc','Vobs_kms','e_Vobs_kms']],raw[g][:,1:4])
        for s,col in enumerate(['Vbar_kms','Vnfw_kms','Vcore_kms']):check(name,p,s,z.r_kpc,z[col])
        observed(name,p,3,obs,'r_kpc','Vobs_kms','e_Vobs_kms')
    name=prefix+'sparc_sample_diagnostics'
    for s,m in enumerate(['core','nfw']):ecdf(name,0,s,fits[m+'_reduced_chi2'])
    check(name,1,0,fits.nfw_reduced_chi2,fits.core_reduced_chi2)
    v=np.r_[fits.nfw_reduced_chi2,fits.core_reduced_chi2]; lim=[v.min()*.7,v.max()*1.4];check(name,1,1,lim,lim)

rows=read('fit_tiers/tier_predictions.csv'); dense=read('fit_tiers/dense_curves.csv')
for tier,tag in [('fixed_inputs','baseline'),('expanded_family','expanded')]:
    for mode in ['inner','full']:
        name=f'tier_{tag}_{"outer" if mode=="inner" else "full"}'
        for p,g in enumerate(chosen):
            ref=select(rows,galaxy=g,tier='fixed_inputs',model='core',fit_mode='inner').sort_values('radius_kpc')
            np.testing.assert_array_equal(ref[['radius_kpc','observed_kms','error_kms']],raw[g][:,1:4])
            s=0
            if tier=='expanded_family':
                z=select(dense,galaxy=g,tier='matched_nuisance',model='core',fit_mode=mode);check(name,p,s,z.radius_kpc,z.prediction_kms);s+=1
            for m in ['core','nfw']:
                z=select(dense,galaxy=g,tier=tier,model=m,fit_mode=mode);check(name,p,s,z.radius_kpc,z.prediction_kms);s+=1
                if m=='core':check(name,p,s,z.radius_kpc,np.sqrt(z.catalogue_baryon_v2));s+=1
            observed(name,p,s,ref,train=~ref.is_outer if mode=='inner' else None)
for p,mode in enumerate(['full','inner']):
    for s,(tier,model) in enumerate([('fixed_inputs','core'),('matched_nuisance','core'),('expanded_family','core'),('fixed_inputs','nfw'),('expanded_family','nfw')]):
        z=select(rows,tier=tier,model=model,fit_mode=mode)
        if mode=='inner':z=z[z.is_outer]
        v=z.assign(loss=((z.prediction_kms-z.observed_kms)/z.error_kms)**2).groupby('galaxy').loss.mean()
        ecdf('tier_population_comparison',p,s,v,True)

base=read('regime_baseline/results/guard_0.001_merged_fits.csv');free=read('population_family/guard_0.001/free_family_fits.csv');groups=[select(base,model='core'),free,select(base,model='nfw')]
for m,g in enumerate(groups):
    for h in [0,1]:
        z=g[g.holdout==bool(h)];ecdf('guarded_family_comparison',h,m,loss(z) if h else z.chi2_data_per_row)
for s,h in enumerate([False,True]):
    z=free[free.holdout==h];check('guarded_family_comparison',2,s,np.arange(3)+(.17 if h else -.17),[(z.eta_boundary==c).sum() for c in ['lower','none','upper']])
z=free[~free.holdout];check('guarded_family_comparison',3,0,z.original_potential_depth_ratio,z.potential_depth_ratio)
check('guarded_family_comparison',3,1,[1e-10,100],[1e-10,100]);check('guarded_family_comparison',3,2,[1e-10,100],[.001,.001])
nu=read('observations/uncertainty/penalized_fits.csv');configs=['fixed','stellar_0.10','geometry','all_0.05','all_0.10','all_0.20']
for s,m in enumerate(['core','nfw']):
    check('nuisance_and_holdout',0,s,np.arange(6)+(-.1 if m=='core' else .1),[select(nu,model=m,configuration=c,holdout=False).chi2_data_per_row.median() for c in configs])
    for j,c in enumerate(['fixed','all_0.10']):ecdf('nuisance_and_holdout',1,s*2+j,loss(select(nu,model=m,configuration=c,holdout=True)))
rad=read('observations/diagnostics/radial_summary.csv');d=select(read('observations/diagnostics/galaxy_diagnostics.csv'),run='ml_0.5')
for s,m in enumerate(['core','nfw']):
    z=select(rad,run='ml_0.5',binning='r_over_rmax',model=m).set_index('radial_bin').loc[['inner','middle','outer']]
    check('residual_physics_diagnostics',0,s,np.arange(3)+(-.08 if m=='core' else .08),z.mean_z_median,y_lower=z.mean_z_q25,y_upper=z.mean_z_q75)
    check('residual_physics_diagnostics',1,s,d.observed_outer_slope,d[m+'_outer_slope'])
v=np.r_[d.observed_outer_slope,d.core_outer_slope,d.nfw_outer_slope];v=v[np.isfinite(v)];lims=[v.min()-.05,v.max()+.05];check('residual_physics_diagnostics',1,2,lims,lims)
for s,b in enumerate([False,True]):
    z=d[(d.baryon_fraction_median>=.5)==b];check('residual_physics_diagnostics',2,s,z.r_max_over_Rdisk,z.log10_core_over_nfw_chi2)
    z=d[d.nfw_bound_hit==b];check('residual_physics_diagnostics',3,s,z.observed_outer_slope,z.nfw_r_max_over_peak)

st=read('significance/primary_summary.csv');paired=read('significance/paired_galaxy_losses.csv')
for i,r in st.iterrows():
    wins=r.right_wins if i<2 else r.left_wins;lo=r.win_ci95_lower;hi=r.win_ci95_upper
    if i<2:lo,hi=1-hi,1-lo
    check('paired_prediction_diagnostics',0,i,[wins/131],[2-i],x_lower=[lo],x_upper=[hi])
delta=np.sort(select(paired,analysis='guard_0.001',split='holdout',left_model='variable_core',right_model='fixed_core').delta_left_minus_right)
for s,b in enumerate([True,False]):
    ix=np.flatnonzero((delta<0)==b);check('paired_prediction_diagnostics',1,s,ix+1,delta[ix])

base=read('regime_baseline/results/guard_0.001_predictions.csv');free=read('population_family/guard_0.001/free_family_predictions.csv');models=[select(base,model='core'),free,select(base,model='nfw')]
ref=models[0][models[0].holdout];galaxies=sorted(ref.galaxy.unique());rank=ref.groupby('galaxy').Vobs_kms.median().sort_values(kind='stable')
assert chosen==rank.index[np.floor((np.arange(6)+.5)*len(rank)/6).astype(int)].tolist()
summary={};outer_count=0
for g in galaxies:
    r=select(ref,galaxy=g).sort_values('r_catalogue_kpc');mask=~r.training_row.to_numpy();outer_count+=mask.sum()
    np.testing.assert_array_equal(r[['r_catalogue_kpc','Vobs_kms','e_Vobs_kms']],raw[g][:,1:4])
    for m in range(3):
        for h in [False,True]:
            z=select(models[m],galaxy=g,holdout=h).sort_values('r_catalogue_kpc')
            np.testing.assert_array_equal(z[['r_catalogue_kpc','Vobs_kms','e_Vobs_kms']],r[['r_catalogue_kpc','Vobs_kms','e_Vobs_kms']])
            q=((z.Vpred_catalogue_kms-z.Vobs_kms)/z.Vobs_kms).to_numpy()[mask];summary[g,m,h]=(q.mean(),abs(q).mean())
assert outer_count==659
for name,displays in [('sparc_guarded_predictions',[(g,True) for g in chosen]),('sparc_fit_prediction_pair_a',[(g,h) for g in chosen[:3] for h in [False,True]]),('sparc_fit_prediction_pair_b',[(g,h) for g in chosen[3:] for h in [False,True]])]:
    for p,(g,h) in enumerate(displays):
        r=select(ref,galaxy=g).sort_values('r_catalogue_kpc')
        for m in range(3):
            z=select(models[m],galaxy=g,holdout=h).sort_values('r_catalogue_kpc');check(name,p,m,z.r_catalogue_kpc,z.Vpred_catalogue_kms)
        c=raw[g][:,4:7];v2=(c*abs(c))@np.array([1,.5,.7]);check(name,p,3,raw[g][:,1],np.sqrt(np.where(v2>=0,v2,np.nan)))
        observed(name,p,4,r,'r_catalogue_kpc','Vobs_kms','e_Vobs_kms',r.training_row if h else None)
for s,(lo,hi) in enumerate([(1,9),(10,11),(0,0)]):
    names=[g for g in galaxies if lo<=morphology[g]<=hi];check('outer_prediction_mass_budget',0,s,[100*summary[g,0,True][0] for g in names],[100*summary[g,0,False][0] for g in names])
for m in range(3):
    for j,h in enumerate([True,False]):ecdf('outer_prediction_mass_budget',1,m*2+j,[100*summary[g,m,h][1] for g in galaxies],True)
tail=read('morphology/tail_results/tail_galaxy_diagnostics.csv')
for m,model in enumerate(['fixed_core','variable_core']):
    z=select(tail,model=model);ecdf('outer_prediction_mass_budget',2,m,100*z.mass_fraction_last,True);ecdf('outer_prediction_mass_budget',3,m,z.mass_fraction_last*(1+z.last_required_additional_enclosed_over_fitted),True)

morph=select(read('morphology/diagnostics/galaxy_model_diagnostics.csv'),epsilon=.001,analysis='outer_holdout');morphs=['spiral_T1_9','irregular_T10_11','S0_T0']
for name in ['morphology_diagnostics','legacy/morphology/morphology_diagnostics']:
    for m,model in enumerate(['fixed_core','variable_core','nfw']):
        for g in range(2):
            z=select(morph,model=model,morphology=morphs[g]);ecdf(name,0,m*2+g,100*z.fractional_residual_mean)
            check(name,3,m*2+g,[m+(g-.5)*.18],[z.slope_residual.median()],y_lower=[z.slope_residual.quantile(.25)],y_upper=[z.slope_residual.quantile(.75)])
    for g in range(3):
        z=select(morph,model='fixed_core',morphology=morphs[g]);check(name,1,g,z.coverage,100*z.fractional_residual_mean)
        if g<2:check(name,2,g,z.median_speed,z.coverage)
g=select(read('outer_mechanism/growth/growth_galaxy_diagnostics.csv'),guard=.001,anchor_offset=0);fixed=select(g,model='fixed_core');variable=select(g,model='variable_core')
for s,v in enumerate([fixed.core_mass_growth,variable.core_mass_growth,fixed.nfw_mass_growth]):ecdf('outer_growth_mechanism',0,s,v)
for s,v in enumerate([fixed.core_original_loss,variable.core_original_loss,fixed.nfw_original_loss,fixed.core_with_nfw_growth_loss,variable.core_with_nfw_growth_loss]):ecdf('outer_growth_mechanism',1,s,v)

bar=read('baryons_only/baseline_predictions.csv');saved=read('baryons_only/galaxy_scores.csv');ids={'bar_fixed':('fixed',''),'bar_profiled':('profiled',''),'A_core':('fixed_inputs','core'),'A_nfw':('fixed_inputs','nfw'),'B_core':('matched_nuisance','core'),'C_core':('expanded_family','core'),'B_nfw':('matched_nuisance','nfw')}
for p,(matched,mode) in enumerate([(False,'full'),(False,'inner'),(True,'full'),(True,'inner')]):
    for s,model in enumerate(['bar_profiled','B_core','C_core','B_nfw'] if matched else ['bar_fixed','A_core','A_nfw']):
        tier,kind=ids[model];z=select(bar,baseline=tier,fit_mode=mode) if not kind else select(rows,tier=tier,model=kind,fit_mode=mode);vals=[]
        for galaxy,part in z.groupby('galaxy'):
            q=part if mode=='full' else part[part.is_outer]
            fail=not np.isfinite(q.prediction_kms).all() or model=='bar_profiled' and part.pipeline_failure.any()
            v=np.inf if fail else np.mean(((q.prediction_kms-q.observed_kms)/q.error_kms)**2)
            old=select(saved,model_id=model,fit_mode=mode,region='all' if mode=='full' else 'outer',galaxy=galaxy).loss.item()
            np.testing.assert_allclose(v,old,rtol=3e-11,atol=2e-10);vals.append(v)
        ecdf('baryons_only_population',p,s,vals,True)
for p,g in enumerate(chosen):
    for s,model in enumerate(['bar_profiled','B_core','B_nfw','C_core']):
        tier,kind=ids[model];z=select(bar,baseline=tier,fit_mode='inner',galaxy=g) if not kind else select(rows,tier=tier,model=kind,fit_mode='inner',galaxy=g)
        v=z.prediction_kms.to_numpy();v=np.full(len(v),np.nan) if model=='bar_profiled' and z.pipeline_failure.any() else v
        check('baryons_only_outer_examples',p,s,z.radius_kpc,v)
    z=select(rows,tier='fixed_inputs',model='core',fit_mode='inner',galaxy=g);observed('baryons_only_outer_examples',p,4,z,train=~z.is_outer)

for folder in ['results','wide_results']:
    name=f'legacy/population_family/{folder}/family_population_comparison';free=read(f'population_family/{folder}/free_family_fits.csv');groups=[select(nu,model='core',configuration='all_0.10'),free,select(nu,model='nfw',configuration='all_0.10')]
    for m,g in enumerate(groups):
        for h in [0,1]:
            z=g[g.holdout==bool(h)];ecdf(name,h,m,loss(z) if h else z.chi2_data_per_row)
    z=free[~free.holdout];check(name,2,0,np.maximum(1e-10,z.eta),z.physical_length_kpc*z.physical_amplitude_kms)
    z=free[free.holdout];nfw=groups[2][groups[2].holdout].set_index('galaxy').loc[z.galaxy];check(name,3,0,loss(nfw),loss(z));check(name,3,1,[.001,1e5],[.001,1e5])
pilot=read('coupled_pilot/predictions.csv')
for p,(g,mode) in enumerate([(g,mode) for g in ['F583-1','UGC00731'] for mode in ['inner','full']]):
    z=select(pilot,galaxy=g,split=mode)
    for s,v in enumerate([z.coupled_kms,z.nfw_kms,np.sqrt(z.baryon_v2)]):check('coupled_galaxy_pilot',p,s,z.r_kpc,v)
    observed('coupled_galaxy_pilot',p,3,z,'r_kpc','observed_kms','error_kms',z.training)

gallery_axes=json.loads((out/'baryons_score_verification.json').read_text())['gallery_axes']
for g in chosen:
    ref=select(rows,tier='fixed_inputs',model='core',fit_mode='inner',galaxy=g)
    halo=select(rows,galaxy=g);curve=select(dense,galaxy=g);b=select(bar,baseline='profiled',fit_mode='inner',galaxy=g)
    old=np.nanmax(np.r_[halo.prediction_kms,curve.prediction_kms,ref.observed_kms+ref.error_kms,np.sqrt(np.where(ref.catalogue_baryon_v2>=0,ref.catalogue_baryon_v2,np.nan))])*1.08
    ymax=max(old,np.nanmax(b.prediction_kms)*1.08);minimum=np.min(ref.observed_kms-ref.error_kms);ymin=min(0,minimum-.025*ymax) if minimum<0 else 0
    np.testing.assert_allclose(gallery_axes[g],[0,ref.radius_kpc.iloc[-1]*1.025,ymin,ymax],rtol=2e-14,atol=1e-12)
pilot_axes=json.loads((out/'pilot_shared_axes.json').read_text())
for g in ['F583-1','UGC00731']:
    ref=select(pilot,galaxy=g,split='inner');span=ref.r_kpc.to_numpy();pad=.03*(span[-1]-span[0]);v=np.r_[ref.coupled_kms,ref.nfw_kms,np.sqrt(ref.baryon_v2),ref.observed_kms-ref.error_kms,ref.observed_kms+ref.error_kms]
    np.testing.assert_allclose(pilot_axes[g],[max(0,span[0]-pad),span[-1]+pad,0,v.max()+.05*(v.max()-v.min())],rtol=2e-14,atol=1e-12)

remaining={}
for name,frame in frames.items():
    keys=set(map(tuple,frame[['panel','series']].drop_duplicates().to_numpy()))
    extra=keys-checked[name]
    if extra:remaining[name]=[list(map(int,v)) for v in sorted(extra)]
if remaining:raise AssertionError(('Unchecked plot series',remaining))
report=dict(passed=True,shared_pilot_axes_verified=True,preserved_baryon_gallery_axes_verified=True,figures=len(frames),series=sum(map(len,checked.values())),numeric_values=numerical_values,max_relative_error_scale_at_least_one=float(max_relative_error),tolerance=dict(relative=3e-11,absolute=2e-10),displayed=chosen,raw_catalogue_galaxies=len(raw),raw_catalogue_radii=sum(map(len,raw.values())),outer_galaxies=len(galaxies),outer_radii=int(outer_count),source_inputs_sha256=inputs,scope='Independent NumPy/Pandas verification of every ObservationPlots numerical series, including error bars, exact display ranks, signed raw SPARC baryonic components, matched outer rows, ECDF failure denominators and reflected win intervals. SVG styling is intentionally independent; no pixel-identity claim.')
(a.output/'observation_plot_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='source_inputs_sha256'},indent=2))
