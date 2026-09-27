"""Development-only empirical refitting, with outer outcomes excluded from the fit process.

Original fit algorithms execute from isolated, provenance-recorded copies. The
perturbation ensemble measures sensitivity to an explicitly stipulated repeated
measurement model; it is neither posterior sampling nor calibrated uncertainty.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import contextlib

sys.dont_write_bytecode = True

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def save(path, value):
    def safe(item):
        if isinstance(item, dict): return {k:safe(v) for k,v in item.items()}
        if isinstance(item, (list,tuple)): return [safe(v) for v in item]
        if hasattr(item, 'item'): return safe(item.item())
        if isinstance(item, float) and not __import__('math').isfinite(item): return None
        return item
    Path(path).write_text(json.dumps(safe(value), indent=2, allow_nan=False)+'\n', encoding='utf-8')

def write_csv(path, rows):
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    p=read(args.protocol);lock=read(args.protocol.with_name('protocol.lock.json'))
    if sha(args.protocol)!=lock['protocol_sha256']:raise ValueError('Changed locked protocol')
    study=args.protocol.resolve().parent;out=args.out.resolve();reference=args.reference.resolve()
    if out.is_relative_to(study) or out.is_relative_to(reference):raise ValueError('Output must be outside study/reference')
    out.mkdir(parents=True,exist_ok=True)
    if (out/'fit_artifact.json').exists() or (out/'backend').exists():raise FileExistsError('Pilot output already used')
    # Check every frozen input hash without opening evaluation outcomes for fitting.
    # The evaluation file is checked by the native caller; it is deliberately never read here.
    permitted={p['files'][k] for k in ['predictors','training','metadata','split']}
    for name in permitted:
        if sha(study/name)!=lock['files'][name]:raise ValueError('Changed locked fit input: '+name)
    if p['study_status']!='retrospective_development':raise ValueError('This pilot backend is development-only; independent samples require a separately reviewed external adapter')
    if p['pilot']['backend']!='empirical_original':raise ValueError('Physical refit adapters are external contracts, not implemented by this empirical pilot')
    if p['shared_parameter_models'] or p['shared_training_galaxies'] or p['ensemble_contract']['kind']!='nuisance_sensitivity':raise ValueError('Empirical pilot supports only nonshared nuisance sensitivity')
    if p['models']!=['core','nfw','profiled_baryons','envelope']:raise ValueError('Empirical backend requires its declared four-model family')
    rules=p['fitting_rules']
    for field,expected in [('stellar_centres',[.5,.7]),('stellar_width_dex',.1),('pull_bound',6),('distance_ratio_floor',.001),('core_depth_guard',.001),('envelope_taper',20),('optimizer_budget','original_multistart_and_tolerances')]:
        if rules[field]!=expected:raise ValueError('Unsupported empirical backend rule: '+field)
    import numpy as np
    import pandas as pd
    import scipy
    chosen=p['evaluation_galaxies'];predictors=[r for r in table(study/p['files']['predictors']) if r['galaxy'] in chosen]
    training={(r['galaxy'],int(r['index'])):r for r in table(study/p['files']['training'])}
    metadata=read(study/p['files']['metadata'])
    replicas=int(p['pilot']['perturbation_replicates'])
    if not 1<=replicas<=1000:raise ValueError('Perturbation replicate count must be 1..1000')
    rng=np.random.default_rng(int(p['pilot']['seed']))
    originals=reference/'reproducibility/source/numerics';stage=out/'backend/numerics'
    paths=['holdout_scoring.py','observations/fit_sparc.py','observations/uncertainty_sensitivity.py','core/core_study.py','regime_baseline/guard_fixed_core.py','baryons_only/fit_baselines.py','conditional_envelope/fit_envelopes.py']
    original_hashes={};adaptations=[]
    for rel in paths:
        source=originals/rel;target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);original_hashes[rel]=sha(source)
    def patch(rel,old,new):
        path=stage/rel;text=path.read_text(encoding='utf-8-sig')
        if old not in text:raise ValueError('Original adapter anchor changed: '+rel+' '+old)
        text=text.replace(old,new);path.write_text(text,encoding='utf-8');adaptations.append(dict(file=rel,old=old,new=new))
    # Replicated catalogue estimates change prior centres, but predictions remain
    # in the original catalogue coordinate. Fixed reference D/i denominators matter.
    ufile='observations/uncertainty_sensitivity.py';gfile='regime_baseline/guard_fixed_core.py';bfile='baryons_only/fit_baselines.py'
    patch(ufile,"d=1+(m['e_D']/m['D'])*dct.get('distance_pull',0.)","d=(m['D']+m['e_D']*dct.get('distance_pull',0.))/m['catalogue_D']")
    patch(ufile,"np.sin(np.deg2rad(m['inc']))","np.sin(np.deg2rad(m['catalogue_inc']))")
    patch(ufile,"(G.get('distance_floor',.001)-1)/(m['e_D']/m['D'])","(G.get('distance_floor',.001)*m['catalogue_D']-m['D'])/m['e_D']")
    patch(ufile,"distance_Mpc=d*m['D']","distance_Mpc=d*m['catalogue_D']")
    patch(gfile,"distance_Mpc=d*self.meta['D']","distance_Mpc=d*self.meta['catalogue_D']")
    patch(gfile,"ratio=1+m['e_D']/m['D']*d['distance_pull']","ratio=(m['D']+m['e_D']*d['distance_pull'])/m['catalogue_D']")
    patch(gfile,"np.sin(np.deg2rad(m['inc']))","np.sin(np.deg2rad(m['catalogue_inc']))")
    patch(gfile,"(.001-1)/(m['e_D']/m['D'])","(.001*m['catalogue_D']-m['D'])/m['e_D']")
    patch(bfile,"self.fd=meta['e_D']/meta['D'];","self.fd=meta['e_D']/meta['catalogue_D'];self.d0=meta['D']/meta['catalogue_D'];")
    patch(bfile,"np.sin(np.deg2rad(self.i0))","np.sin(np.deg2rad(meta['catalogue_inc']))")
    patch(bfile,"(.001-1)/self.fd","(.001-self.d0)/self.fd")
    patch(bfile,"d=1+self.fd*x[-2]","d=self.d0+self.fd*x[-2]")
    sys.path.insert(0,str(stage));u=load(stage/ufile,'research_u');guard=load(stage/gfile,'research_guard');baryon=load(stage/bfile,'research_baryon');fitter=load(stage/'observations/fit_sparc.py','research_fitter');core=load(stage/'core/core_study.py','research_core')
    solution=core.solve_core(24.,2e-11);mass=float(24**2*solution.sol(24)[3]);phi0=abs(float(solution.sol(0)[2]))
    envelope_original=(stage/'conditional_envelope/fit_envelopes.py').read_text(encoding='utf-8-sig')
    # Use the unchanged original optimization body; replace input parsing and omit
    # post-fit historical 131-galaxy summaries/figures. No evaluation outcomes enter.
    start=envelope_original.index('meta={}');end=envelope_original.index('def shape(')
    envelope_script=envelope_original[:start]+'meta=json.loads((ROOT/"inputs/metadata.json").read_text())\n\n'+envelope_original[end:]
    envelope_script=envelope_script[:envelope_script.index('summary=[]')]
    adaptations.append(dict(file='conditional_envelope/fit_envelopes.py',change='Metadata supplied by locked study, sample-size assertion removed for declared subset, stop after optimizer/prediction CSV outputs; optimization equations unchanged'))
    save(out/'adapter_provenance.json',dict(original_sources=original_hashes,changes=adaptations,adapted_source_sha256={r:sha(stage/r) for r in paths},adapter_sha256=sha(__file__),transformed_envelope_sha256=hashlib.sha256(envelope_script.encode('utf-8')).hexdigest()))
    results=[];fit_records=[];draw_records=[]
    for rep in range(replicas+1):
        print('Empirical refit replicate',rep,'of',replicas,flush=True)
        folder=out/('replicate_'+str(rep));folder.mkdir()
        frames=[];model_fits={}
        for name in chosen:
            rows=sorted((r for r in predictors if r['galaxy']==name),key=lambda r:int(r['index']))
            n=len(rows);nt=sum(r['partition']=='inner' for r in rows)
            if nt!=n-max(2,int(np.ceil(.2*n))):raise ValueError('Original empirical adapter requires the declared outer-20% suffix: '+name)
            meta=dict(metadata[name]);meta.setdefault('catalogue_D',meta['D']);meta.setdefault('catalogue_inc',meta['inc'])
            if rep:
                for _ in range(10000):
                    d=float(rng.normal(metadata[name]['D'],meta['e_D']));inc=float(rng.normal(metadata[name]['inc'],meta['e_inc']))
                    if d>0 and .01<inc<90:break
                else:raise RuntimeError('Physical-domain Gaussian nuisance draw failed')
                meta.update(D=d,inc=inc)
            noise=rng.standard_normal(nt) if rep else np.zeros(nt)
            a=np.zeros((n,9));a[:,0]=meta['catalogue_D']
            for i,r in enumerate(rows):
                a[i,1]=float(r['radius_kpc']);a[i,3]=float(r['sigma_kms']);a[i,4:7]=[float(r[k]) for k in ['gas_kms','disk_kms','bulge_kms']]
                a[i,2]=float(training[(name,i)]['observed_kms'])+noise[i]*a[i,3] if i<nt else np.nan
            # Equal-galaxy objective = mean inner z^2 + same squared prior pulls.
            # Outer errors and evaluation data are never modified.
            if rules['training_weight']=='equal_galaxy':a[:nt,3]*=np.sqrt(nt)
            draw_records.append(dict(replicate=rep,galaxy=name,prior_centre_D=meta['D'],prior_centre_inc=meta['inc'],velocity_standard_normal=noise.tolist()))
            u.G.clear();u.G.update(data={name:a},meta={name:meta},sol=solution,mass=mass,fitter=fitter,warm={},pull_limit=6.,distance_floor=.001)
            for model in ['core','nfw']:
                try:
                    row,_=u.solve((name,model,'all_0.10',True))
                    factor=np.sin(np.deg2rad(row['inclination_deg']))/np.sin(np.deg2rad(meta['catalogue_inc']))
                    depth=phi0*row['amplitude_catalogue_kms']**2/(factor**2*299792.458**2)
                    guard_records=[]
                    if model=='core' and depth>.001:
                        row,guard_records=guard.Problem(u,row,.001,phi0).refit()
                        factor=np.sin(np.deg2rad(row['inclination_deg']))/np.sin(np.deg2rad(meta['catalogue_inc']))
                        depth=phi0*row['amplitude_catalogue_kms']**2/(factor**2*299792.458**2)
                    bar=row['distance_ratio']*factor**2*((a[:,4:7]*abs(a[:,4:7]))@np.array([1.,row['ml_disk'],row['ml_bulge']]))
                    x=a[:,1]/row['length_catalogue_kpc'];shape=u.core_shape(x) if model=='core' else fitter.nfw_shape(x);halo=row['amplitude_catalogue_kms']**2*shape
                    squared=bar+halo;pred=np.full(n,np.nan);valid=squared>=0;pred[valid]=np.sqrt(squared[valid]);failed=not row['success'] or not np.isfinite(pred).all() or (model=='core' and depth>.001*(1+1e-7))
                    fit_records.append(dict(replicate=rep,galaxy=name,model=model,failed=failed,fit=row,guard_records=guard_records))
                    model_fits[(name,model)]=(pred,failed)
                    if model=='core' and not failed:
                        for i,r in enumerate(rows):frames.append(dict(galaxy=name,r=a[i,1],obs=a[i,2],error=a[i,3],pred=pred[i],baryon_v2=bar[i],halo_v2=halo[i],halo_mass_fraction=x[i]*shape[i]/mass,distance_ratio=row['distance_ratio'],inclination_factor=factor,guard_depth_ratio=depth,is_training_fit=i<nt,is_outer=i>=nt,r_train_edge=a[nt-1,1],current_display=False))
                except Exception as error:
                    model_fits[(name,model)]=(np.full(n,np.nan),True);fit_records.append(dict(replicate=rep,galaxy=name,model=model,failed=True,error=str(error)))
            try:
                problem=baryon.Problem(a,meta,True);state,audit=problem.solve();pred=np.full(n,np.nan) if state is None else problem.prediction(state)[0];failed=state is None or not np.isfinite(pred).all();model_fits[(name,'profiled_baryons')]=(pred,failed);fit_records.append(dict(replicate=rep,galaxy=name,model='profiled_baryons',failed=failed,audit=audit))
            except Exception as error:
                model_fits[(name,'profiled_baryons')]=(np.full(n,np.nan),True);fit_records.append(dict(replicate=rep,galaxy=name,model='profiled_baryons',failed=True,error=str(error)))
        # Run the original envelope grid separately for each galaxy: one failed
        # auxiliary configuration must not erase other galaxies or replicates.
        for galaxy_index,name in enumerate(chosen):
            galaxy_frames=[row for row in frames if row['galaxy']==name]
            if not galaxy_frames:
                fit_records.append(dict(replicate=rep,galaxy=name,model='envelope',failed=True,error='No successful guarded-core input'))
                continue
            ep=folder/('envelope_'+str(galaxy_index));(ep/'inputs').mkdir(parents=True)
            pd.DataFrame(galaxy_frames).to_csv(ep/'inputs/fixed_core_frozen_profile_inputs.csv',index=False)
            pd.DataFrame(galaxy_frames).to_csv(ep/'inputs/outer_predictions.csv',index=False)
            save(ep/'inputs/metadata.json',metadata);script=ep/'fit_envelopes.py';script.write_text(envelope_script,encoding='utf-8')
            try:
                with (ep/'fit.log').open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
                    exec(compile(envelope_script,str(script),'exec'),{'__file__':str(script),'__name__':'research_envelope'})
                fitted=pd.read_csv(ep/'results/galaxy_fits.csv');predicted=pd.read_csv(ep/'results/predictions.csv')
                rec=fitted[(fitted.galaxy==name)&(fitted.model=='envelope_a1_q20')]
                if len(rec)!=1:raise ValueError('Primary envelope fit missing or duplicated')
                fit=rec.iloc[0];pred=predicted[(predicted.galaxy==name)&(predicted.model=='envelope_a1_q20')].sort_values('r').pred.to_numpy()
                failed=not bool(fit.optimizer_success) or not np.isfinite(pred).all()
                model_fits[(name,'envelope')]=(pred,failed)
                fit_records.append(dict(replicate=rep,galaxy=name,model='envelope',failed=failed,fit=fit.to_dict()))
            except Exception as error:
                # Conservative policy: any failure in the retained grid makes this
                # galaxy's envelope slot fail. Never publish a partially completed grid.
                fit_records.append(dict(replicate=rep,galaxy=name,model='envelope',failed=True,error=str(error)))
        for name in chosen:
            rows=sorted((r for r in predictors if r['galaxy']==name),key=lambda r:int(r['index']))
            for model in p['models']:
                pred,failed=model_fits.get((name,model),(np.full(len(rows),np.nan),True))
                for i,r in enumerate(rows):
                    if r['partition']=='outer':results.append(dict(replicate=rep,galaxy=name,index=i,model=model,predicted_kms='' if failed else float(pred[i]),fit_failed=failed))
    write_csv(out/'predictions.csv',results);save(out/'fit_records.json',fit_records);save(out/'perturbations.json',draw_records)
    if any(sha(originals/r)!=h for r,h in original_hashes.items()):raise RuntimeError('Original scientific sources changed')
    save(out/'fit_artifact.json',dict(protocol_sha256=sha(args.protocol),predictions_sha256=sha(out/'predictions.csv'),fitting_rules_json=lock['fitting_rules_json'],read_evaluation_outcomes=False,training_keys=sorted(f'{g}:{i}' for g,i in training if g in chosen),shared_training_keys=[],has_shared_parameters=False,replicates=replicas+1,ensemble_kind='nuisance_sensitivity',backend='empirical_original',source_hashes=original_hashes,adapter_sha256=sha(__file__),scope='Baseline plus paired measurement/prior-centre sensitivity refits; empirical independent galaxy parameters, no shared physical fit, no posterior draws or calibrated interval',software=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__)))
    print('PASS: empirical pilot refitted and produced blind outer predictions; evaluation delegated to native process',flush=True)

if __name__=='__main__':main()
