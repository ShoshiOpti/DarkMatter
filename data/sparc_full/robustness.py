"""Additional declared robustness checks; never select the primary by these outcomes."""
from pathlib import Path
import os,sys,json,csv
os.environ['OPENBLAS_NUM_THREADS']='1'
import test_sparc as t
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit,logit
HERE=t.HERE

def fit_nuisance(d,p,st=None):
 r=d['r'];mask=~d['outer'];sigi=max(d['incerr'],.1);sigD=max(d['Derr']/d['D'],.001);sigS=.1*np.log(10)
 lo=np.log(r.min()/100);hi=np.log(r.max()*100)
 def pred(v,full=False):
  if p is None:ld,inc,ls=v
  else:ll,eta,ld,inc,ls=v
  dr=np.exp(ld);ti=np.sin(np.deg2rad(inc))/np.sin(np.deg2rad(d['inc']));sf=np.exp(ls)
  B=dr*ti*ti*(d['gas']+sf*d['star']);pen=0.
  if p is None:a=0.;Lcat=0.;pr=np.sqrt(np.maximum(B,0))
  else:
   Lcat=np.exp(ll);x=r/Lcat;H=p.mass(x)/x;lower=max(0.,float(np.max(-B/H)))*(1+1e-12)
   upper=ti*ti*.001*t.C2/p.depth
   pen=max(lower/upper-1,0)*1000
   a=lower+max(upper-lower,0)*expit(eta);pr=np.sqrt(np.maximum(B+a*H,0))
  prior=np.array([(dr-1)/sigD,(inc-d['inc'])/sigi,ls/sigS,pen])
  if full:return pr,{'distance_ratio':float(dr),'inclination_deg':float(inc),'stellar_factor':float(sf),'halo_L_kpc':float(Lcat*dr),'halo_mass_Msun':float((a/ti**2)*Lcat*dr/t.G),'guard_penalty':float(pen)}
  return np.r_[(pr-d['y'])[mask]/d['err'][mask],prior]
 if p is None:
  lb=[np.log(.05),5.,np.log(.1)];ub=[np.log(5.),90.,np.log(4.)];starts=[[0.,d['inc'],0.]]
 else:
  x=r/st['L_kpc'];hh=p.mass(x)/x;lower=max(0.,float(np.max(-d['B']/hh)))*(1+1e-12);upper=.001*t.C2/p.depth
  frac=np.clip((st['amplitude_GM_over_L']-lower)/max(upper-lower,1),1e-12,1-1e-12);e=logit(frac)
  lb=[lo,-35.,np.log(.05),5.,np.log(.1)];ub=[hi,35.,np.log(5.),90.,np.log(4.)]
  starts=[[np.clip(np.log(st['L_kpc'])+dl,lo+1e-8,hi-1e-8),e,0.,d['inc'],0.] for dl in [-.5,0.,.5]]
 fits=[least_squares(pred,np.clip(v,np.array(lb)+1e-9,np.array(ub)-1e-9),bounds=(lb,ub),xtol=2e-9,ftol=2e-9,gtol=2e-9,max_nfev=900) for v in starts]
 f=min(fits,key=lambda z:sum(z.fun*z.fun));pr,rec=pred(f.x,True);rec.update({'objective':float(sum(f.fun*f.fun)),'success':bool(f.success),'optimality':float(f.optimality),'nfev':int(f.nfev)})
 return pr,rec

def main():
 data=t.load();profiles=t.profiles();p=profiles['new_analytic'];rows=list(csv.DictReader((HERE/'fit_states.csv').open()));state={r['galaxy']:{k:float(r[k]) for k in ['L_kpc','amplitude_GM_over_L']} for r in rows if r['phase']=='inner' and r['model']=='new_analytic'}
 scores=[];params=[];predrows=[]
 for i,(n,d) in enumerate(data.items()):
  for name,prof in [('baryons',None),('new_analytic',p)]:
   pr,rec=fit_nuisance(d,prof,state[n]);res=(pr-d['y'])[d['outer']];err=d['err'][d['outer']]
   params.append({'galaxy':n,'model':name,**rec});scores.append({'galaxy':n,'model':name,'loss':float(np.mean((res/err)**2)),'mae':float(np.mean(abs(res)))})
   for j in range(len(pr)):predrows.append({'galaxy':n,'model':name,'index':j,'predicted':float(pr[j])})
  if i%25==0:print('nuisance sensitivity',i+1,'/175',flush=True)
 t.csvout('full_nuisance_states.csv',params);t.csvout('full_nuisance_scores.csv',scores);t.csvout('full_nuisance_predictions.csv',predrows)
 lookup={(r['galaxy'],r['model']):r for r in scores};record={'scope':'Additional inner-only profile fits; catalogue Gaussian D and inclination priors, shared stellar 0.1 dex prior. Conditional sensitivity, not a marginal likelihood.', 'failed_optimizer_count':sum(not s['success'] for s in params),'max_guard_penalty':max(s['guard_penalty'] for s in params)}
 for sample,names in [('all175',list(data)),('quality131',[n for n,d in data.items() if d['quality']])]:
  diff=[lookup[n,'new_analytic']['loss']-lookup[n,'baryons']['loss'] for n in names]
  record[sample]={'baryon_mean_loss':float(np.mean([lookup[n,'baryons']['loss'] for n in names])), 'new_mean_loss':float(np.mean([lookup[n,'new_analytic']['loss'] for n in names])),**t.inference(diff,[data[n]['group'] for n in names])}
 t.savejson('nuisance_sensitivity.json',record);print(json.dumps(record,indent=2),flush=True)
 # Tighter weak-field guard: the radius/mass extrapolation cannot hide bound dependence.
 bound=[]
 original=t.C2
 for guard in [1e-4,1e-5]:
  t.C2=original*guard/.001
  for n,d in data.items():
   st,pr=t.fit(d,p,~d['outer']);rr=(pr-d['y'])[d['outer']]/d['err'][d['outer']]
   bound.append({'galaxy':n,'guard':guard,'loss':float(np.mean(rr*rr)),'M_total_Msun':st['M_total_Msun'],'L_kpc':st['L_kpc'],'bound_hit':st['weak_field_bound']})
  print('guard',guard,'finished',flush=True)
 t.C2=original;t.csvout('guard_sensitivity.csv',bound)
 t.savejson('guard_sensitivity.json',{str(g):{'mean_outer_loss_all175':float(np.mean([r['loss'] for r in bound if r['guard']==g])), 'mean_outer_loss_quality131':float(np.mean([r['loss'] for r in bound if r['guard']==g and data[r['galaxy']]['quality']])), 'bound_hits':sum(r['bound_hit'] for r in bound if r['guard']==g)} for g in [1e-4,1e-5]})

if __name__=='__main__':main()
