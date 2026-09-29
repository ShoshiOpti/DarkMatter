"""Full SPARC full-profile tests. Conditional templates, not shared-microphysics fits."""
from pathlib import Path
import os,sys,json,csv,hashlib,time
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
if not (ROOT/'Cubic_Scalar_GR_Intersection.tex').exists():ROOT=HERE
sys.path.insert(0,str(ROOT/'output/scientific_runtime'))
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq,minimize_scalar,minimize
from scipy.integrate import quad
from scipy.special import exp1
from scipy.stats import norm,binomtest
G=4.300917270e-6;C2=299792.458**2

def savejson(name,x): (HERE/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csvout(name,rows):
 with (HERE/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def load():
 meta={}
 for l in (HERE/'inputs/SPARC_Lelli2016c.mrt').read_text().splitlines():
  a=l.split()
  if len(a)!=19:continue
  try:meta[a[0]]={'D':float(a[2]),'Derr':float(a[3]),'method':int(a[4]),'inc':float(a[5]),'incerr':float(a[6]),'Rd':float(a[11]),'Q':int(a[17])}
  except ValueError:continue
 data={n:[] for n in meta}
 for l in (HERE/'inputs/MassModels_Lelli2016c.mrt').read_text().splitlines():
  a=l.split()
  if len(a)==10 and a[0] in meta:data[a[0]].append(list(map(float,a[1:])))
 out={}
 for n in sorted(meta):
  a=np.array(data[n]);a=a[np.argsort(a[:,1])];m=meta[n];N=len(a);no=max(2,int(np.ceil(.2*N)))
  gas=a[:,4]*abs(a[:,4]);star=.5*a[:,5]*abs(a[:,5])+.7*a[:,6]*abs(a[:,6]);B=gas+star
  assert np.all(a[:,1]>0) and np.all(a[:,3]>0)
  out[n]={'r':a[:,1],'y':a[:,2],'err':a[:,3],'B':B,'gas':gas,'star':star,'outer':np.arange(N)>=N-no,
   'quality':m['Q']<3 and m['inc']>=30 and N>=8,'group':'distance:UrsaMajor' if m['method']==4 else n,**m}
 assert len(out)==175 and sum(len(d['r']) for d in out.values())==3391
 assert sum(d['quality'] for d in out.values())==131
 return out

class Density:
 """Positive PCHIP volume density, analytic segment integrals, nonzero exterior tail."""
 def __init__(self,x,radial,k,join=30.):
  keep=x<join;xx=x[keep];rho=radial[keep]/xx**2
  # Regular central extrapolation in r^2.
  rho0=max(0.,(rho[0]*xx[1]**2-rho[1]*xx[0]**2)/(xx[1]**2-xx[0]**2))
  full=PchipInterpolator(x,radial/x**2)
  self.x=np.r_[0,xx,join];self.rho=np.r_[rho0,rho,max(0,float(full(join)))];self.k=k;self.join=join
  p=PchipInterpolator(self.x,self.rho);self.pol=p.c[::-1].T
  self.masscoef=[];self.depthcoef=[]
  for l,co in zip(self.x[:-1],self.pol):
   a=np.polynomial.polynomial.polymul(co,[l*l,2*l,1]);self.masscoef.append(np.r_[0,a/np.arange(1,len(a)+1)])
   a=np.polynomial.polynomial.polymul(co,[l,1]);self.depthcoef.append(np.r_[0,a/np.arange(1,len(a)+1)])
  self.masscoef=np.array(self.masscoef);self.depthcoef=np.array(self.depthcoef)
  width=np.diff(self.x)
  def ev(co,h):return np.sum(co*h[:,None]**np.arange(co.shape[1])[None,:],axis=1)
  seg=ev(self.masscoef,width);self.prefix=np.r_[0,np.cumsum(seg)]
  self.tailmass=self.rho[-1]*join*join/(2*k);self.total=self.prefix[-1]+self.tailmass
  self.depth=(sum(ev(self.depthcoef,width))+self.rho[-1]*join**2*np.exp(2*k*join)*exp1(2*k*join))/self.total
 def mass(self,z):
  z=np.asarray(z);flat=z.ravel();j=np.searchsorted(self.x,flat,side='right')-1;j=np.clip(j,0,len(self.pol)-1)
  h=np.maximum(flat-self.x[j],0);co=self.masscoef[j];v=self.prefix[j]+np.sum(co*h[:,None]**np.arange(7)[None,:],axis=1)
  outside=flat>=self.join;v[outside]=self.prefix[-1]+self.tailmass*(-np.expm1(-2*self.k*(flat[outside]-self.join)))
  return (v/self.total).reshape(z.shape)

def f0(z):
 z=np.asarray(z);v=2/np.pi*(np.arctan(z)-z/(1+z*z));small=z<1e-3;t=z[small]
 v[small]=2/np.pi*sum((-1)**(j+1)*2*j/(2*j+1)*t**(2*j+1) for j in range(1,7));return v

def f2(z,q):return f0(z)+8*q/(3*np.pi)*np.asarray(z)**3/(1+np.asarray(z)**2)**3

class Mix:
 def __init__(self,core,pars,exact=False):
  self.core=core;self.sc=pars['core_dilation'];self.f=pars['envelope_fraction'];a=pars['a_over_L'];b=pars['b_over_L'];self.s=(a+b)/2;self.e=(b-a)/(b+a);self.q=self.e**2;self.exact=exact
  self.nodes,self.weights=np.polynomial.legendre.leggauss(8)
  dep=2*np.arctanh(self.e)/(np.pi*self.s*self.e) if exact else 2/(np.pi*self.s)*(1+self.q/3)
  self.depth=(1-self.f)*core.depth/self.sc+self.f*dep
 def mass(self,z):
  z=np.atleast_1d(z);x=z/self.s
  if self.exact:env=np.sum(self.weights[:,None]*f0((x[None,:]/(1+self.e*self.nodes[:,None])).ravel()).reshape(8,-1),axis=0)/2
  else:env=f2(x,self.q)
  return (1-self.f)*self.core.mass(z/self.sc)+self.f*env

def profiles():
 raw=np.genfromtxt(HERE/'field_profiles.csv',delimiter=',',names=True);d=json.loads((HERE/'profile_diagnostics.json').read_text());fine=d['fine']
 core=Density(raw['x'],raw['core_radial_density'],np.sqrt(2*abs(fine['mu'])))
 resp=Density(raw['x'],raw['response_radial_density'],d['response_tail_decay_rate'])
 pars=json.loads((HERE/'inputs/compression.json').read_text())['parameters']['envelope_mode_0']
 return {'scalar':core,'new_analytic':Mix(core,pars),'old_exact':Mix(core,pars,True),'mean_response':resp}

def fitamp(B,H,y,err,mask,upper):
 lower=max(0.,float(np.max(-B/H)));lower=lower*(1+1e-12)+(1e-10 if lower else 0.)
 if lower>=upper:return float('inf'),lower
 b=B[mask];hh=H[mask];yy=y[mask];ww=1/err[mask]**2
 def der(a):return np.sum(ww*hh*(1-yy/np.sqrt(np.maximum(b+a*hh,1e-30))))
 if der(lower)>=0:a=lower
 elif der(upper)<=0:a=upper
 else:a=brentq(der,lower,upper,xtol=1e-9,rtol=2e-13)
 cost=float(np.sum(((np.sqrt(np.maximum(b+a*hh,0))-yy)/err[mask])**2))
 return cost,a

def fit(d,p,mask,ngrid=65,starfactor=1.):
 r=d['r'];B=d['gas']+starfactor*d['star'];lo=np.log(r.min()/100);hi=np.log(r.max()*100);upper=.001*C2/p.depth
 def objective(s,ret=False):
  L=np.exp(s);x=r/L;H=p.mass(x)/x;cost,a=fitamp(B,H,d['y'],d['err'],mask,upper)
  return (cost,a,L,H) if ret else cost
 grid=np.linspace(lo,hi,ngrid);vals=np.array([objective(z) for z in grid]);ids=[i for i in range(1,ngrid-1) if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]]
 candidates=[(vals[0],lo),(vals[-1],hi)]
 for i in ids:
  z=minimize_scalar(objective,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':2e-9});candidates.append((z.fun,z.x))
 cost,s=min(candidates);cost,a,L,H=objective(s,True);raw=B+a*H
 assert np.all(raw>=-1e-7)
 pred=np.sqrt(np.maximum(raw,0));M=a*L/G
 rec={'L_kpc':float(L),'M_total_Msun':float(M),'amplitude_GM_over_L':float(a),'fit_chi2':cost,
  'length_bound':bool(abs(s-lo)<1e-6 or abs(s-hi)<1e-6),'weak_field_bound':bool(a>=upper*(1-1e-8)),
  'mass_outside_last_observation_fraction':float(1-p.mass(np.array([r.max()/L]))[0]),
  'central_potential_over_c2':float(a*p.depth/C2),'starfactor':starfactor}
 return rec,pred

def stellar_sensitivity(d,p,mask,start):
 # Same single relative stellar M/L factor and 0.1-dex prior in both models.
 sigma=.1*np.log(10);bounds=[np.log(.1),np.log(4.)]
 if p is None:
  def fun(q):return np.sum(((np.sqrt(np.maximum(d['gas']+np.exp(q)*d['star'],0))-d['y'])[mask]/d['err'][mask])**2)+(q/sigma)**2
  z=minimize_scalar(fun,bounds=bounds,method='bounded',options={'xatol':1e-9});sf=np.exp(z.x)
  return {'starfactor':float(sf),'penalized_objective':float(z.fun)},np.sqrt(np.maximum(d['gas']+sf*d['star'],0))
 r=d['r'];lo=np.log(r.min()/100);hi=np.log(r.max()*100);upper=.001*C2/p.depth
 def fun(v,ret=False):
  L=np.exp(v[0]);sf=np.exp(v[1]);x=r/L;hh=p.mass(x)/x;B=d['gas']+sf*d['star'];cost,a=fitamp(B,hh,d['y'],d['err'],mask,upper);cost+=(v[1]/sigma)**2
  return (cost,a,L,sf,hh,B) if ret else cost
 fits=[minimize(fun,[np.log(start['L_kpc'])+dl,0.],method='Powell',bounds=[(lo,hi),tuple(bounds)],options={'xtol':2e-6,'ftol':2e-8,'maxiter':120}) for dl in [-.8,0,.8]]
 z=min(fits,key=lambda z:z.fun);cost,a,L,sf,hh,B=fun(z.x,True)
 return {'starfactor':float(sf),'penalized_objective':float(cost),'L_kpc':float(L),'M_total_Msun':float(a*L/G)},np.sqrt(np.maximum(B+a*hh,0))

def main():
 data=load();ps=profiles();states=[];predictions=[];scores=[];cache={};refine=[]
 for ix,(name,d) in enumerate(data.items()):
  for phase in ['inner','all']:
   mask=~d['outer'] if phase=='inner' else np.ones(len(d['r']),bool)
   bp=np.sqrt(np.maximum(d['B'],0));models={'baryons':({},bp)}
   for mod,p in ps.items():
    st,pred=fit(d,p,mask);models[mod]=(st,pred)
    states.append({'galaxy':name,'phase':phase,'model':mod,**st})
   if phase=='inner':
    for mod,p in [('profiled_baryons',None),('profiled_new',ps['new_analytic'])]:
     st,pred=stellar_sensitivity(d,p,mask,models['new_analytic'][0]);models[mod]=(st,pred)
    cache[name]=models
   for mod,(st,pred) in models.items():
    for j in range(len(pred)):predictions.append({'galaxy':name,'phase':phase,'model':mod,'index':j,'r_kpc':float(d['r'][j]),'observed':float(d['y'][j]),'sigma':float(d['err'][j]),'predicted':float(pred[j]),'outer':bool(d['outer'][j]),'quality_sample':bool(d['quality']),'negative_catalogue_baryon_force':bool(d['B'][j]<0)})
    scoremask=d['outer'] if phase=='inner' else np.ones(len(pred),bool)
    res=(pred-d['y'])[scoremask];sig=d['err'][scoremask]
    scores.append({'galaxy':name,'phase':phase,'model':mod,'quality_sample':bool(d['quality']),'group':d['group'],'rows':int(sum(scoremask)),'loss':float(np.mean((res/sig)**2)),'mae':float(np.mean(abs(res))),'bias':float(np.mean(res)),'chi2':float(np.sum((res/sig)**2)), 'starfactor':float(st.get('starfactor',1.))})
  if ix%20==0:print('fitted',ix+1,'/',len(data),name,flush=True)
 csvout('fit_states.csv',states);csvout('predictions.csv',predictions);csvout('galaxy_scores.csv',scores)
 # Independent denser scale scan on original six examples and every bound-active fit.
 targets={'DDO154','UGC11557','F568-3','NGC4559','UGC06614','NGC2955'}|{s['galaxy'] for s in states if s['length_bound'] or s['weak_field_bound']}
 for name in sorted(targets):
  d=data[name];st,pr=fit(d,ps['new_analytic'],~d['outer'],ngrid=129);old=cache[name]['new_analytic'];refine.append({'galaxy':name,'chi2_change':float(st['fit_chi2']-old[0]['fit_chi2']),'max_speed_change':float(np.max(abs(pr-old[1])))})
 csvout('optimization_refinement.csv',refine)
 z=np.geomspace(1e-7,1e6,3001);checks={'all_galaxies':175,'all_rows':3391,'quality_galaxies':131,'quality_rows':sum(len(d['r']) for d in data.values() if d['quality']),
  'negative_baryon_force_rows':sum(int(sum(d['B']<0)) for d in data.values()),'no_galaxy_dropped':True,
  'optimization_refinement_max_speed_change':max(r['max_speed_change'] for r in refine),
  'profiles':{n:{'minimum_cdf':float(min(p.mass(z))),'cdf_at_1e6':float(p.mass(np.array([1e6]))[0]),'monotone':bool(np.all(np.diff(p.mass(z))>=-1e-14)),'unit_potential_depth':float(p.depth)} for n,p in ps.items()},
  'new_vs_exact_max_relative_mass_error':float(np.max(abs(ps['new_analytic'].mass(z)/ps['old_exact'].mass(z)-1))),
  'core_explicit_tail_mass_fraction':float(ps['scalar'].tailmass/ps['scalar'].total),'response_explicit_tail_mass_fraction':float(ps['mean_response'].tailmass/ps['mean_response'].total)}
 assert all(x['monotone'] for x in checks['profiles'].values())
 savejson('verification.json',checks)
 summarize(data,scores,states)
 print(json.dumps(checks,indent=2),flush=True)

def holm(vals):
 v=np.array(vals);order=np.argsort(v);out=np.zeros(len(v));cur=0.
 for j,i in enumerate(order):cur=max(cur,(len(v)-j)*v[i]);out[i]=min(1.,cur)
 return out

def inference(diff,groups,draws=99999,seed=2026092951):
 diff=np.asarray(diff);unique=sorted(set(groups));mem=np.array([[int(x==g) for x in groups] for g in unique]);sums=mem@diff;counts=mem.sum(1);obs=float(diff.mean());rng=np.random.default_rng(seed);samples=[]
 for j in range(0,draws,2048):
  m=rng.multinomial(len(unique),np.full(len(unique),1/len(unique)),size=min(2048,draws-j));samples.append((m@sums)/(m@counts))
 samp=np.concatenate(samples);k=int(np.sum(abs(samp-obs)>=abs(obs)));p=(k+1)/(draws+1)
 gd=sums/counts;tol=1e-7;wins=int(sum(gd < -tol));losses=int(sum(gd>tol));signp=binomtest(wins,wins+losses,.5).pvalue
 return {'mean_difference':obs,'bootstrap95':np.quantile(samp,[.025,.975]).tolist(),'two_sided_p':float(p),'equivalent_Z':float(norm.isf(p/2)),
  'bootstrap_exceedances':k,'replications':draws,'minimum_reportable_p':1/(draws+1),'groups':len(unique),'group_wins':wins,'group_losses':losses,'group_ties':len(unique)-wins-losses,
  'group_sign_p':float(signp),'group_sign_Z':float(norm.isf(signp/2)),'galaxy_wins':int(sum(diff < -tol)),'galaxy_losses':int(sum(diff>tol))}

def summarize(data,scores,states):
 record={'scope':'Conditional retrospective halo-template comparison. The added fields still gravitate through GR. This is not a test against all GR models.', 'summaries':{},'inference':{},'sensitivity':{}}
 lookup={(s['galaxy'],s['phase'],s['model']):s for s in scores}
 for sample,names in [('all175',list(data)),('quality131',[n for n,d in data.items() if d['quality']]),('admissible174',[n for n,d in data.items() if np.all(d['B']>=0)]),('admissible_quality130',[n for n,d in data.items() if d['quality'] and np.all(d['B']>=0)])]:
  sums=[]
  for phase in ['inner','all']:
   mods=['baryons','scalar','new_analytic','old_exact','mean_response']+(['profiled_baryons','profiled_new'] if phase=='inner' else [])
   for model in mods:
    ss=[lookup[(n,phase,model)] for n in names];sums.append({'phase':phase,'model':model,'galaxies':len(names),'rows':sum(s['rows'] for s in ss),'mean_loss':float(np.mean([s['loss'] for s in ss])),'median_loss':float(np.median([s['loss'] for s in ss])),'mae':float(np.mean([s['mae'] for s in ss])),'total_chi2':float(sum(s['chi2'] for s in ss))})
  record['summaries'][sample]=sums
  pairs=[('new_analytic','baryons'),('mean_response','baryons'),('new_analytic','scalar'),('mean_response','scalar')]
  tests=[]
  for a,b in pairs:
   diff=[lookup[(n,'inner',a)]['loss']-lookup[(n,'inner',b)]['loss'] for n in names]
   tests.append({'model':a,'baseline':b,**inference(diff,[data[n]['group'] for n in names])})
  for test,p,sp in zip(tests,holm([x['two_sided_p'] for x in tests]),holm([x['group_sign_p'] for x in tests])):
   test['Holm_p']=float(p);test['Holm_Z']=float(norm.isf(p/2));test['Holm_group_sign_p']=float(sp);test['Holm_group_sign_Z']=float(norm.isf(sp/2))
  record['inference'][sample]=tests
  diff=[lookup[(n,'inner','profiled_new')]['loss']-lookup[(n,'inner','profiled_baryons')]['loss'] for n in names]
  record['sensitivity'][sample]={'common_stellar_ML_profiling':inference(diff,[data[n]['group'] for n in names])}
 record['bounds']={mod:{'length_hits':sum(s['length_bound'] for s in states if s['phase']=='inner' and s['model']==mod),'weak_field_hits':sum(s['weak_field_bound'] for s in states if s['phase']=='inner' and s['model']==mod)} for mod in ['scalar','new_analytic','old_exact','mean_response']}
 savejson('RESULTS.json',record)
 for sample in ['all175','quality131']:
  print(sample,json.dumps(record['inference'][sample],indent=2),flush=True)

if __name__=='__main__':main()
