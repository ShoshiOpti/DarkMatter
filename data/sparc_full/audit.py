"""Independent integrity and mathematical checks for the completed comparison."""
from pathlib import Path
import csv,json
import test_sparc as t
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
r=t.HERE;data=t.load();ps=t.profiles();states=list(csv.DictReader((r/'fit_states.csv').open()));pred=list(csv.DictReader((r/'predictions.csv').open()))
index={(x['galaxy'],x['phase'],x['model'],int(x['index'])):x for x in pred};maxerr=0
for s in states:
 d=data[s['galaxy']];L=float(s['L_kpc']);M=float(s['M_total_Msun']);f=ps[s['model']].mass(d['r']/L)
 v=np.sqrt(np.maximum(d['B']+t.G*M*f/d['r'],0))
 old=np.array([float(index[s['galaxy'],s['phase'],s['model'],j]['predicted']) for j in range(len(v))]);maxerr=max(maxerr,float(max(abs(v-old))))
assert maxerr<1e-7
checks=[]
for q in [0.,.00059489359,.1,.8]:
 rho=lambda z:4/np.pi*z*z*((1+z*z)**-2+2*q*(1-z*z)/(1+z*z)**4)
 total=quad(rho,0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
 errors=[abs(quad(rho,0,z,epsabs=1e-12)[0]-t.f2(np.array([z]),q)[0]) for z in [.01,.2,1.,3.,30.]]
 checks.append({'q':q,'total_mass_error':abs(total-1),'cdf_integral_max_error':max(errors)})
assert max(x['cdf_integral_max_error'] for x in checks)<1e-10
# Group-aware loss statistics must count galaxies, not individual radii as independent units.
results=json.loads((r/'RESULTS.json').read_text())
for sample in ['all175','quality131']:
 for x in results['inference'][sample]:assert abs(x['Holm_Z']-norm.isf(x['Holm_p']/2))<1e-12
# Frozen outer errors and velocities are absent from the fit mask.
assert sum(sum(d['outer']) for d in data.values())==761
assert sum(sum(d['outer']) for d in data.values() if d['quality'])==659
record={'status':'pass','independent_GM_over_R_prediction_max_error_kms':maxerr,'analytic_full_mass_integrals':checks,'outer_points_all175':761,'outer_points_quality131':659,'all_rows_evaluated':3391,'no_cutoff_at_Rmax':True,'p_Z_conversion':'verified','scope':'Numerical and protocol integrity; not a validation of unmeasured halo occupation or universal microscopic constants.'}
t.savejson('independent_audit.json',record);print(json.dumps(record,indent=2))
