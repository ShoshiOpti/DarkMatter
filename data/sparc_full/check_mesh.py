import json,csv,numpy as np
import test_sparc as t
from build_profiles import build
x,c,r,k,diag=build(599);core=t.Density(x,c,np.sqrt(2*abs(diag['mu'])));resp=t.Density(x,r,k)
data=t.load();saved=list(csv.DictReader((t.HERE/'galaxy_scores.csv').open()));lookup={(s['galaxy'],s['model']):float(s['loss']) for s in saved if s['phase']=='inner'}
rows=[]
for i,(n,d) in enumerate(data.items()):
 losses=[]
 for p in [core,resp]:
  st,pr=t.fit(d,p,~d['outer']);losses.append(float(np.mean(((pr-d['y'])[d['outer']]/d['err'][d['outer']])**2)))
 fine=lookup[n,'mean_response']-lookup[n,'scalar'];rows.append({'galaxy':n,'coarse_difference':losses[1]-losses[0],'fine_difference':fine,'difference_change':losses[1]-losses[0]-fine})
 if i%50==0:print('mesh audit',i+1,flush=True)
t.csvout('mesh_sensitivity.csv',rows)
record={}
for sample,names in [('all175',list(data)),('quality131',[n for n,d in data.items() if d['quality']])]:
 rr=[r for r in rows if r['galaxy'] in names];record[sample]={'coarse_mean_response_minus_scalar':float(np.mean([r['coarse_difference'] for r in rr])),'fine_mean_response_minus_scalar':float(np.mean([r['fine_difference'] for r in rr])),'change_in_mean_difference':float(np.mean([r['difference_change'] for r in rr])),'max_galaxy_difference_change':float(max(abs(r['difference_change']) for r in rr))}
t.savejson('mesh_sensitivity.json',record);print(json.dumps(record,indent=2))
