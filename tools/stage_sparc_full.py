"""Developer-only staging of the 29 September study; not used by .NET at runtime.

Run with the study's NumPy/SciPy environment. Existing unrelated manifest entries
are preserved and verified, never silently refreshed.
"""
from pathlib import Path
import sys,json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'research/sparc_full_profile_2026_09_29'
sys.path.insert(0,str(SRC))
import test_sparc as t
import numpy as np
DST=ROOT/'DarkMatter/data/sparc_full'
DST.mkdir(parents=True,exist_ok=True)
for name in ['field_profiles.csv','profile_diagnostics.json','PROTOCOL.json','README.md',
             'FULL_FIT_SIGNIFICANCE.md','build_profiles.py','test_sparc.py','robustness.py',
             'check_mesh.py','audit.py','render_report.py','full_fit_significance.py',
             'mesh_sensitivity.json','guard_sensitivity.json','guard_sensitivity.csv','requirements.txt']:
    shutil.copy2(SRC/name,DST/name)
shutil.copytree(SRC/'inputs',DST/'inputs',dirs_exist_ok=True)
(DST/'reference').mkdir(exist_ok=True)
for name in ['predictions.csv','fit_states.csv','galaxy_scores.csv','RESULTS.json',
             'full_nuisance_states.csv','full_nuisance_predictions.csv','full_nuisance_scores.csv',
             'nuisance_sensitivity.json','full_fit_significance.json']:
    shutil.copy2(SRC/name,DST/'reference'/name)
seed=2026092951; gen=np.random.default_rng(seed)
initial=gen.bit_generator.state['state']
streams=[]
for groups in [148,111]:
    gen=np.random.default_rng(seed);digest=hashlib.sha256();first=None
    for start in range(0,99999,2048):
        counts=gen.multinomial(groups,np.full(groups,1/groups),size=min(2048,99999-start)).astype(np.uint8)
        if first is None:first=counts[0].tolist()
        digest.update(counts.tobytes())
    streams.append({'groups':groups,'draws':99999,'sha256_uint8_counts':digest.hexdigest(),'first_counts':first})
rng={'seed':seed,'state':str(initial['state']),'increment':str(initial['inc']),
     'numpy_version':np.__version__,'distribution':'sequential conditional binomial inversion; PCG64 double=(raw>>11)*2^-53',
     'streams':streams}
(DST/'random_contract.json').write_text(json.dumps(rng,indent=2)+'\n')
x=np.geomspace(1e-6,1e5,97)
gold={'radii':x.tolist(),'profiles':{name:{'depth':p.depth,'mass':p.mass(x).tolist()} for name,p in t.profiles().items()}}
(DST/'profile_contract.json').write_text(json.dumps(gold,indent=2)+'\n')
manifest_path=ROOT/'DarkMatter/data/input_manifest.json';manifest=json.loads(manifest_path.read_text())
for name,digest in manifest.items():
    if not name.startswith('sparc_full/'):
        assert hashlib.sha256((manifest_path.parent/name).read_bytes()).hexdigest()==digest,name
for p in sorted(DST.rglob('*')):
    if p.is_file():manifest[p.relative_to(manifest_path.parent).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
manifest_path.write_text(json.dumps(dict(sorted(manifest.items())),indent=2)+'\n')
print(json.dumps({'staged':str(DST),'files':len(list(DST.rglob('*'))),'random':rng},indent=2))
