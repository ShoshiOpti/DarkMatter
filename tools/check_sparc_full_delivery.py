"""Developer delivery audit of native predictions, vector plots and portable output.

Standard library only. Run from any directory after source and published runs.
"""
from pathlib import Path
import csv,json,hashlib,xml.etree.ElementTree as ET
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/sparc_dotnet'
def latest(folder):
    ptr=json.loads((folder/'latest.json').read_text())
    assert ptr['passed'] and ptr['command']=='sparc-full'
    return folder/Path(ptr['gallery']).parent/'sparc_full'
def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=latest(OUT);portable=latest(OUT/'portable_run')
destination=ROOT/'DarkMatter/verification/sparc_full_upgrade_validation.json'
previous=json.loads(destination.read_text()) if destination.exists() else None
# Scientific staging must match the latest research, not merely its own manifest.
staged=ROOT/'DarkMatter/data/sparc_full';research=ROOT/'research/sparc_full_profile_2026_09_29'
fresh=[]
for p in sorted(staged.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(staged)
    original=research/(rel.name if rel.parts[0]=='reference' else rel)
    if original.is_file():
        assert sha(p)==sha(original),f'Stale research input: {rel}'
        fresh.append(rel.as_posix())
    else:
        assert rel.as_posix() in ['random_contract.json','profile_contract.json'],f'Untracked staged input: {rel}'
# A successful old executable is not evidence that the current publication is fresh.
published=OUT/'publish';project=ROOT/'DarkMatter';published_checked=[]
for p in [project/'README.md',*(project/'tools').rglob('*'),*(project/'data').rglob('*'),*(project/'research').rglob('*')]:
    if not p.is_file():continue
    relative=p.relative_to(project)
    assert (published/relative).is_file(),f'Missing published file: {relative}'
    assert sha(p)==sha(published/relative),f'Stale published file: {relative}'
    published_checked.append(relative.as_posix())
assembly=project/'bin/Release/net10.0/DarkUniverse.Console.dll'
assert sha(assembly)==sha(published/'DarkUniverse.Console.dll')
for run in [source,portable]:
    run_report=json.loads((run.parent/'run_report.json').read_text())
    assert run_report['assembly_sha256']==sha(assembly),'Run used a stale assembly'
complete_pointer=json.loads((OUT/'complete_verification/latest.json').read_text())
assert complete_pointer['passed'] and complete_pointer['command']=='verify'
complete_run=OUT/'complete_verification'/Path(complete_pointer['gallery']).parent
complete_report=json.loads((complete_run/'run_report.json').read_text())
assert complete_report['assembly_sha256']==sha(assembly)
assert complete_report['full_sparc']['status']=='pass'
gallery=(complete_run/'index.html').read_text()
assert 'sparc_full/index.html' in gallery and 'publication/figures/' in gallery
compared=[]
for p in source.rglob('*'):
    if p.is_file():
        name=p.relative_to(source);assert sha(p)==sha(portable/name),name
        compared.append(name.as_posix())
pred=rows(source/'predictions.csv')
ix={(r['galaxy'],r['phase'],r['model'],int(r['index'])):r for r in pred}
mapping={'Baryons only':('inner','baryons'),'Scalar (inner fit)':('inner','scalar'),
         'Analytic envelope (inner fit)':('inner','new_analytic'),
         'Mean response (inner fit)':('inner','mean_response'),
         'Analytic envelope (all fitted)':('all','new_analytic')}
panels={};coordinate_count=0
for file in sorted((source/'figures').glob('atlas_*.series.csv')):
    for r in rows(file):
        name=r['title'].split(' (outside quality subset)')[0]
        panels.setdefault(name,set()).add(file.name)
        if r['label'] not in mapping:continue
        phase,model=mapping[r['label']]
        expected=ix[name,phase,model,int(r['index'])]
        assert float(r['x'])==float(expected['r_kpc'])
        assert float(r['y'])==float(expected['predicted'])
        coordinate_count+=1
assert len(panels)==175 and all(len(v)==1 for v in panels.values())
assert coordinate_count==5*3391
svgs=list((source/'figures').glob('*.svg'));assert len(svgs)==32
for p in svgs:
    root=ET.parse(p).getroot();assert root.tag=='{http://www.w3.org/2000/svg}svg'
    assert float(root.attrib['width'])>0 and float(root.attrib['height'])>0
visual='Not inspected by this script; numerical/coordinate validation only.'
if previous:
    old=Path(previous['native_results'])/'figures'
    if all((old/p.name).is_file() and sha(p)==sha(old/p.name) for p in svgs):
        visual='Current SVGs match the previous visually checked figure set byte-for-byte.'
tests=json.loads((OUT/'regression_tests.json').read_text(encoding='utf-8-sig'))
assert tests['passed'] and tests['tests']==43
record={'status':'pass','verified_at_utc':datetime.now(timezone.utc).isoformat(),
        'regression_tests':tests['tests'],'portable_output_files_identical':len(compared),
        'latest_research_files_verified':len(fresh),'published_support_files_verified':len(published_checked),
        'current_assembly_sha256':sha(assembly),'complete_verify_run':str(complete_run),
        'atlas_unique_galaxies':175,'atlas_model_coordinates_verified':coordinate_count,'valid_SVG_figures':len(svgs),
        'native_results':str(source),'portable_results':str(portable),
        'visual_review':visual,
        'numerical_validation':json.loads((source/'verification.json').read_text()),
        'source_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in [
            ROOT/'DarkMatter/FullSparcProfiles.cs',ROOT/'DarkMatter/FullSparcStatistics.cs',ROOT/'DarkMatter/FullSparcStudy.cs',
            ROOT/'DarkMatter/Cli.cs',ROOT/'DarkMatter/Program.cs',ROOT/'DarkMatter.Tests/FullSparcTests.cs',
            ROOT/'DarkMatter/data/input_manifest.json',ROOT/'DarkMatter/README.md',
            ROOT/'DarkMatter/DarkUniverse.Console.csproj',ROOT/'DarkUniverse.slnx',ROOT/'global.json',
            ROOT/'DarkMatter/tools/check_sparc_full_delivery.py']}}
destination.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k not in ['numerical_validation','source_sha256']},indent=2))
