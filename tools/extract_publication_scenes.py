"""Extract immutable publication scenes and scientific coordinates from reference renderers.

This maintenance tool is not used by the .NET application. It copies the supplied
reference tree to an isolated temporary directory, runs only figure renderers,
captures their typed SVG vector primitives and all plotted artist coordinates,
and records source hashes. The original reference tree is never written.

Requires Python 3.12 with numpy, scipy, pandas, matplotlib and sympy.
Run: python extract_publication_scenes.py --reference <reference-root> --output <data/publication/figures>
Use --packages <directory> only if the packages are not on Python's normal path.
"""
from pathlib import Path
import argparse, gzip, hashlib, io, json, math, os, re, runpy, shutil, subprocess, sys, tempfile
import xml.etree.ElementTree as ET

COMMANDS = [
    ['figure_sources/plot_review_diagnostics.py'],
    ['figure_sources/physical_mechanism_tests.py'],
    ['figure_sources/plot_pointwise.py'],
    ['figure_sources/plot_extension.py'],
    ['numerics/bridge_halo_2026_09_26/plot_bridge_results.py','--gallery','--statistics','--diagnostics'],
    ['numerics/persistent_halo_population_2026_09_26/plot_population.py'],
    ['consolidation_2026_09_26/plot_consolidation.py'],
    ['numerics/inverse_field_configuration_2026_09_26/plot_configuration.py'],
    ['numerics/inverse_closure_2026_09_26/plot_closure.py'],
    ['completeness_review_2026_09_26/reproduction/plot_compact_physical.py'],
    ['completeness_review_2026_09_26/plot_inverse_compatibility.py','--root','.','--output','figures'],
]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def scene_node(node):
    return dict(tag=node.tag, attributes=dict(node.attrib), text=node.text,
                tail=node.tail, children=[scene_node(child) for child in node])

def worker(args):
    if args.packages: sys.path.insert(0, str(Path(args.packages).resolve()))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.collections import PathCollection, LineCollection, PolyCollection, QuadMesh
    from matplotlib.container import ErrorbarContainer
    root=Path(args.stage).resolve(); output=Path(args.output).resolve()
    manifest=json.loads((root/'PUBLICATION_EXTRACTION_INPUTS.json').read_text())
    expected=set(manifest['figures']); opened=set(); written=set(); captured={}
    os.chdir(root)
    command=COMMANDS[args.worker]
    script=root/command[0]
    sys.path.insert(0,str(script.parent))
    def audit(event, params):
        def relative(value):
            if not isinstance(value,(str,bytes,os.PathLike)): return None
            try:
                p=Path(value).resolve()
                return p.relative_to(root).as_posix() if p.is_relative_to(root) else None
            except (OSError,ValueError): return None
        if event=='shutil.copyfile':
            source_key=relative(params[0]); destination_key=relative(params[1])
            if source_key: opened.add(source_key)
            if destination_key: written.add(destination_key)
        elif event=='open':
            key=relative(params[0]);mode=params[1];flags=params[2]
            if not key: return
            if (mode and any(flag in mode for flag in ['w','a','+'])) or flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND): written.add(key)
            else: opened.add(key)
    sys.addaudithook(audit)

    def clean_number(x):
        try:
            value=float(x)
            return value if math.isfinite(value) else None
        except (ValueError,TypeError): return None

    def xy_rows(values):
        z=np.asarray(values)
        return [[clean_number(v) for v in row[:2]] for row in z.reshape((-1,z.shape[-1]))] if z.size else []

    def artist_label(artist):
        value=artist.get_label()
        return '' if value.startswith('_') else value

    def axis_description(axis):
        scale=axis._scale
        return dict(type=axis.get_scale(),base=clean_number(getattr(scale,'base',None)),
                    linthresh=clean_number(getattr(scale,'linthresh',None)),
                    linscale=clean_number(getattr(scale,'linscale',None)))

    def coordinates(fig):
        axes=[];series=[]
        for i,ax in enumerate(fig.axes):
            desc=dict(index=i,title=ax.get_title() or ax.get_title(loc='left'),
                      xlabel=ax.get_xlabel(),ylabel=ax.get_ylabel(),
                      xlim=list(map(float,ax.get_xlim())),ylim=list(map(float,ax.get_ylim())),
                      xscale=axis_description(ax.xaxis),yscale=axis_description(ax.yaxis),
                      position=list(map(float,ax.get_position().bounds)),
                      xticks=[dict(value=clean_number(v),label=t.get_text()) for v,t in zip(ax.get_xticks(),ax.get_xticklabels())],
                      yticks=[dict(value=clean_number(v),label=t.get_text()) for v,t in zip(ax.get_yticks(),ax.get_yticklabels())],
                      annotations=[t.get_text() for t in ax.texts],axis_visible=ax.axison)
            axes.append(desc)
            def add(artist,kind,points,transform=None,**extra):
                # Keep raw source coordinates and declare their original transform;
                # mixed data/axes reference lines must not be mislabelled as data.
                coord='data' if transform is not None and transform==ax.transData else 'artist'
                record=dict(axis=i,kind=kind,label=artist_label(artist),coordinate_system=coord,
                            transform=str(transform) if coord=='artist' else 'transData',
                            points=xy_rows(points),**extra)
                series.append(record)
            for line in ax.lines:
                add(line,'line',line.get_xydata(),line.get_transform(),
                    linestyle=str(line.get_linestyle()),marker=str(line.get_marker()))
            for collection in ax.collections:
                if isinstance(collection,PathCollection):
                    extra={}
                    if collection.get_array() is not None:
                        extra['color_values']=[clean_number(x) for x in collection.get_array()]
                    add(collection,'scatter',collection.get_offsets(),collection.get_offset_transform(),**extra)
                elif isinstance(collection,LineCollection):
                    for j,segment in enumerate(collection.get_segments()):
                        add(collection,'segment',segment,collection.get_transform(),part=j)
                elif isinstance(collection,QuadMesh):
                    # Colorbar cells are included exactly in the vector scene.
                    add(collection,'mesh',collection.get_coordinates().reshape((-1,2)),collection.get_transform())
                else:
                    for j,path in enumerate(collection.get_paths()):
                        add(collection,'polygon',path.vertices,collection.get_transform(),part=j)
            for j,patch in enumerate(ax.patches):
                vertices=patch.get_transform().transform(patch.get_path().vertices)
                values=ax.transData.inverted().transform(vertices)
                add(patch,'patch',values,ax.transData,part=j)
            for j,container in enumerate(ax.containers):
                if not isinstance(container,ErrorbarContainer): continue
                data_line,cap_lines,bar_lines=container.lines
                # Error bars also appear as line/segment artists. This explicit
                # record makes their relationship inspectable without SVG parsing.
                if data_line is not None:
                    add(data_line,'errorbar_centers',data_line.get_xydata(),data_line.get_transform(),
                        has_x_errors=container.has_xerr,has_y_errors=container.has_yerr)
        return axes,series

    original_save=Figure.savefig
    def save(self, fname, *positional, **kwargs):
        path=Path(fname) if isinstance(fname,(str,os.PathLike)) else None
        name=path.stem if path else ''
        if name=='shared_calibration': name='physical_population'
        if name in expected and name not in captured:
            # Use the same save bounding box and draw/layout pass as the source.
            stream=io.BytesIO()
            opts=dict(kwargs)
            opts.pop('format',None);opts['metadata']={'Date':None}
            original_save(self,stream,format='svg',**opts)
            node=ET.fromstring(stream.getvalue())
            axes,series=coordinates(self)
            artifact=dict(schema_version=1,name=name,
                scope='Frozen figure replay: retained vector scene and plotted artist coordinates; no refitting, inference, or evolution.',
                svg=scene_node(node),axes=axes,series=series,
                source_script=command[0],source_command=command,
                extraction_versions=dict(python=sys.version.split()[0],matplotlib=matplotlib.__version__,numpy=np.__version__))
            captured[name]=artifact
        # Some archived scripts hash/copy their output. Preserve that behavior in
        # the isolated stage only. Never direct their outputs at the reference.
        return original_save(self,fname,*positional,**kwargs)
    Figure.savefig=save
    sys.argv=[str(script),*command[1:]]
    runpy.run_path(str(script),run_name='__main__')
    # Renderers reopen their own SVG/PDF/PNG outputs to hash or copy them.
    # Those are rendered artifacts, never scientific inputs to these commands.
    provenance={key:manifest['files'][key] for key in sorted(opened-written)
                if key in manifest['files'] and Path(key).suffix.lower() not in {'.svg','.pdf','.png'}}
    for key,expected_hash in provenance.items():
        if sha(root/key)!=expected_hash:
            raise RuntimeError('Renderer read a staged input changed from the reference: '+key)
    provenance[command[0]]=manifest['files'][command[0]]
    results=[]
    for name,artifact in captured.items():
        artifact['inputs_sha256']=provenance
        artifact['archived_pdf_sha256']=manifest['figures'][name]
        dest=output/(name+'.scene.json.gz')
        payload=json.dumps(artifact,separators=(',',':'),allow_nan=False).encode('utf-8')
        with dest.open('wb') as target:
            with gzip.GzipFile(filename='',fileobj=target,mode='wb',mtime=0) as gz: gz.write(payload)
        results.append(dict(name=name,file=dest.name,sha256=sha(dest),axes=len(artifact['axes']),
             series=len(artifact['series']),coordinate_rows=sum(len(s['points']) for s in artifact['series']),
             source_script=command[0],archived_pdf_sha256=manifest['figures'][name]))
    write_json(output/('worker_'+str(args.worker)+'.json'),results)
    print('CAPTURED: '+', '.join(captured))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--packages',type=Path)
    ap.add_argument('--worker',type=int)
    ap.add_argument('--stage',type=Path)
    args=ap.parse_args()
    if args.worker is not None: worker(args);return
    if args.reference is None: ap.error('--reference is required')
    source=args.reference.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    if output==source or output.is_relative_to(source): ap.error('Output must be outside the reference tree')
    texts='\n'.join((source/n).read_text(encoding='utf-8') for n in ['Cubic_Scalar_GR_Intersection.tex','Cubic_Scalar_GR_Supplement.tex'])
    names=sorted(set(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{figures/([^}]+)\.pdf\}',texts)))
    if len(names)!=25: raise ValueError(f'Expected 25 current manuscript figures; found {len(names)}')
    script=Path(__file__).resolve()
    with tempfile.TemporaryDirectory(prefix='publication-extraction-') as temp:
        stage=Path(temp)/'reference'
        shutil.copytree(source,stage,ignore=shutil.ignore_patterns('__pycache__'))
        files={p.relative_to(source).as_posix():sha(p) for p in source.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        metadata=dict(figures={n:sha(source/'figures'/(n+'.pdf')) for n in names},files=files)
        write_json(stage/'PUBLICATION_EXTRACTION_INPUTS.json',metadata)
        env=os.environ.copy();env.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')
        for i,command in enumerate(COMMANDS):
            call=[sys.executable,str(script),'--worker',str(i),'--stage',str(stage),'--output',str(output)]
            if args.packages: call+=['--packages',str(args.packages.resolve())]
            log=output/('extraction_'+str(i)+'.log')
            print('Extracting '+command[0],flush=True)
            with log.open('w',encoding='utf-8') as stream:
                result=subprocess.run(call,cwd=stage,env=env,stdout=stream,stderr=subprocess.STDOUT)
            if result.returncode: raise RuntimeError(f'Renderer failed; inspect {log}')
        figures=[]
        for i in range(len(COMMANDS)):
            p=output/('worker_'+str(i)+'.json');figures+=json.loads(p.read_text());p.unlink()
        if sorted(row['name'] for row in figures)!=names: raise ValueError('Coverage mismatch')
        write_json(output/'manifest.json',dict(schema_version=1,scope='Frozen publication vector scenes with plotted coordinates. Source scripts regenerate scenes after changed science inputs. The .NET replay does not refit or rerun inference.',
          figure_count=len(figures),figures=sorted(figures,key=lambda r:r['name']),extractor_sha256=sha(script),
          manuscript_sha256={n:sha(source/n) for n in ['Cubic_Scalar_GR_Intersection.tex','Cubic_Scalar_GR_Supplement.tex']}))
        # Verify extraction never changed any source file.
        for key,expected in files.items():
            if sha(source/key)!=expected: raise RuntimeError('Reference changed during extraction: '+key)
        print(f'Extracted all {len(figures)} figures; reference hashes unchanged.',flush=True)

if __name__=='__main__': main()
