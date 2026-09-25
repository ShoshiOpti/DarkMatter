from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/numerics/pointwise_comparison_2026_09_25/random_streams.json'
OUT.parent.mkdir(parents=True,exist_ok=True)
configs=[(2026092501,131,2000,'int64'),(2026092502,2,5000,'int8'),(2026092503,130,5000,'int64')]
streams=[]
for seed,upper,batch,dtype in configs:
    rng=np.random.default_rng(seed)
    state=rng.bit_generator.state['state']
    raw=rng.bit_generator.random_raw(16).tolist()
    rng=np.random.default_rng(seed)
    h=hashlib.sha256();first=[];last=[];boundaries=[]
    for start in range(0,99999,batch):
        count=min(batch,99999-start)* (131 if seed!=2026092503 else 130)
        values=rng.integers(0,upper,size=count,dtype=np.dtype(dtype))
        h.update(values.astype('u1').tobytes())
        if start==0:first=values[:64].tolist()
        if start in (0,batch,99999//batch*batch):
            boundaries.append(dict(start_draw=start,first=values[:16].tolist(),last=values[-16:].tolist()))
        last=values[-64:].tolist()
    streams.append(dict(seed=seed,state=str(state['state']),increment=str(state['inc']),raw_uint64=[str(x) for x in raw],upper=upper,dtype=dtype,batch_draws=batch,galaxies=131 if seed!=2026092503 else 130,draws=99999,first=first,last=last,boundaries=boundaries,full_stream_sha256_uint8=h.hexdigest()))
checks=[]
for dtype,upper,lengths in [('int8',127,[3,5,17,64,1,7]),('int64',1073741825,[1,3,2,64])]:
    rng=np.random.default_rng(2026092501)
    batches=[rng.integers(0,upper,size=n,dtype=np.dtype(dtype)).tolist() for n in lengths]
    checks.append(dict(seed=2026092501,dtype=dtype,upper=upper,batches=batches))
result=dict(numpy_version=np.__version__,generated_by='verification/pointwise_rng_reference.py',provenance='Initial states from numpy.random.default_rng(seed).bit_generator.state; draws from Generator.integers with the exact analysis batch sizes.',implementation_reference='https://github.com/numpy/numpy/blob/v'+np.__version__+'/numpy/random/src/distributions/distributions.c',streams=streams,rejection_and_buffer_checks=checks)
OUT.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
