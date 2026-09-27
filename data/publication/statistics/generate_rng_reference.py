"""Regenerate NumPy PCG64 golden streams used by the native publication statistics.
Requires NumPy; not needed to run the .NET application.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

root = Path(__file__).resolve().parent
streams = []
for seed, galaxies, dtype, batch in [
    (2026092601,131,"int64",2000), (2026092602,131,"int8",5000),
    (2026092604,130,"int64",2000), (2026092611,131,"int64",2000),
    (2026092612,131,"int8",5000), (2026092614,130,"int64",2000),
]:
    rng = np.random.default_rng(seed)
    state = rng.bit_generator.state["state"]
    raw = [str(int(v)) for v in rng.bit_generator.random_raw(12)]
    rng = np.random.default_rng(seed)
    digest = hashlib.sha256()
    upper = 2 if dtype == "int8" else galaxies
    for start in range(0,99999,batch):
        values = rng.integers(0,upper,size=(min(batch,99999-start),galaxies),dtype=dtype)
        digest.update(values.astype(np.uint8).tobytes())
    streams.append(dict(seed=seed,state=str(state["state"]),increment=str(state["inc"]),
                        galaxies=galaxies,dtype=dtype,batch_draws=batch,draws=99999,
                        upper=upper,raw_uint64=raw,full_stream_sha256_uint8=digest.hexdigest()))
(root/"numpy_rng_reference.json").write_text(json.dumps(dict(numpy_version=np.__version__,streams=streams),indent=2)+"\n")
