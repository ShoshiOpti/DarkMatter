# Declared bootstrap seed

The source analysis `numerics/baryons_only/analyze_saved.py` declares `SEED=20260924` and `BOOT=20000`. Each paired sample resets `np.random.default_rng(SEED)` and draws indices with `rng.integers(0,n,size=(1000,n))` in twenty batches.

The adjacent JSON records the initial PCG64 state produced by NumPy SeedSequence for that fixed seed and independent raw and bounded test sequences. The state and test vectors were extracted during development with the bundled Python runtime and NumPy 2.3.5, using fresh generators for each expression:

```
np.random.default_rng(20260924).bit_generator.state
np.random.default_rng(20260924).bit_generator.random_raw(12)
np.random.default_rng(20260924).integers(0,153,12)
```

`BaryonBootstrap.cs` stores this fixed seed initialization and independently advances PCG64 in native C#. It uses all 20,000 resamples for each of the twenty pair/sample combinations; no saved random indices, bootstrap intervals, Python process, or Python runtime are used by the application. It checks both native sequences before computing intervals and compares every output field with the original table.

Algorithm references: [NumPy PCG64](https://github.com/numpy/numpy/blob/main/numpy/random/src/pcg64/pcg64.h), [NumPy bounded-integer generation](https://github.com/numpy/numpy/blob/main/numpy/random/src/distributions/distributions.c), and [Lemire, Fast Random Integer Generation in an Interval](https://arxiv.org/abs/1805.10941).

