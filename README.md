# Dark Universe — C# console reconstruction

One .NET 10 console application using Math.NET Numerics 5.0.0. No Python runtime requirement.

It imports SPARC, recomputes statistics from frozen predictions, and writes SVG plots with CSV sidecars. The current analysis is the **25 September 2026 pointwise comparison**; earlier manuscript analyses remain available.

## Run

From this folder, reproduce the updated analysis and its three figures:

```powershell
dotnet run -c Release -- pointwise
```

Run everything, including all checks and 47 figures:

```powershell
dotnet run -c Release -- verify
```



## Results

Open these generated files after a run:

Output\Contents


Additional tables and bootstrap arrays are in `output/pointwise_comparison/results/`. Earlier sign tests are in `output/significance/`; they concern win frequency rather than mean error magnitude.

SPARC is Lelli, McGaugh & Schombert (2016). Gas already includes helium, negative components represent signed force, and undefined baryonic circular speeds remain undefined. Observations and quoted errors are never rescaled a second time.

The application reconstructs saved-fit results; it performs no new population fit, nuisance refit, model selection, or field evolution. SVG typography differs from Matplotlib; numerical agreement is the reproduction target. Original manuscripts and source analyses remain unchanged.
