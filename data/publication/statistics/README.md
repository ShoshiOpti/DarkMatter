# Frozen publication statistics inputs

These inputs make the native PublicationStatistics workflow independent of the retired source-and-software folder. Run "dotnet run --project DarkMatter -- publication" from the project root, or the documented statistics command for current and historical statistical families.

- pressure/predicted_rows.csv and charged/predicted_rows.csv retain all 3,034 observation keys for all seven required models. Each family includes its exact original protocol.
- envelope/saved_galaxy_scores.csv preserves the original score serialization solely for the historical exact-zero sign sensitivity. The .NET workflow first reconstructs and checks these scores against data/numerics/pointwise_comparison_2026_09_25/data/pointwise_all.csv. The primary 1e-4 and sensitivity 1e-8 signs use freshly reconstructed losses.
- The reference subdirectories contain independent original Python outputs and evaluator source. Output tables are read only after native calculation, to verify results; they do not supply calculated effects or bootstrap intervals.
- numpy_rng_reference.json records six NumPy PCG64 initial states and SHA-256 digests of every bootstrap draw. generate_rng_reference.py independently regenerates them with NumPy; Python is unnecessary for normal .NET runs.
- input_manifest.json hashes every retained input here. The original observational reference is also pinned by hash in the native implementation.

The pressure pairs/wild seeds are 2026092601 and 2026092602; its feasible-profiled-baryon subset uses 2026092604. Charged pairs/wild seeds are 2026092611 and 2026092612; its subset uses 2026092614. The latter seeds follow the protocol's base-plus-zero-based-comparator-index rule.

Both physical comparisons retain all six Holm slots. Explicit training failures, nonfinite outer predictions and negative outer speeds invalidate the galaxy/model block. Signs retain every galaxy; no continuous primary p is assigned to a failed full-sample contrast. Subset effects and intervals are labeled as conditional on the 130 jointly feasible galaxies.

All inference conditions on saved predictions, uses whole galaxies as clusters, and remains retrospective. No significance here establishes the cubic mechanism, halo formation, physical persistence, or a test of GR.
