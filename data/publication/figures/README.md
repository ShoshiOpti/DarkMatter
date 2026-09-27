# Frozen manuscript figure scenes

This directory contains the 25 distinct figures referenced by the consolidated main paper and supplement. The compressed JSON files contain a typed SVG scene tree (paths, glyphs, axes, collections, annotations and legends), scientific artist coordinates, axis limits and transforms, source script identity, hashes of the actual source files read, and the hash of the corresponding archived publication PDF.

`PublicationPlots.Run` renders these scenes in .NET without Python, emits one `.series.csv` and `.axes.json` sidecar per figure, checks all input hashes and internal references, and writes `publication/figure_replay_report.json`. Output is isolated in `publication/figures` so frozen publication figures cannot overwrite figures drawn from newly computed pointwise statistics. There are 92 axes, 871 artist series and 71,902 coordinate rows. One original colourbar contains its original self-contained PNG; all other scene graphics are vector primitives.

This is **frozen figure replay**, not a new fit, bootstrap, equilibrium computation or field evolution. Changing predictions or numerical results does not automatically update a frozen scene. Fresh statistics and physical diagnostics are separate native application workflows.

To regenerate the scenes after changing scientific inputs, use the maintained Python reference renderers through the extraction tool:

```text
python tools/extract_publication_scenes.py --reference ../reference --output data/publication/figures
```

Run from `DarkMatter`. Python needs NumPy, SciPy, pandas, Matplotlib and SymPy; `--packages <directory>` is optional for an existing package directory. The extractor copies the reference tree to an isolated temporary folder, runs only the figure renderers, verifies all original source hashes remain unchanged, and deletes the temporary copy. It refuses output inside the reference directory. Its normal workflow never reads any former project directory.

`manifest.json` binds all 25 compressed scenes to their hashes, original PDFs, source scripts, counts, manuscript hashes and extraction recipe hash. The small `extraction_*.log` files preserve the figure renderers' own checks. Source hashes inside each scene refer to paths relative to the retained `reference` tree. For scientific review, inspect the uncompressed numeric sidecars emitted by the native application and the corresponding source CSV/NPZ files in `reference`.

Coordinates marked `data` are in scientific axis units. Coordinates marked `artist` retain their original Matplotlib transform, recorded in the series metadata; these include reference lines that mix data and axis coordinates. Error-bar segments and polygon vertices are explicit artists. A few paths contain missing values that deliberately break a curve; these are empty numeric CSV cells. The exact vector scene preserves all transformations, error bars, bands, colour values, labels and limits independently of the sidecar representation.
