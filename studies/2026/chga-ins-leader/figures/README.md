# Manuscript figures

Four original figures generated from cached prediction outputs and curated public donor data. Every figure is 7 inches wide. Vector PDF and SVG, and 300 dpi PNG versions are supplied; editable SVG text and embedded PDF fonts are retained.

## Reproduce the plots

Run `python plot_figures.py` in an environment with Matplotlib 3.11.2 and NumPy. The script consumes only the adjacent `figure*.csv` source-data files.

To independently rebuild and verify the plotted source data, run `python build_source_data.py --bundle-only` from the final package. This resolves cached inputs through `../reproducibility/provenance/bundle_origins.json` and the three small supplementary files in `source_inputs/`. It does not require the original research tree. Without `--bundle-only`, the script can also read the original investigation cache. This checks peptide coordinates, complete prediction grids, recomputed counts against saved summaries, and the exact donor detection mapping. It performs no network access or model inference.

`source_data_manifest.json` records input paths and checksums. `figure_manifest.json` records rendered output checksums, dimensions and software versions. Figure heights are 6.5, 7.2, 5.3 and 6.0 inches, respectively. `captions.md` uses image paths relative to the parent manuscript directory for insertion into the manuscript builder.

No synthetic observations, intensity-derived protection estimates, error bars or biological replicate claims are used.
