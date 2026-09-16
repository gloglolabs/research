# L9C and T17C in the CHGA leader of engineered insulin

Research note, 15 September 2026. GloGlo Labs.

[Read the note](research-note.md) or [download the PDF](chga-ins-research-note.pdf).

## Full paper

[Read the manuscript](article.md) or [download the full paper](chga-ins-leader-variants.pdf):
**Computational prioritization of chromogranin A signal peptide variants for engineered insulin**,
GloGlo Labs, working paper version 1.0, 15 September 2026.

The paper provides the experimental context from Kobaisi et al., full computational
methods, results, limitations, four figures and eleven references. It reports the
L9C single-substitution control, full-precursor counts, residual affinity tradeoff
and the 3% threshold reversal alongside the selected double-variant result.

`manuscript.md` is the authoring source. `build_manuscript.py` combines it with
`references.json` and `figures/captions.md` to generate `article.md`,
`references.bib` and an editable DOCX. With the dependencies in
`requirements-manuscript.txt`, run from this study directory:

```sh
python3 figures/build_source_data.py --bundle-only
python3 figures/plot_figures.py
python3 build_manuscript.py
```

PDF conversion requires a separate document renderer. The website serves a copy
of the final PDF at `/research/papers/chga-ins-leader-variants.pdf`; its build does
not depend on this submodule. See [figure documentation](figures/README.md) for
source-data provenance and rendering details.

This computational analysis screens substitutions in the reconstructed
CHGA1-18 + INS25-110 precursor. L9C + T17C reduces predicted leader/junction
presentation counts across the selected comparisons while retaining the predicted
18/19 cleavage boundary. The exact synthesized construct used by the source
study was unavailable; all sequence calculations refer to the stated reconstruction.

## Run

From this study directory, using Python 3.10 or later:

```sh
python3 reproducibility/reproduce.py --output recalculated
```

Expected result: `PASS`, 13 constructs and 34,762 final prediction rows. At the
strict 2% cutoff, leader/junction counts change 17 to 8 for NetMHCpan EL and
27 to 14 for MHCflurry with flanks. The runner also verifies both flank settings,
the 42-comparison selection, sequence identity and cleavage outputs.

This command reanalyses cached model predictions. It does not rerun model inference.
The bundle also retains related donor-ligand and inflammation analyses; the note's
central variant-selection result is computational, not a donor intervention result.

## Materials

| Path                                | Contents                                                             |
| ----------------------------------- | -------------------------------------------------------------------- |
| `reproducibility/reproduce.py`      | Supported offline analysis entry point                               |
| `reproducibility/constructs.fasta`  | Final precursor sequences                                            |
| `reproducibility/cached/`           | Saved predictor outputs                                              |
| `reproducibility/metadata/`         | Scan definitions and canonical sequence references                   |
| `reproducibility/tables/`           | Input tables and attributed experimental excerpts                    |
| `reproducibility/results/`          | Reference results for comparison                                     |
| `reproducibility/provenance/`       | Source URLs, versions, requests and hashes                           |
| `reproducibility/archived_scripts/` | Original producer and analysis scripts, preserved for methods review |

See the [bundle README](reproducibility/README.md),
[data dictionary](reproducibility/DATA_DICTIONARY.md) and
[provenance](reproducibility/PROVENANCE.md) for details. Historical scripts depend
on the original working layout and, in some cases, external predictors, weights
or source workbooks. They are not the portable execution interface.

## Publication history

Initial public package: original analysis scripts, cached inputs and reference
CSV results preserved byte-for-byte from the reviewed bundle. Packaging README
and provenance licensing text were updated for public distribution; input hashes
were regenerated. The research-note PDF presents the analysis with a formal source-paper citation. An old private
packaging receipt was omitted; current verification is executable from the root.

The source study is [Kobaisi et al., Cell Reports 45, 117941 (2026)](https://doi.org/10.1016/j.celrep.2026.117941).
Computational analysis and drafting used Codex agents. Cite this study's exact
commit or release together with the relevant predictor and source-study references.
