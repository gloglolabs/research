# Reproducibility bundle

**Canonical CHGA–INS leader and junction analysis.** This bundle accompanies the manuscript and contains generated prediction outputs, attributed source-table excerpts, construct sequences, and executable analysis code. This is the public companion package for the September 2026 research note.

## Run without installing packages

Requires Python 3.10 or later. From this directory:

```sh
python3 reproduce.py
```

The command uses only the Python standard library, makes no network requests, verifies the bundled inputs against `checksums.sha256`, and recreates the clean tables in `results/`. It finishes by printing `PASS` in a structured report and writes `results/validation.json`.

To choose a separate output directory:

```sh
python3 reproduce.py --output recalculated
```

The default computation was also tested after extracting the bundle into a separate directory, using Python isolated mode. Run the repository-level `python3 scripts/verify.py` to verify the published package against its reference tables.

## Main outputs

- `results/final_peptide_HLA_predictions.csv`: **34,762 rows**, 13 full precursor sequences × 382 windows × 7 HLA alleles. Includes NetMHCpan affinity/EL predictions and MHCflurry affinity/presentation outputs with and without full five-residue flanks. NetMHCpan predictions outside each 29-residue variant request are reused from the independently scanned full canonical precursor only after exact peptide identity is verified; `net_score_origin` records this distinction for every row.
- `results/reference_junction_comparison.csv`: 21 selected native/junction peptide–HLA comparisons, recalculated from both models’ raw reference outputs.
- `results/prediction_totals.csv`: all 13 construct totals by scope/model/context. The 504 signal/junction pairs per construct are 72 overlapping windows × 7 alleles.
- `results/objectives42.csv`, `results/criteria42_summary.csv`: seven alleles × two presentation models × three strict cutoffs, separately for both MHCflurry flank settings.
- `results/expanded_objectives42.csv`, `results/expanded_criteria42_summary.csv`: the 616-pair processing-flank sensitivity check.
- `results/leader_3_9_netmhcpan_summary.csv`, `results/position17_netmhcpan_summary.csv`, `results/single_mhcflurry_summary.csv`: both models’ staged single-residue screens, including baseline controls.
- `results/position17_objectives42.csv`, `results/position17_criteria42_summary.csv`: the position-17 comparison over 266 core-containing pairs and 546 core-or-flank-affected pairs.
- `results/donor_peptide_matrix.csv`: six published donor preparations × four curated INS peptides, with six positive donor–peptide observations.
- `results/donor_sensitivity.csv`, `results/donor_leave_one_out.csv`: blank/depth sensitivity and leave-one-donor-out calculations.
- `constructs.fasta`, `tables/constructs.csv`: the 13 final/control precursor sequences, substitutions, and SignalP cleavage outputs.

## What the runner validates

1. Complete peptide/HLA grids, unique keys, inclusive sequence coordinates, raw-model input/output alignment and exact processing flanks.
2. The 133 position 3–9 single substitutions, 19 position 17 alternatives and 8 selected doubles; L9C/G are the two early leader-only alternatives without NetMHCpan EL<0.5% hits.
3. L9C+T17C signal/junction counts: NetMHCpan 4/8/17→0/2/8 and MHCflurry full-flank 9/19/27→3/6/14 at 0.5%, 1%, 2%.
4. Its 42-criterion comparisons: 16 improved/26 unchanged/0 worse with flanks; 15/27/0 without. It is the only tested double passing both settings.
5. The 616-pair expanded analysis reproduces NetMHCpan counts 4/8/19 → 0/2/10 and MHCflurry counts 9/19/28 → 3/6/15, with the same 42-criterion comparison. Full-precursor and unchanged-INS checks include every aligned pair crossing the specified cutoffs.
6. Exact SignalP processed sequences equal INS25–110; the candidate differs from baseline only at CHGA positions 9 and 17.
7. Published donor-coordinate mapping, background/depth sensitivity, and the recorded B10–18 IFNα fold changes 11.25 and 1.44.

## Scope of reproducibility

The default runner **reanalyses saved predictor outputs**. It does not execute neural-network inference or claim to retrain the models. It uses small attributed excerpts of the donor and inflammation source tables; it does not independently reconstruct the authors' full mass-spectrometry search or background-filtering pipeline.

`cached/` contains the generated predictions used for calculation. `provenance/requests/` retains the successful NetMHCpan API requests for the staged scans. `provenance/mhcflurry_manifest.json` records code/runtime versions and per-model-file hashes, with the official weights URL. `cached/signalp/` retains the final 16-sequence batch, while `cached/signalp_position17/` retains the initial 21-sequence batch. Both include full submitted sequences and returned predictions. `archived_scripts/` preserves the original producer and analysis code for methodological inspection. Those historical scripts use the original task layout; **`reproduce.py` is the self-contained entry point**.

The 133 position-3–9 single alternatives were scored with a 33-residue input prefix sufficient to provide complete five-residue flanks for all retained signal/junction windows. The 20-choice position-17 scan and final 13-member panel were scored by MHCflurry using full 104-residue sequences.

The reference-proteome exact-match result is included with its query sequences, release, source checksums and matched controls. The default runner checks that recorded result; it does not rescan the omitted 82.7 million-residue reference corpus. `audit_reference_matches.py` can repeat the scan after those inputs are retrieved separately.

## Reading the result correctly

All designs start from the explicit reconstruction CHGA1–18 + INS25–110. The exact published synthesized construct was unavailable. Predictions are indexed to this reconstructed precursor, not to an experimentally verified sequence. Native endogenous CHGA is a separate protein and is not changed by these INS construct variants.

The 42 criteria concern **presentation-count thresholds**, not all binding affinities or all possible cutoffs. For example, the candidate creates a NetMHCpan A02 affinity prediction of 202.39 nM for VLACLLCAGQV, versus 654.57 nM for its aligned baseline; it remains outside the selected EL cutoffs. The candidate is not a global optimum, an experimental secretion result, or a demonstrated immune-escape outcome. Only 8 selected double mutants were tested.

Donor HLA labels from the source supplement are computational assignments. Positive detection is a curated ligand observation, not a direct killing assay. Recorded ion intensities are kept for provenance and are not pooled across donors or treated as molecule counts.

## Source material and licensing

Generated predictions, analysis code, derived tables and compact provenance are included. Original papers/PDFs, complete third-party experimental workbooks, downloaded environments, model weights and full reference proteomes are excluded. `provenance/retrieval_manifest.json` provides URLs, exact source hashes and purposes; `PROVENANCE.md` gives attribution and version details.

See the repository LICENSE.md for reuse terms. Third-party resources retain their own terms; no new license is assigned to them. Original papers, complete workbooks, model weights and reference proteomes are not redistributed.

See `DATA_DICTIONARY.md` and `data_dictionary.json` for schemas, units, missing values and table provenance.
