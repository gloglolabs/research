# Provenance and analysis boundaries

## Study identity

This is an original computational sensitivity analysis prepared 15 September 2026. It uses a canonical human CHGA signal sequence joined to native proinsulin. The analysis concerns that defined reconstruction; it does not claim access to the exact experimentally synthesized construct.

The baseline is CHGA1–18 + INS25–110, 104 residues. The final 13-member panel comprises baseline, L9C, L9G, T17C, T17G, and eight combinations of position 9 C/G with position 17 C/G/D/P. The preceding scans covered every non-native residue at CHGA positions 3–9 (133 alternatives) and every alternative at position 17 (19).

## Predictions generated for this analysis

| Component                       | Recorded implementation                                                        | Included provenance                                                                                             |
| ------------------------------- | ------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------- |
| NetMHCpan affinity              | Official IEDB API, `netmhcpan_ba-4.1`                                          | Cached TSVs, successful request payloads, input sequence manifests                                              |
| NetMHCpan eluted-ligand score   | Official IEDB API, `netmhcpan_el-4.1`                                          | Cached TSVs, successful request payloads, input sequence manifests                                              |
| MHCflurry affinity/presentation | Code 2.2.1; model distribution 2.2.0, presentation weights dated 20200611; CPU | Complete final-panel outputs with/without flanks; runtime versions, model-file hashes, original producer script |
| Signal peptide cleavage         | SignalP 6.0, eukaryote, fast mode                                              | Full submitted FASTA, returned output JSON/text and processed FASTA, retrieval timestamp, source hashes         |

NetMHCpan requests used peptide lengths 8,9,10,11 and seven alleles: A*02:01,A*03:01,A*24:02,B*40:01,B*49:01,C*03:04,C\*07:01. The final NetMHCpan variant scan used 29-residue contexts covering every signal and junction window. Complete clean grids reuse predictions from the independently scanned full canonical precursor beyond that context only after verifying exact peptide identity and absence of mutated residues. The 616-pair sensitivity analysis therefore includes every peptide affected through five-residue flanks. Full 104-residue MHCflurry sequences included those complete flanks. The position-3–9 MHCflurry scan used a 33-residue prefix and retained all 504 signal/junction pairs, preserving complete five-residue flanks. The position-17 scan used full 104-residue sequences and the expanded single-position analysis includes 546 affected pairs. These are model predictions of presentation or binding, not new experimental measurements.

The 42 criteria are 7 alleles × 2 presentation models × 3 strict percentile cutoffs, 0.5%, 1%, 2%. Full-flank and no-flank comparisons are separate sensitivity analyses, not independent experiments. Their percentile calibrations are not averaged. Criteria and nested windows are not statistically independent trials. Affinity predictions are recorded separately and are not included in these 42 presentation-count criteria.

Model references: [NetMHCpan 4.1](https://doi.org/10.1093/nar/gkaa379), [official IEDB API](https://tools.iedb.org/main/tools-api/), [MHCflurry 2.0](https://doi.org/10.1016/j.cels.2020.06.010), [SignalP 6.0](https://doi.org/10.1038/s41587-021-01156-3). Provider software, weights and services retain their own terms; this bundle does not redistribute model weights or a licensed predictor executable.

## Public experimental sources

### Nanaware 2025 donor ligand projection

[Nanaware et al., The antigen presentation landscape of cytokine-stressed human pancreatic islets](https://pmc.ncbi.nlm.nih.gov/articles/PMC12624573/), PMID 40684438. Relevant source files are Table S4 (`NIHMS2107237-supplement-4.xlsx`) and Table S5 (`NIHMS2107237-supplement-5.xlsx`). Their full checksums and retrieval URLs are in `provenance/retrieval_manifest.json`.

`tables/donor_observations.csv` is a nine-row INS/CHGA excerpt derived from the authors' curated class-I table, joined to its matching raw-hormone/blank entry. Six rows concern INS, corresponding to four unique sequences. Source row numbers, donor labels, recorded intensities and published predicted HLA assignments are preserved. `tables/donor_depth_and_projection.csv` records the six-preparation denominator and detection depth. The complete 3,766-row class-I table and original workbooks are not bundled.

The standard-library runner recomputes sequence mapping and sensitivity from those attributed excerpts. It does not replay the authors' mass-spectrometry database search. The original workbook-extraction code is retained in `archived_scripts/analyze_peptidome.py`; it requires separately retrieved source workbooks and OpenPyXL. Its original relative-layout requirements are documented in the script itself.

### Carré 2025 inflammation summaries

[Carré et al., Interferon-α promotes neo-antigen formation and preferential HLA-B-restricted antigen presentation in pancreatic β-cells](https://www.nature.com/articles/s41467-025-55908-9), PMID 39824805. The included target-level excerpts derive from Supplementary Data 2 and Source Data worksheets Fig. 4d/e, Fig. S3c and Fig. 8c.

Main intensities are the source-reported condition means across four biological replicates. The independent experiment provides peptide-level median summaries across a three-replicate proteasome experiment on another mass-spectrometry platform. The main and independent B10–18 IFNα/basal ratios reproduce 11.25281 and 1.43803. The source provides no per-replicate B10–18 intensities in these worksheets; no B10 confidence interval is computed. A discrepant SP15–24 summary between Fig. 4 and Fig. 8 is preserved in the supplied consistency table.

The original workbooks are excluded; attributed target-level excerpts, source-cell addresses, exact hashes and retrieval URLs are included.

### Functional evidence used to interpret the model

[Whalley et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7058665/) supplies independent published human β-cell killing evidence for INS B10–18/HLVEALYLV and the InsB4/InsB6 TCR. This is prior experimental evidence, not a new experiment in this computational study. [Kobaisi et al., Cell Reports, 3 September 2026](https://doi.org/10.1016/j.celrep.2026.117941) motivates the canonical CHGA signal replacement. No original article/PDF is redistributed.

## Sequence and database references

Canonical protein references are UniProt INS/P01308 and CHGA/P10645; entry and sequence versions are recorded in `metadata/reference_proteins.json`.

The exact/I-L-equivalent search used 169,634 sequence records and 82,674,231 residues from UniProt human reference proteome UP000005640 plus additional sequences, release 2026_03. Published MD5 and locally recorded SHA256 checksums were verified. The included `proteome_exact_matches.json` records all query results, matched positive controls and source hashes. Full reference FASTA files are excluded. The optional standalone `audit_reference_matches.py --directory PATH` rescans separately retrieved matching files and validates the complete result.

The IEDB exact query was live and unfiltered at the recorded retrieval time. Query URL, compact response and timestamp are included. Database non-detection establishes absence from the queried records, not a novel TCR specificity. Retrieval URLs using `current_release` are mutable: the stored checksums identify the required snapshot, and mismatching later downloads must not be silently substituted.

## File lineage

`provenance/bundle_origins.json` links copied files to their original research-relative paths and source hashes. `checksums.sha256` identifies the included executable code, cached inputs and compact metadata. `results/` contains regenerated derived outputs. The original producer scripts are preserved for inspection, while `reproduce.py` supplies the portable cached-data calculation.

See the repository LICENSE.md for reuse terms. Third-party sources retain their original terms; no blanket license is applied to them.
