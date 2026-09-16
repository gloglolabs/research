# Data dictionary

## Common conventions

- CSV files are UTF-8 with a header. Boolean outputs are `True`/`False`. Empty fields mean unmeasured, unscanned, unknown or not applicable as defined below; they are never numeric zero by default.
- All **curated** coordinates are one-based and inclusive. Raw MHCflurry `pos` is zero-based; `reproduce.py` converts it explicitly. CHGA leader positions and whole reconstructed-precursor positions coincide through 18. INS coordinate 25 maps to precursor 19; an INS coordinate maps to graft coordinate `INS_position−6`.
- The precursor is 104 residues. The final NetMHCpan variant request covers its first 29 residues. Beyond that context, the clean table uses independently cached canonical-precursor scores only for exactly identical peptides. Peptide lengths are 8–11; seven HLA alleles are listed in the manifests.
- `signal_or_junction=True` means peptide start≤18: the window contains at least one leader residue. Leader-only windows end≤18; junction windows start≤18<end. An unchanged INS core starts≥19. Overlapping windows and multiple allele assignments are counted as separate window/HLA pairs.
- Affinity fields are predicted nanomolar values; smaller means stronger predicted binding. EL/presentation percentile ranks are percentages, not probabilities; smaller means stronger model output. All count thresholds use strict `<`, including API-rounded NetMHCpan values.

## `results/final_peptide_HLA_predictions.csv`

34,762 unique rows keyed by `(construct,HLA,start_1based,end_1based)`. This is the principal clean per-peptide table.

| Field | Definition |
|---|---|
| `construct` | ID matching `constructs.fasta` and the 13-row construct ledger |
| `HLA` | One explicitly scored allele, not an experimentally assigned restriction |
| `start_1based`, `end_1based` | Inclusive coordinates in the 104-aa reconstructed precursor |
| `peptide` | Exact uppercase amino-acid sequence |
| `n_flank_5aa`, `c_flank_5aa` | Up to five actual precursor residues immediately before/after the peptide; shorter at a terminus |
| `signal_or_junction` | Whether the peptide contains any leader residue |
| `net_score_origin` | `direct_variant_context` for the original 29-residue request; `exact_peptide_reused_canonical_full_scan` for a verified identical downstream peptide |
| `net_EL_percentile` | NetMHCpan 4.1 EL percentile, with complete coverage after identity-verified reuse |
| `net_BA_nM`, `net_BA_percentile` | NetMHCpan 4.1 predicted affinity and binding percentile; origin follows `net_score_origin` |
| `mhc_affinity_nM`, `mhc_affinity_percentile` | MHCflurry affinity output and its calibrated percentile; affinity is peptide/HLA-specific |
| `mhc_presentation_percentile_with_flanks` | Presentation prediction with actual five-residue flanks |
| `mhc_presentation_percentile_without_flanks` | Presentation prediction using the no-flank model |
| `mhc_processing_with_flanks`, `mhc_processing_without_flanks` | Model processing scores, dimensionless; not experimental cleavage efficiency |

## Count and objective outputs

`results/prediction_counts_by_allele.csv` and `results/objectives42.csv`:

| Field | Definition |
|---|---|
| `construct` | Scored construct |
| `model` | `NetMHC` or `MHCflurry`; EL output is used for NetMHC counts |
| `context` | `with_flanks` or `without_flanks` for the MHCflurry sensitivity setting. NetMHC values repeat across settings because it is peptide-only |
| `scope` | `signal_or_junction` (start ≤ 18), `full_precursor`, `unchanged_INS_core` (start ≥ 19), or `flank_affected` (start ≤ 22) |
| `HLA` | Allele counted separately |
| `cutoff_percentile` | 0.5, 1 or 2, with strict less-than |
| `baseline_count`, `variant_count` | Number of qualifying overlapping peptide/HLA windows within that scope |
| `delta` | Variant minus baseline; negative means fewer qualifying predictions |

`results/criteria42_summary.csv` contains `better`, `unchanged`, `worse`: the numbers of 42 fixed criteria with negative, zero or positive delta. `results/prediction_totals.csv` sums allele counts by construct/model/context/scope, with `count_lt_0.5`, `count_lt_1`, `count_lt_2`. Full-precursor NetMHC totals include exact-identity reuse from the independently scanned full canonical precursor. `results/expanded_objectives42.csv` and `results/expanded_criteria42_summary.csv` use the same schema for the 616-pair `flank_affected` scope.

`results/{leader_3_9,position17,double}_netmhcpan_summary.csv` adds `stage`, `substitution`, `scope`, `window_HLA_pairs`, `EL_lt_*` and `BA_lt_500nM`. The `position17_containing` scope counts only windows containing residue 17; it must not be confused with the complete signal/junction region.

`results/single_mhcflurry_summary.csv` summarizes both single-residue stages by construct, flank context, scope and strict presentation cutoff. The position-3–9 raw cache contains 67,536 rows per mode, restricted to complete signal/junction windows after scoring a 33-residue prefix. The position-17 cache contains 53,480 rows per mode from full sequences.

`results/position17_objectives42.csv` records 42 criteria for each of 20 choices in both flank contexts, separately for `position17_containing` (266 pairs) and `position17_core_or_flank` (546 pairs). The latter includes any peptide whose core or five-residue flanks contain residue 17. Its summary reports improved, unchanged and worse criterion counts.

## Construct sequences and SignalP

`constructs.fasta` stores 13 full precursor sequences. `tables/constructs.csv` adds:

- `edited_positions_1based`: semicolon-delimited leader positions; blank for baseline.
- `leader_sequence`: 18 amino acids.
- `precursor_length`: 104.
- `predicted_cleavage_after`: last leader coordinate retained before model cleavage.
- `cleavage_probability`: SignalP output for that predicted boundary; not a measured functional probability.
- `native_INS25_110_unchanged`: deterministic sequence identity.
- `sequence_status`: canonical-design status.

`cached/signalp/input.fasta` and processed output contain 16 entries: the 13 final HLA-scan constructs plus native INS, T17D and T17P controls. This is a different denominator from the 13-member final prediction panel. `cached/signalp_position17/` contains the initial 21-sequence batch: native INS and all 20 position-17 choices. Batch identifiers and probabilities are retained separately.

## Donor data

`tables/donor_observations.csv` is the nine-row INS/CHGA excerpt from Nanaware 2025:

| Field group | Definition |
|---|---|
| `excel_row`, `protein`, `gene`, `peptide`, `length` | Original Table S5 row and molecular identification |
| `ion_intensity` | Original recorded peptide ion intensity in arbitrary instrument units; no cross-donor normalization or pooling |
| `donor` | Published donor label, including original spacing |
| `predicted_hla_as_published` | Source computational assignment; multi-allele strings are preserved |
| `start`, `end` | Coordinates in the source protein, not the reconstructed precursor |
| `footprint` | Native INS signal, native INS boundary, retained downstream INS, native CHGA signal, or downstream CHGA |
| `sequence_in_canonical_CHGA1_18_INS25_110` | Exact sequence inclusion, not predicted presentation after editing |
| `S4_exact_blank_intensity` | Matched source blank intensity. `0` is explicitly zero; empty is unknown |
| `same_donor_blank_overlap_at_least_6aa` | Precomputed source-table blank peptide(s) sharing at least six residues; empty means none identified by that original extraction |

`tables/donor_depth_and_projection.csv` records all six preparations. `all_curated_mhci_species` is the number of source-curated rows in that donor; zero INS observations are retained. Source fractions describe detected distinct species only and are not protection estimates.

`results/donor_peptide_matrix.csv` has 24 rows for all six donors × four INS peptides. `detected_in_curated_table=False` means absent from the curated table, not biological absence. Intensity and HLA are blank for those rows. `donor_label` removes spaces only for consistent plotting; `donor_as_published` preserves provenance.

`results/donor_sensitivity.csv` and `results/donor_leave_one_out.csv` count observations, distinct sequences, INS-positive donors and donors with retained sequences. The “explicit zero blank” rule excludes unknown blanks, whereas “drop positive blank, keep unknown” does not.

## Inflammation excerpts

Files beginning `tables/inflammation_` retain source worksheet addresses. `mean_basal_intensity`/`mean_ifna_intensity` are main-study condition means across four biological replicates; `independent_median_basal`/`independent_median_ifna` are separate-experiment summaries across three biological replicates. `ifna_fold_change` is within-peptide IFN/basal. Identification counts are MS detection counts, not peptide-level intensity replicates. No B10–18 replicate variance is inferred from these summaries.

The Fig. 8c run-intensity and consistency files document the source SP15–24 discrepancy. Three non-INS peptides in the larger source-coordinate audit had ambiguous alternative-protein coordinate mapping; those do not affect the 29-row canonical INS projection supplied here.

`results/reference_junction_comparison.csv` contains three native or reconstructed junction peptides × seven alleles. It uses the same rank, affinity and coordinate conventions as the principal prediction table. `metadata/reference_mhcflurry_sequences.json` records five initial reference/context sequences, a separate panel from the final 13-member FASTA.

## Raw predictor files and metadata

`cached/*/netmhcpan_{el,ba}.tsv` preserves provider fields: `seq_num` indexes that stage's input manifest; `start/end` are one-based; `core/icore` are model alignment fields; EL `score` is a raw model score; BA `ic50` is predicted nanomolar affinity; `percentile_rank` belongs to the indicated model output.

Raw MHCflurry CSVs preserve `sequence_name`, zero-based `pos`, `peptide`, flanks, `sample_name`, `affinity`, `best_allele`, `affinity_percentile`, `processing_score`, `presentation_score`, `presentation_percentile`. Each sample contains a single HLA allele, so `best_allele` is not a mixed-genotype inference.

JSON manifests retain method names, sequence IDs, runtime/model versions, timestamps and source hashes. `data_dictionary.json` inventories each bundled CSV/TSV header and row count; this document supplies the scientific semantics.
