# L9C and T17C in the CHGA leader of engineered insulin

GloGlo Labs | Computational follow-up to Kobaisi et al. (2026) | 15 September 2026 | v1.1

**Question:** can two substitutions reduce HLA-I presentation from the replacement CHGA leader while preserving insulin function?

Kobaisi et al. showed that a CHGA leader can support functional insulin production while preventing recognition by the tested 1E6 TCR. We examined peptides within that replacement leader and across its junction with proinsulin. The reconstructed junction yields the HLA-C\*03:04 candidate `VTAFVNQHL` (NetMHCpan EL rank 0.25%; MHCflurry presentation rank 0.1112%). L9C targets leader peptides; T17C changes the junction. Presentation of these candidates remains unmeasured.

## Defined sequence comparison

Parental CHGA leader: `MRSAAVLALLLCAGQVTA` + INS25-110\
L9C + T17C candidate: `MRSAAVLACLLCAGQVCA` + INS25-110

Coordinates refer to the CHGA leader. The 104-residue precursor is a canonical CHGA1-18 + INS25-110 reconstruction; the exact synthesized construct was unavailable. Confirm that sequence before making constructs or interpreting junction predictions.

## Predicted presentation counts

Counts below each strict percentile cutoff, parental to double variant:

| Predictor and context    | <0.5%  | <1%     | <2%      | Reduction at <2% |
| ------------------------ | ------ | ------- | -------- | ---------------- |
| NetMHCpan EL             | 4 to 0 | 8 to 2  | 17 to 8  | 52.9%            |
| MHCflurry with flanks    | 9 to 3 | 19 to 6 | 27 to 14 | 48.1%            |
| MHCflurry without flanks | 4 to 1 | 12 to 5 | 21 to 12 | 42.9%            |

Denominator: 504 overlapping leader/junction peptide-window/HLA pairs (72 windows of 8-11 residues x seven alleles). These are prediction counts, not immune-escape percentages. Across the full precursor, <2% counts decrease 59 to 50 in NetMHCpan (15.3%) and 72 to 59 in MHCflurry with flanks (18.1%); downstream insulin sequences remain unchanged.

**L9C single control:** It also passes the reported no-increase criterion. Adding T17C lowers <2% counts from 11 to 8 in NetMHCpan and 20 to 14 in MHCflurry with flanks. That incremental reduction is the reason to test the second substitution.

**Selection scope:** 152 single substitutions and eight shortlisted doubles were screened. L9C + T17C was the only tested double with no increase across 42 allele/model/cutoff comparisons under both flank settings. Selection and this audit were exploratory. At a strict 3% MHCflurry cutoff without flanks, A\*03:01 increases from zero to one candidate (`QVCAFVNQH`, rank 2.9106%; aligned parental peptide 4.7322%).

<!-- pagebreak -->

## What still needs testing

SignalP 6.0 predicts the same 18/19 cleavage boundary and INS25-110 product (probability 0.979955 parental; 0.979581 double). It does not test the added cysteines' effects on processing, secretion, ER stress or insulin activity.

Three HLA-C*03:04 pairs remain below the MHCflurry 0.5% cutoff: `RSAAVLACL`, `SAAVLACLL` and `CAFVNQHL`. The A*02:01 peptide `VLACLLCAGQV` has stronger predicted NetMHCpan binding than its parental counterpart (654.57 to 202.39 nM), although EL rank remains 5.7%.

## Proposed experiment and decision

**Constructs:** parental CHGA-INS, L9C, T17C and L9C + T17C, using the confirmed fusion sequence in the same cellular background and expression system. The single variants separate each substitution's contribution.

**Presentation:** confirm the target cells' HLA genotype. Compare HLA-I-associated leader/junction peptides, including parental `VTAFVNQHL`, variant `VCAFVNQHL` and the residual candidates above, using targeted mass spectrometry with appropriate peptide standards. Include the A\*02:01 affinity-tradeoff peptide where that allele is expressed. Normal cleavage cuts through junction peptides, so their cellular production must be established.

**Function:** measure proinsulin/insulin ratio, cellular insulin content, glucose-stimulated secretion, secreted insulin activity and ER stress, with viable-cell and expression measurements to interpret changes.

**Recognition:** test construct-expressing cells with HLA-matched, peptide-specific T cells where available; use peptide pulsing to confirm specificity. Retain the original 1E6 assay as a control; it cannot establish reduced recognition of the new leader/junction candidates.

**Decision:** advance the double only if it reduces measured presentation beyond L9C alone while preserving processing and function, with no increased recognition in the tested T-cell panel. If L9C preserves function and T17C adds no measurable benefit, favor L9C.

## Methods and materials

HLA panel: A*02:01, A*03:01, A*24:02, B*40:01, B*49:01, C*03:04 and C\*07:01. Predictors: [NetMHCpan 4.1](https://doi.org/10.1093/nar/gkaa379), [MHCflurry 2.2.1 / model release 2.2.0](https://doi.org/10.1016/j.cels.2020.06.010) and [SignalP 6.0](https://doi.org/10.1038/s41587-021-01156-3). MHCflurry was evaluated with available five-residue precursor flanks and without flanks.

[Full paper: methods, figures and sensitivity analyses](https://github.com/gloglolabs/research/blob/main/studies/2026/chga-ins-leader/chga-ins-leader-variants.pdf)\
[Code, sequences, 34,762 prediction rows and provenance](https://github.com/gloglolabs/research/tree/main/studies/2026/chga-ins-leader)

From the repository root, `python3 scripts/verify.py` checks input hashes and reproduces the reported results from cached predictions, without rerunning model inference.

Source study: [Kobaisi et al., Cell Reports 45, 117941 (2026)](https://doi.org/10.1016/j.celrep.2026.117941). Computational analysis and drafting used Codex agents; scripts and predictor outputs are available for review. No new laboratory experiments were performed.
