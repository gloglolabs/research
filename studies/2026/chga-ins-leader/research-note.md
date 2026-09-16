# L9C and T17C in the CHGA leader of engineered insulin

Computational analysis of Kobaisi et al. (2026) | 15 September 2026

**Candidate:** introduce L9C + T17C into the CHGA-derived leader of the engineered INS construct. The two substitutions lower predicted HLA-I presentation counts while preserving the predicted cleavage site and the proinsulin sequence.

## Sequence

Baseline: `MRSAAVLALLLCAGQVTA` + INS25-110  
Candidate: `MRSAAVLACLLCAGQVCA` + INS25-110

Coordinates refer to the CHGA leader. Calculations assume CHGA1-18 + INS25-110, a 104-residue precursor reconstructed from canonical sequences. The exact synthesized construct in Kobaisi et al. was unavailable; sequence agreement requires confirmation.

## Quantitative results

Counts below each strict percentile cutoff, baseline → candidate:

| Predictor and context    | <0.5% |    <1% |     <2% | Reduction at <2% |
| ------------------------ | ----: | -----: | ------: | ---------------: |
| NetMHCpan EL             | 4 → 0 |  8 → 2 |  17 → 8 |            52.9% |
| MHCflurry with flanks    | 9 → 3 | 19 → 6 | 27 → 14 |            48.1% |
| MHCflurry without flanks | 4 → 1 | 12 → 5 | 21 → 12 |            42.9% |

Denominator: 504 leader/junction peptide-window/HLA pairs, comprising 72 overlapping 8-11mer windows starting at precursor positions 1-18 × seven HLA alleles. These are predicted-presentation counts, not measured immune-escape percentages.

**Full precursor:** at <2%, NetMHCpan counts decrease 59 → 50 (**15.3%**) and MHCflurry with flanks 72 → 59 (**18.1%**), across 2,674 pairs. Downstream INS peptide sequences remain unchanged.

**Selection:** 152 single substitutions and eight selected double variants were screened. L9C + T17C was the only tested double with no count increase in any of 42 comparisons (seven alleles × two models × three cutoffs) under both flank settings: 16 decreased / 26 unchanged with flanks; 15 / 27 without.

**Cleavage:** SignalP 6.0 predicts the same 18/19 boundary and INS25-110 product. Cleavage probabilities are 0.979955 → 0.979581. The added cysteines have not been assessed for folding, secretion or insulin activity.

## Residual predictions

Three HLA-C*03:04 pairs remain below the MHCflurry 0.5% cutoff: `RSAAVLACL`, `SAAVLACLL` and `CAFVNQHL`. The A*02:01 peptide `VLACLLCAGQV` has stronger predicted NetMHCpan binding than its baseline counterpart, 654.57 → 202.39 nM, while EL rank remains 5.7%.

**Suggested comparison:** parental versus L9C + T17C constructs, measuring HLA-I ligand presentation alongside proinsulin processing, secretion and insulin activity.

## Methods and accompanying files

HLA panel: A*02:01, A*03:01, A*24:02, B*40:01, B*49:01, C*03:04, C\*07:01. Predictors: [NetMHCpan 4.1](https://doi.org/10.1093/nar/gkaa379), [MHCflurry 2.2.1 with 2.2.0 models](https://doi.org/10.1016/j.cels.2020.06.010), [SignalP 6.0](https://doi.org/10.1038/s41587-021-01156-3). Flank comparisons use available five-residue flanks or none. Accompanying data contain 34,762 prediction rows for 13 constructs, sequences and cached outputs; `python3 reproducibility/reproduce.py` reproduces the reported counts offline.

Source study: [Kobaisi et al., Cell Reports 45, 117941 (2026)](https://doi.org/10.1016/j.celrep.2026.117941).

Computational analysis and drafting used Codex agents; predictor outputs and analysis scripts are supplied for review.
