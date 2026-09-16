# GloGlo Research

Code, data and reports accompanying research findings shared by GloGlo Labs.
Each study is a self-contained, versioned package with its own methods, inputs,
source attribution and reproducibility command.

## Studies

| Study                                                          | Date       | Materials                                                                                                |
| -------------------------------------------------------------- | ---------- | -------------------------------------------------------------------------------------------------------- |
| [CHGA-INS leader: L9C and T17C](studies/2026/chga-ins-leader/) | 2026-09-15 | Full working paper, research note, sequences, cached predictions, analysis scripts and reference results |

## Reproduce

Clone this repository and run from its root with Python 3.10 or later:

```sh
python3 scripts/verify.py
```

The current study uses only the Python standard library. Verification checks
input hashes, runs the analysis in isolated mode, and compares every regenerated
CSV with the published reference results. No network requests or predictor
installation are needed. New outputs go into a temporary directory.

## Organization

- `studies/<year>/<study-slug>/`: stable home for a paper, note or related analysis.
- Each study's `README.md`: question, scope, report links, methods and run instructions.
- Each study's `study.json`: identity, reproducibility entry point and reference outputs.
- Study-local data and provenance: inputs, retrieval sources, checksums and schemas.
- `scripts/`: shared publication validation, without study-specific scientific logic.

Keep studies independent. Add shared analysis code only when multiple studies
use the same scientific contract. Large datasets belong in a durable public
archive or release asset, linked from the study with a retrieval command,
checksum, version and reuse terms. Small inputs needed for offline reproduction
can live alongside the study.

See [CONTRIBUTING.md](CONTRIBUTING.md) for adding a study and publishing a release.
Original code is MIT-licensed; original reports and generated research data are
CC BY 4.0. See [LICENSE.md](LICENSE.md) for the boundaries and third-party terms.

## Citation

Cite the specific study, its authors or organization, and the exact Git commit or
study release tag used. Repository-level citation metadata is in [CITATION.cff](CITATION.cff).
Corrections and reproducibility questions can be raised in GitHub Issues.
