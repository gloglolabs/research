# Contributing research

## Add a study

1. Create `studies/<year>/<descriptive-slug>/` and add it to the root study table.
2. Write a README with the research question, evidence status, scope, report links,
   exact run command, runtime requirements and expected outputs.
3. Include scripts, compact inputs and reference results sufficient to reproduce
   the reported findings. Keep original inputs immutable; write new results to a
   separate output directory.
4. Record source URLs, versions, retrieval dates, checksums, units and missing-value
   meanings. Attribute third-party data and state their applicable reuse terms.
   For externally hosted inputs, provide exact retrieval instructions and hashes.
5. Add `study.json` using the existing study as the schema example. Its
   `reproduction` paths are relative to that file. The entry point must accept
   `--output`, write a `validation.json` with `status: PASS`, and support offline
   verification from a fresh checkout. Declare the expected CSV filenames.
6. Run `python3 scripts/verify.py`. CI verifies all registered studies. If a future
   study needs additional dependencies, declare and pin them locally and extend
   the verification workflow explicitly.

Preserve negative results and analysis-selection history. Distinguish model
predictions, retrospective analyses and experimental observations. Corrections
must explain which result changed and why in the study's changelog.

## Publish a release

Review the exact files being published for private correspondence, credentials,
participant information and redistribution rights. Exclude original third-party
papers, proprietary tools and model weights unless their terms permit inclusion.
Use study-scoped tags such as `chga-ins-leader-v1.0.0` and attach large public data
only with a manifest, checksums and source attribution. Preserve previous tags.
Record substantive corrections in a new release, keeping the old result traceable.

Contributions retain the code/content licensing split in LICENSE.md. Identify
exceptions next to the relevant files and never relicense third-party material.
