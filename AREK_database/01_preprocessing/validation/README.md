# Validation entries and evidence scope

Validation programs live here and write new evidence into the specified work directories. Historical logs, snapshots, and full validation outputs remain in the original source archive.

## Reproducible checks

Run from `AREK_database/01_preprocessing/`:

```bash
# Current manifest, syntax, internal paths, and import checks.
python validation/audit_package.py .

# Repository relocation, route dry-runs, and independent matrix numerics.
# Does not require QIIME2, R, or reference databases.
python validation/check_repository.py --report /tmp/arek_repository_checks.json

# Real matrix, fault, and statistical-design tests; requires analysis and R dependencies.
python validation/validate_release.py --work-root /tmp/arek_validation_matrix --matrix-only

# Also tests real PE and explicit R1-SE; requires amplicon environment and classifier.
python validation/validate_release.py --work-root /tmp/arek_validation_full \
  --classifier "$AREK_CLASSIFIER"

# Additional data-contract and statistical-design checks.
python validation/check_contract_transitions.py /tmp/arek_contract
python validation/check_extra_statistics.py /tmp/arek_statistics

# Linux lock test; terminates only processes created by the test itself.
python validation/check_orphan_lock.py /tmp/arek_lock
```

`validate_release` writes evidence under `--work-root/evidence/`. Contract and extra-statistics checks also use their work directories' `evidence/`. The lock test writes `orphan_lock_test.json` under its work root. Historical source-package evidence is not overwritten.

## Historical source-package validation

The source archive records real PRJDB10075 PE processing, an explicitly selected R1-SE chain, real matrices and relative inputs, resume/fault recovery, independent numerics, and statistical-design checks. Consult its `validation/VALIDATION_zh.md` for detailed provenance and unverified scope. Historical execution is not current-repository rerun evidence. The repository's [public examples](../examples/README.md) preserve accession/extraction/hash records without cohort phenotypes.

## Repository checks and English-language update

[repository_checks_20261031.json](repository_checks_20261031.json) records the original repository-preparation checks. A fresh [english_update_checks.json](english_update_checks.json) records checks after the English-language update. These checks cover entries, internal paths, configuration, discovery, matrix calculations, and source preservation. Additional language, link, configuration-equivalence, and source-comparison results are recorded in [english_content_audit.json](english_content_audit.json).

QIIME2, Rscript, and large reference databases are unavailable in the editing environment. The full FASTQ and R statistical chains were not rerun here. Source translation is limited to comments and messages; metric implementations and input data are preserved.
