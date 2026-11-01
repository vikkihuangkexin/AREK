# Public verification inputs and software templates

`real_fastq/PRJDB10075/` contains public subsets from runs DRR234531–DRR234534 in project PRJDB10075, retaining the first 2,000 read pairs per run. Original reads were recovered from imported FASTQ QZA archives without synthetic or altered sequences. Accessions, extraction details, and SHA256 values are recorded in [provenance.json](provenance.json).

`matrix/counts.tsv` and `matrix/taxonomy.tsv` are feature/taxonomy verification fixtures for these subsets. They contain no full study cohort or verified drinking/disease/clinical labels. `real_samples.tsv` explicitly maps PE reads; `real_se_samples.tsv` explicitly selects R1 from the same data to verify the SE entry.

The RUN_A/RUN_B IDs, groups, ages, and sexes in `metadata_template.tsv` illustrate field format and are not participant records. Shotgun manifest paths are placeholders to be supplied by the user; no real shotgun dataset is included.

These inputs support installation and software-entry verification. They cannot support alcohol-related biological conclusions. Full study datasets and large references are external.
