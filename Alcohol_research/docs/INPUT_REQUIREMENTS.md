# Inputs and environments

## Read-only authorities and aliases

| Alias | Required layout | Purpose |
|---|---|---|
| `M1` | `reproducibility/inputs/`, comparison folders, `evidence_supplemented/` | Frozen M1 profiles, memberships, families, completed results |
| `CLIN` | `reproducibility/inputs/metadata.tsv.gz`, `harmonised_clinical_values.tsv.gz`, `evidence/`, `interaction_models/`, `presentation/source/` | Accepted clinical evidence, original clinical results, immutable presentation |
| `V3` | `provenance/source_inputs/<analysis_set>/TSS.tsv.gz`, `CLR_unadjusted.tsv.gz`, `provenance/taxon_crosswalk.tsv.gz` | Source-to-canonical profile reconstruction and original full-composition alpha |
| `V2` | `inputs_v1/<body>/layers/<analysis_set>/TSS.tsv.gz`, `functional/faprotax/`, `<body>/faprotax/` | Official cached annotation states and fresh source-scale function reconstruction |
| `DB` | `faprotax.R`, `FAPROTAX_ref_data/` | Optional original mapper re-execution; not performed in validation |

`configs/external_inputs.tsv` pins the two ZIP identities and key provenance assets. `history/source_lineage.tsv` maps the current implementations to recovered originals. `history/recovered_code_index.tsv` indexes additional historical scripts without pretending they are supported active entry points. Original ZIP manifests remain in their original archives.

## Table conventions and required fields

All analysis tables are tab-separated, optionally gzip compressed. Specimen-by-feature matrices use the first column as the exact unique specimen index; feature columns are **canonical IDs**, not display genus labels. `NA` is missing. Do not let spreadsheet tools round numeric columns or convert IDs.

| Input | Required schema / constraints |
|---|---|
| `metadata_4603.tsv.gz` | Exact specimen index; `Study`, `analysis_set`, `assay_analysis`, `material_analysis`, `technical_batch`, `working_cluster` and frozen clinical/support metadata. Preserve all columns. |
| `comparison_membership_4603.tsv.gz` | Same ordered index; five exact comparison keys; values `case`, `reference`, `excluded`; `clinical_complete_age_sex_BMI`. |
| `panels/<panel>/TSS.tsv.gz` and `CLR_unadjusted.tsv.gz` | Identical index/columns and structural-missingness mask; TSS fractions, not percentages. |
| `panels/<panel>/family.tsv` | `taxon_key`, `selected`, `genus_name`, `taxon_status` and original family metadata; unique keys and complete panel coordinates. |
| `diversity/<panel>/alpha_per_specimen.tsv.gz` | Exact specimen index; `Shannon`, `Simpson`, `Observed_genera`; full-source measures originate from individual source compositions. |
| `functions/<body>/contribution.tsv.gz`, `family.tsv` | Original-scale trait contributions and fixed trait family. Coverage branches additionally use `global_coverage.tsv.gz` / `known_mass.tsv.gz`. |
| Clinical `metadata.tsv.gz` | Frozen technical metadata, `body`, `analysis_group`, `working_cluster`, `explicit_disease_reference`; exact index. |
| `harmonised_clinical_values.tsv.gz` | Exact index; original variable keys, levels and units. Original string levels such as numeric-coded categories must retain their model meaning. |
| `individual_value_evidence.tsv.gz` | `variable`, `sample_id`, `status`, `value` and original provenance columns; only accepted exact-ID evidence is promoted. |
| Source canonical crosswalk | `analysis_set`, `original_feature_id`, `taxonomy_lineage`, `resolution_status`, `genus_name`, `taxon_key`, `excluded_from_genus_inference`. |
| Optional rank aggregation | Matrix plus taxonomy with `feature_id` and the requested source-supported rank column. No lineage backfilling or species inference. |
| Candidate contributions | `candidate_taxon_trait_mapping.tsv.gz` with `analysis_set`, `original_feature_id`, `taxon_key`, `trait`, `accepted_state`; original source TSS tables; frozen expected per-ID contribution table. |

Use original files, not a hand-recreated subset, for a full reproduction. The original Biomni archive's AC/demo contains the frozen fields/files for audit/profile/membership/main-M1 checks. It is external to this code release. Some modules compare newly computed values against completed result tables; those tables must remain available under their original relative paths.

## Raw reads and source adapters

The original Biomni archive retains history/raw_acquisition_manifest.tsv.gz with acquisition URLs/bytes/MD5 and history/source_file_manifest.tsv with original source objects. These large inventories are not bundled in the GitHub increment. URL-only entries are **not downloaded reads**. Neither a downstream profile nor a completed result proves raw-read processing has been reproduced.

The included adapters are recovered substantive source-specific programs:

- `amplicon_v5`: mate-aware QIIME 2 / VSEARCH, 97% OTUs, biological-mate selection or paired-read joining as defined in the source, median read-length gate 100, QC read threshold 3,000, SILVA 138-99 classifier. No new primer-trimming step is added.
- `shotgun_v2`: fastp paired-end adapter detection, minimum length 50; hg38 Bowtie2 very-sensitive host depletion; samtools `-f 12 -F 256`; paired Kraken2 and Bracken S/G/F/P ranks, read length 150 and threshold 10.

Other historical source routes are indexed but not all ported. Do not apply one adapter indiscriminately to every analysis set.

Copy and edit the example JSON. Paths resolve relative to the **working directory**, not the JSON directory. Input paths must be original assets; output paths must be inside the exact module node. For `--output-root run/raw`, use:

| Module | Amplicon `RESULT_ROOT` | Shotgun `OUT_ROOT` / `MISS_ROOT` |
|---|---|---|
| 03 | `run/raw/03/results` | `run/raw/03/results`, `run/raw/03/missing` |
| 04 | `run/raw/04/results` | Not applicable |
| 05 | Not applicable | `run/raw/05/results`, `run/raw/05/missing` |

The shipped examples are templates, not ready-to-run configurations. Set the original project, `RAW_ROOT`, classifier or host/Kraken indexes and every destination. Module 03 stops after the adapter's layout/QC stage; modules 04/05 execute their entire source adapter, including its own QC. Module 03's output is not automatically used by modules 04/05.

## Environments and missing pins

- Python: exact installed versions in `requirements-python.lock.txt`; validated Python 3.11.13.
- R: validated R 4.4.2, vegan 2.7.6, cluster 2.1.8.3. Install those exact versions using an appropriate R package repository/archive. `requirements-R.lock.tsv` is an environment record, not an automatic installer.
- Shell programs: `requirements-shell.tsv` distinguishes required programs from unrecovered original versions. Original classifier/index checksums and raw tool versions are **not recovered**. A modern replacement database is not an exact reproduction.
- Official FAPROTAX mapper: castor and optparse were absent; original optional package pins are not established. No originally used castor version is inferred from a “tested with” comment.
- No package installation, read download or HPC job is launched automatically.

## Earliest verified reproducible stage

With `M1` and `CLIN`: frozen specimen profiles, model inputs and accepted clinical evidence. With `V3`: original source TSS/CLR integration and full-composition alpha. With `V2`: accepted official annotation replay and fresh source-scale trait contributions. The earliest raw-read stage remains **BLOCKED** until original reads, compatible tools and matching databases are supplied and separately validated.
