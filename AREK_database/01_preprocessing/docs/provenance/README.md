# Source provenance and repository preparation

The original input is `AREK_FASTQ_to_Abundance_Diversity_pipeline_20261030.zip`. Repository preparation took place on 2026-10-31 (America/Chicago).

| File | Purpose |
| --- | --- |
| [reorganization.json](reorganization.json) | Original archive SHA256, source-file count, preparation scope, and adaptation notes |
| [file_relocation.tsv](file_relocation.tsv) | Retained source paths, current locations, original hashes, and identical/adapted status |
| Module-root `FILE_MANIFEST.tsv` | Current module paths, sizes, and hashes, excluding the manifest itself |

The initial preparation retained computational source, configuration, environments, and example inputs byte-for-byte where marked `identical`. Six historical entries were moved under `legacy/`; three validation programs received evidence-output-path adaptations. Documentation and source maps were updated and historical server paths removed from inventories.

The English-language update translates all repository documentation, renames current language-suffixed documents, and updates links, source notes, comments, and error messages. Legacy missing-group aliases in JSON use Unicode escapes so decoded values remain unchanged. Original source hashes are preserved, while files changed by translation are marked `adapted`. Numerical methods, thresholds, environment locks, and biological fixtures are unchanged.

Historical logs, validation snapshots, cohort/result inventories, project scans, QZA intermediates, and test tables/figures remain in the original archive. The repository contains source, configuration, dependencies, documentation, public verification fixtures, and explicitly scoped check summaries.

Some `historical_*` modules remain dependencies of the main workflow, including Bracken lineage parsing and metric formulas. They are kept in `scripts/`, separately from independently executable `legacy/` branches. Historical source filenames in provenance tables identify the original archive and are retained as evidence.
