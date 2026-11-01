# Standalone historical branches

These branches preserve source provenance and study-specific functionality. They do not participate in the main `run_pipeline.py` checkpoint contract. Use the module-root entry for routine processing.

| Location | Purpose | Current limits |
| --- | --- | --- |
| `amplicon/selected_se_batch.sh` | Historical batched SE dereplication route | Isolated research entry; not executed during repository preparation |
| `amplicon/PRJCA022064.sh` | Special naming and best-mate branch | Confirm the study's input design; not executed during repository preparation |
| `solid/` | SOLiD preparation, BFAST indexing, and chunked host removal | Partial route; no complete downstream classification evidence |
| `pooled_diversity/integrated_diversity.R` | Diversity of an existing merged Genus matrix | Does not construct/merge matrices; full cohort not rerun |

Shell entries require explicit `AREK_LEGACY_SOURCE_MODE=isolated_research` and real `RAW_DATA_ROOT`, `RESULT_ROOT`, `DATABASE_ROOT`, and other route-specific paths. SOLiD also needs a separate BFAST environment and references. See the [user guide](../docs/user_guide.md) and [source map](../docs/source_script_map.tsv).

Run the merged-matrix entry from the module root:

```bash
AREK_POOLED_SCOPE=verified_same_host_site_pipeline \
Rscript legacy/pooled_diversity/integrated_diversity.R /path/to/new_diversity_output /path/to/existing_integration_root
```

Confirm the same host/body site, comparable processing, supported Genus identity, and inclusion scope. This branch retains historical unfiltered Alpha and filtered Beta inputs, which differ from the main workflow's defaults.
