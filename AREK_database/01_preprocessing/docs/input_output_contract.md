# Input/output contract and downstream interfaces

## Inputs

| Input | Requirement |
| --- | --- |
| Configuration JSON | Merged with `config/historical.json`; unknown keys rejected; relative paths resolved against the configuration directory |
| FASTQ sample table | TSV with `dataset_id`, `sample_id`, `run_id`, `layout`, `r1`, `r2`; read paths relative to the table; empty `r2` for SE and both mates for PE |
| Biological identity | Populate `specimen_id`/`subject_id` from evidence; resolve multiple lanes/runs upstream rather than merging automatically |
| Metadata | Unique sample IDs in the first column; configured additional fields; unknown groups remain missing |
| Counts matrix | Feature-by-sample matrix with unique IDs and finite nonnegative values; structural missingness is rejected; Chao1 is NA for noninteger estimates |
| Relative matrix | Complete 0–1 denominator with columns summing to 1; percentages or truncated tables require explicit upstream treatment |
| Taxonomy | Unique feature IDs and a `Taxon` column, or explicit `native_level` for a native-rank matrix |
| External references | SILVA classifier for amplicon; Kraken/Bracken database and configured Bowtie2 host index for shotgun; see [databases.md](databases.md) |

A `dataset_id` identifies a study/host/body-site/processing unit. Different body sites in one study should be registered separately. Discovery proposes file mappings; metadata establishes biological identity.

## Standard outputs

Outputs are written under `--out/DATASET_ID/`.

| Relative path | Contents and units | Downstream use |
| --- | --- | --- |
| `upstream/` or matrix-route `input/` | QZA/source features/taxonomy or imported matrix copies | Source tracing; cross-study OTU IDs do not establish equivalence |
| `abundance/LEVEL/counts.tsv` | Rank counts; Bracken values are estimated reads | Count-based analysis; absent for relative-only input |
| `abundance/LEVEL/relative.tsv` | 0–1 abundance using all original features as denominator | Full abundance profile |
| `abundance/LEVEL/filtered_counts.tsv` | Counts after configured filtering | Default diversity input; absent for relative-only input |
| `abundance/LEVEL/filtered_original_denominator.tsv` | Retained taxa using the original total | Audit filtering loss; column sums may be below 1 |
| `abundance/LEVEL/filtered_tss.tsv` | Retained taxa renormalized to 0–1 | Analyses explicitly requiring the retained-taxon denominator |
| `abundance/LEVEL/taxonomy.tsv`, `feature_to_taxon.tsv` | Full-lineage keys and source-feature mapping | Alignment by supported taxonomic identity |
| `abundance/LEVEL/contract.json`, `mass_balance.tsv`, `filter_audit.tsv` | Units, totals, and filtering/classification audits | Verification before downstream import |
| `diversity/LEVEL/` | Alpha/Gamma, three distances, PCoA coordinates/eigenvalues, tests, and sample audits | Results and downstream interpretation |
| `diversity/LEVEL/figures/{pdf,png}/`, `captions.tsv` | Figures and source-consistent captions | Presentation |
| `effective_config.json`, `state/*.json`, `execution.jsonl`, `logs/` | Parameters, hashes, source/tool fingerprints, commands, and exit states | Reproducibility and execution evidence |

Ranks without named assignments receive `not_applicable.json`. This does not establish computed abundance/diversity for that rank. Zero libraries are explicitly excluded from diversity with a recorded reason.

Later modules should read tables together with contracts, taxonomy, and metadata, and distinguish original denominators from filtered TSS. Relative abundance cannot reconstruct counts. Missing reports are not zero-abundance reports. Drinking/disease labels require explicit metadata.
