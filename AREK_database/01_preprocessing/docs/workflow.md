# Workflow dependencies and checkpoint contract

## Source responsibilities

| Source | Responsibility |
| --- | --- |
| `run_pipeline.py` | Unified command-line entry |
| `workflow/cli.py` | Arguments, configuration, batches, metadata, and stage ranges |
| `workflow/scanner.py` | Naming discovery, FASTQ structure, and pairing checks |
| `workflow/common.py` | Dataset locks, execution records, content validation, and checkpoints |
| `workflow/routes.py` | Amplicon QIIME2/vsearch and shotgun fastp/host/Kraken2/Bracken commands |
| `scripts/abundance/historical_bracken_merge.py` | Recovered lineage parsing with exact-rank handling |
| `scripts/abundance/levels.py` | Taxonomic aggregation, filtering, and mass balance |
| `scripts/diversity/compute.py` | Diversity computation and output export |
| `scripts/diversity/historical_metrics.py` | Recovered metric formulas |
| `scripts/diversity/historical_gamma.py` | Gamma definitions and independent permutations |
| `scripts/diversity/statistics.R` | Corrected statistical tests |
| `scripts/diversity/plot.py` | Figures and standalone TSV-based redraw |

Original hashes and packaged entries are listed in [source_script_map.tsv](source_script_map.tsv). Historical logs and source execution paths remain in the original archive. Required local Python/R/shell modules are packaged; external dependencies are command-line tools, input data, and references supplied by the user. Tool environments are recorded in explicit locks.

## Route dependencies

The amplicon route runs preparation, import, optional PE merge, dereplication, clustering, chimera detection, and table/sequence filtering. Filtered tables produce features; filtered sequences produce classification/taxonomy. Features and taxonomy feed each rank's abundance, metadata, diversity, statistics, and plot stages.

The shotgun route runs fastp, configured host removal, Kraken2, and independent Bracken P/F/G/S estimates per sample. Same-rank reports are merged into native-rank matrices, followed by the shared diversity workflow. Species estimates are not summed to represent genus/family Bracken estimates.

The matrix route strictly reads source features, accepts taxonomy or explicit `native_level`, and proceeds to abundance/diversity. Genus-only inputs provide no unsupported Phylum/Species output.

## Successful-stage contract

A successful stage requires matching inputs, method parameters, source, and tool versions; nonempty outputs; valid structure/numerics/IDs; and matching hashes for every declared output. QZA validation checks ZIP CRC/data structure and feature/representative-sequence/taxonomy identity consistency. Directory existence or a zero exit status alone is insufficient.

Long tasks record parent/child PIDs and wait synchronously for completion. A dataset-level `flock` prevents duplicate launches. Tool children inherit the lock.

`resume` uses stage-relevant source fingerprints. Plot-only changes do not invalidate FASTQ stages. Shared orchestration or input-reader changes can conservatively invalidate multiple stages. A damaged output triggers its stage and explicit downstream dependencies while allowing unrelated ranks to remain reusable. Historical source-package testing deleted a Genus distance and recovered only the Genus diversity/statistics/plot branch.

Upstream stages before `--from-stage` are checked even when not executed. Existing external tables should enter through the matrix route. `validate` checks the selected dependency graph, including upstream content, without recomputing products or submitting jobs. `check` preserves actionable preflight errors. Batch execution catches errors per dataset and returns 2 if any fail.

`effective_configuration_sha256` records the complete resolved configuration. Each stage fingerprint includes its relevant parameters, input content, source, tools, and resources. Versions come from actual Python packages, conda metadata, executable paths, and R sessions.

Standalone merged-Genus R and legacy shell branches do not participate in this checkpoint contract. Their applicability and unverified scope are documented separately.
