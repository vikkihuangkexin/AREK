# Processing steps and source modules

Run commands from `AREK_database/01_preprocessing/`. Production execution uses `run_pipeline.py`; computational Python modules with internal dependencies are called through the workflow.

## Step reference

| Step | Input, operation, and output | Entry point or source | Stage selector |
| --- | --- | --- | --- |
| 00 | Raw directory, explicit naming/pair discovery, and file/sample/conflict TSVs | `scan`; `workflow/scanner.py` | Standalone command |
| 01 | Configuration, sample table, references, preflight checks, and execution plan | `check`, `dry-run`; `workflow/cli.py` | Standalone commands |
| 02 | FASTQ read processing and classification preparation, yielding QZA or sample reports | `workflow/routes.py`; shotgun calls `scripts/preprocess/host_filter.sh` | `preprocess` |
| 03 | Filtered amplicon table/sequences to feature TSV/BIOM/FASTA; Bracken reports to rank matrices; or existing matrix import | `workflow/routes.py`, `workflow/cli.py`, and lineage helpers in `scripts/abundance/historical_bracken_merge.py` | `features` |
| 04 | Amplicon representative sequences and classifier to taxonomy QZA/TSV | `workflow/routes.py` | `taxonomy` |
| 05 | Feature/taxonomy or native-level matrix to rank aggregation, filtering, and denominator audits | `scripts/abundance/levels.py` | `levels` |
| 06 | Selected rank matrix and metadata to Alpha, distances, PCoA, Gamma, and group tests | `workflow/cli.py`; `scripts/diversity/compute.py`, `historical_metrics.py`, `historical_gamma.py`, `statistics.R` | `diversity`, including metadata and statistics substages |
| 07 | Diversity TSVs to PDF/PNG and captions | `scripts/diversity/plot.py` | `plot`; standalone `redraw` also available |
| 08 | Existing outputs, original inputs, configuration, and tools to validation status | `validate`; `workflow/common.py`, `workflow/cli.py` | Standalone command with optional stage range |

`workflow/common.py` handles command execution, logging, content fingerprints, and checkpoints. `run` re-executes selected stages, `resume` reuses valid fingerprints, and `validate` checks existing outputs.

## Main routes

**Amplicon:** explicit PE/SE mapping, QIIME2 import, PE merge where applicable, dereplication, 97% OTU clustering, UCHIME, table/sequence filtering, feature export, SILVA sklearn classification, abundance, diversity/statistics, and figures. The `preprocess` selector includes import through chimera filtering; the six selectors are not six individual external commands.

**Shotgun:** fastp, configured host removal, Kraken2, independent Bracken rank estimates, same-rank merging, abundance, diversity/statistics, and figures. Sample-level external processing belongs to `preprocess`; rank merging belongs to `features`. P/F/G/S are estimated separately. `variant=v2` requires PE; `original` supports PE/SE. This route does not use SILVA.

**Matrix:** strict import, preservation of source features/taxonomy, rank aggregation or explicit `native_level`, diversity/statistics, and figures. A Genus-only table produces no unsupported Phylum or Species matrix. Matrix import establishes no FASTQ preprocessing reproduction.

## Downstream interfaces

Later modules should read verified outputs under `--out/DATASET_ID/` and retain dataset ID, host/body site, rank, units, filtering, sample identities, and method provenance. See the [input/output contract](input_output_contract.md).

Each dataset is processed separately. Cross-study feature IDs are not interchangeable. A future Genus integration module should explicitly reconcile taxonomic identity, host, body site, processing, and original denominators.

The machine-readable index is [step_registry.tsv](step_registry.tsv). Separate historical entries are documented in [legacy/README.md](../legacy/README.md).
