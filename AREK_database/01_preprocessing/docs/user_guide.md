# AREK preprocessing user guide

This module implements the first stage of the AREK database: FASTQ processing, taxonomic abundance, and diversity. It was prepared from the 2026-10-30 source archive and organized for the repository on 2026-10-31 (America/Chicago). Start with the [quick start](../README.md) and [step reference](steps.md). Inputs, references, and outputs are supplied through configuration or command-line arguments. Historical execution evidence is distinguished from checks performed on the current repository.

## Scope and source

The main entry point is `python run_pipeline.py`. It connects recovered computational modules using local, sequential dataset scheduling; individual tools use the configured thread count. It does not submit cluster jobs or stop unrelated processes. By default, it refuses writes into legacy production directories named `04_export` or `03_vsearch_otu97`.

| Route or entry | Method | Support and evidence |
| --- | --- | --- |
| `amplicon_vsearch97`, `layout=PE` | QIIME2 import, vsearch merge-pairs, dereplication, 97% de novo OTUs, UCHIME, table/sequence filtering, and sklearn SILVA classification | Complete FASTQ-to-diversity route; a real PE subset was validated in the source package |
| `amplicon_vsearch97`, `layout=SE` | Explicit biological read, import, and the same OTU97/classification steps | Complete SE route; historical software verification explicitly selected R1 from real PE data |
| `shotgun_kraken_bracken`, `shotgun.variant=v2` | fastp, Bowtie2 host removal, Kraken2, and independent Bracken estimates by rank | Recovered PRJNA613834 PE route; source, dependency, parameter, and database checks, without a full rerun |
| `shotgun.variant=original` | Recovered PRJCA021196 fastp settings, including PE correction, and PE/SE host-removal routes | Source restored; full chain not rerun |
| `matrix` | Existing features plus taxonomy, or an explicit `native_level`, followed by abundance/diversity | Counts and complete 0–1 relative abundance supported; real matrices validated |
| `legacy/pooled_diversity/integrated_diversity.R` | Standalone diversity analysis of an existing merged Genus matrix | Reads a verified integration root; full merged cohort not rerun |
| `legacy/amplicon/selected_se_batch.sh`, `legacy/amplicon/PRJCA022064.sh` | Historical batch SE, naming, and best-mate branches | Isolated research use; no main-workflow checkpoint guarantee; not executed during repository preparation |
| `legacy/solid/` | ABI SOLiD 5500 color-space diagnostics and BFAST 0.5.6 hg38 host removal | Partial branch; successful downstream base-space classification has not been established |

[source_script_map.tsv](source_script_map.tsv) records original hashes, current entries, dependencies, and adaptations. Historical logs, QZA provenance, full study inventories, and result snapshots remain in the original archive.

Version labels such as `v3` or `latest` do not establish a successful execution. PRJDB10075 has a successfully terminated v5 log, but the recovered table/taxonomy QZA provenance dates to 2026-06-26/27, after the May log. PRJNA909793 has a v5 SIGSEGV log and a historical batch containing both WGS and amplicon inputs. Confirm the actual sample strategy in the manifest. The recovered OTU97 method is retained; no new denoising or primer-removal method is introduced.

Primers, amplified regions, platforms, or database builds without recovered evidence remain unspecified. The main OTU97 route has no separate cutadapt or quality-truncation step. Its PE quality settings are the recovered merge defaults; the SE route adds no unsupported quality-filtering step. If primers still need removal, establish the study's preprocessing separately.

Differential abundance, regression modules, FAPROTAX, PICRUSt2, MICOM, patient subtyping, and prediction models are outside this module. Covariate adjustment within diversity is limited to the configured group-distance test.

## Installation and resources

The environment locks target Linux x86_64. Each `envs/*.linux-64.explicit.txt` file contains the required seed packages and their transitive dependencies recovered from installed conda metadata, including URLs, builds, and checksums. Corresponding `*.versions.json` files record versions/builds and package SHA256/MD5. BFAST and QIIME2 use separate environments.

Run these commands from the repository root:

```bash
cd AREK_database/01_preprocessing
conda create -n arek_amplicon --file envs/amplicon.linux-64.explicit.txt
conda activate arek_amplicon

# For existing matrices and diversity only:
# conda create -n arek_analysis --file envs/analysis.linux-64.explicit.txt
# For shotgun processing:
# conda create -n arek_shotgun --file envs/shotgun.linux-64.explicit.txt
python run_pipeline.py --help
```

Recorded installed versions include Python 3.10.20, QIIME2 2024.10.1, q2-vsearch/q2-feature-table/q2-feature-classifier 2024.10.0, scikit-learn 1.4.2, NumPy 1.26.4, pandas 2.2.2, SciPy 1.13.0, matplotlib 3.8.4, R vegan 2.7-1, and permute 0.9-8. Consult the lock/version files for exact route-specific versions. A 2026.4 URL in an old environment-creation script was not treated as evidence of the installed environment.

These locks were checked against recovered metadata; clean network installations were not rerun during repository preparation. Missing remote packages or network failures should cause an explicit conda error. Other platforms require dependency resolution and validation.

Reference databases are external. See [databases.md](databases.md) and [database_inventory.tsv](database_inventory.tsv) for resource identities, hashes, and missing provenance.

```bash
export AREK_CLASSIFIER=/path/to/silva-138-99-nb-classifier.qza
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
```

Example resource settings are new portable defaults: two threads, an 8 GiB memory preflight request, and at least 1 GiB free output space. They do not enforce process RSS and are not measured peaks. Historical subset timings are in `validation/final_tests/real_chain_measurements.json` in the original archive; peak RSS was not measured. Full amplicon clustering may require tens to hundreds of GiB depending on unique sequences. Shotgun analysis also loads Kraken and host references. Pilot the workload and adjust resources using actual measurements.

## Inputs, discovery, and sample identity

The scanner recognizes `.fastq`/`.fq` and `.gz`/`.bz2`/`.xz`. Actual compression support depends on the downstream tool; uncompressed or gzip FASTQ is recommended. Mate suffixes include `_R1_001`/`_R2_001`, `_R1`/`_R2`, `_1`/`_2`, `.1`/`.2`, `_read1`/`_read2`, and `_f1`/`_r2`. Files without a mate suffix are SE candidates. Multiple lanes are not merged automatically. Orphan mates, duplicate IDs, and simultaneous SE/PE assignments produce conflicts. Discovery cannot establish whether runs belong to the same specimen or subject; confirm that mapping from metadata.

From the module root:

```bash
python run_pipeline.py scan --root /path/to/data_root --out /path/to/scan_report
python run_pipeline.py scan --root examples/real_fastq --out /tmp/arek_scan
```

Outputs are `files.tsv`, `samples.tsv`, `conflicts.tsv`, `artifacts.tsv`, and `step_status.tsv`. Historical artifacts are marked as present but unverified. Recognized checkpoints have their output hashes checked, but still require validation against the current configuration. Run `validate` before accepting completion. Conflicts return exit code 2 while preserving discovery results for independent datasets.

The sample manifest is TSV. Relative read paths are resolved against the manifest location. A `dataset_id` identifies a verified study/host/body-site/processing unit; saliva and feces from one study should use separate units.

```tsv
dataset_id	sample_id	run_id	specimen_id	subject_id	layout	r1	r2
PRJDB10075	DRR234531	DRR234531			PE	real_fastq/PRJDB10075/DRR234531_1.fastq.gz	real_fastq/PRJDB10075/DRR234531_2.fastq.gz
```

See [real_samples.tsv](../examples/real_samples.tsv) for a tab-delimited file. Required columns are `dataset_id`, `sample_id`, `run_id`, `layout`, `r1`, and `r2`. SE requires an empty `r2`; PE requires both mates. Populate `specimen_id` and `subject_id` only when supported by evidence. Empty values do not establish independence. Multiple runs or lanes from one specimen require an explicit upstream technical-replicate policy. Run IDs do not substitute for subject IDs in independent statistical analysis.

Safe identifiers can contain letters, numbers, dots, hyphens, and underscores. PRJNA, PRJDB, PRJCA, SRP, CRA, and custom IDs are supported.

Metadata is TSV with unique sample IDs in its first column. Additional columns are selected through configuration, for example `group`, `subject_id`, `study`, `host_species`, `body_site`, `age`, and `sex`. Unmatched and extra samples are audited; group labels need not exist for all samples. The metadata template is a formatting example. Conflicting run/specimen/subject mappings, missing configured columns, duplicate IDs, or no shared samples produce an error.

Statistical fields are under `diversity`: `group_column`, `subject_column`, `block_column`, `covariates`, and `design`. Designs are `descriptive`, `independent`, or `within_subject`. `missing_group_values` explicitly lists unknown labels; their original values are retained in `source_group_value` and normalized to missing. Unknown labels are never assigned to Control. Legacy missing-value aliases remain encoded as Unicode escapes in the English configuration to preserve existing metadata behavior.

Configuration-relative paths use the configuration directory; `--out` and command-line reference paths use the working directory; internal modules use the package root. Defaults merge from `config/historical.json`; unknown keys are rejected. Historical method settings are distinguished from portable resource defaults in [parameters.md](parameters.md). Fixed metric definitions are recorded in `config/metric_contract.tsv`.

## Commands and execution

Activate the appropriate environment, set `AREK_CLASSIFIER`, and run:

```bash
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
python run_pipeline.py dry-run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
python run_pipeline.py run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
python run_pipeline.py validate --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
```

`dry-run` prints external commands, internal stages, and parameters without launching tools, creating output directories, or submitting jobs. Rank applicability that depends on future taxonomy is determined during execution.

A manifest may contain multiple datasets. Use `--dataset all`, a comma-separated list, or a single ID. Datasets run separately. A dataset failure does not cancel independent datasets; the batch returns 2 if any fail and 0 if all succeed. Use separate configurations for incompatible layouts, routes, or study designs.

```bash
# Reuse only stages whose content and current fingerprints are valid.
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe

# Execute a range after validating its upstream checkpoints.
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe \
  --from-stage levels --to-stage diversity

# Explicitly select R1 as SE input and use a separate output directory.
python run_pipeline.py run --config config/real_se_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_se
```

Stage order is `preprocess`, `features`, `taxonomy`, `levels`, `diversity`, `plot`. Skipped upstream stages still require valid checkpoints. Use the matrix route when importing existing tables. `run` re-executes selected stages; `resume` reuses valid content.

Each dataset uses a kernel `flock` to prevent duplicate execution. Running-state JSON records parent and child PIDs; a success checkpoint is written only after success. Tool children inherit the lock, so a scheduler crash does not permit a duplicate launch while a child still holds it. Interactive interruption terminates the workflow's separate tool process group. `validate` refuses a locked output. Stale PID files alone do not establish a live job. Confirm the full child process tree has exited before resuming an interrupted task.

Existing matrices require no classifier:

```bash
python run_pipeline.py check --config config/matrix_example.json --dataset REAL_MATRIX
python run_pipeline.py run --config config/matrix_example.json --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py validate --config config/matrix_example.json --dataset REAL_MATRIX --out /tmp/arek_matrix
```

Point `matrix` and `taxonomy` to your tables. With taxonomy, aggregate only provided ranks. For a Genus-only matrix, set `native_level: "Genus"` and `levels: ["Genus"]` while retaining source taxonomic keys. For relative-only input, set `matrix_unit: "relative"` and `diversity.input` to `"filtered_tss"` or `"relative"`. Input must represent a complete 0–1 denominator with columns summing to 1. Percentage or truncated tables require explicit upstream handling.

Shotgun examples are `config/shotgun_v2_example.json` and `config/shotgun_original_example.json`. Supply the actual manifest and references:

```bash
python run_pipeline.py check --config config/shotgun_v2_example.json \
  --dataset PRJCA021196 --kraken-db /path/to/kraken_db --host-index /path/to/hg38/hg38
python run_pipeline.py dry-run --config config/shotgun_v2_example.json \
  --dataset PRJCA021196 --kraken-db /path/to/kraken_db --host-index /path/to/hg38/hg38 \
  --out /path/to/new_shotgun_output
# After a successful check, replace dry-run with run, then validate.
# The user_fastq paths in the example manifest must be supplied by the user.
```

Figures can be redrawn from existing TSV outputs:

```bash
python run_pipeline.py redraw --folder /tmp/arek_matrix/REAL_MATRIX/diversity/Genus \
  --dataset REAL_MATRIX --level Genus
# To color by existing metadata, add:
# --metadata /path/to/metadata.tsv --group-column group
```

`redraw` updates figures/captions without changing metric tables and can invalidate the original plot checkpoint. Use the original configuration with `resume --from-stage plot` to restore workflow-managed figures and their fingerprint.

## Outputs and metrics

Outputs are under `--out/DATASET_ID/`:

- `upstream/`: amplicon QZA artifacts, feature BIOM/TSV, representative FASTA, and taxonomy; shotgun clean/nonhost reads, reports, and independent P/F/G/S Bracken tables. The matrix route preserves imported inputs under `input/`.
- `abundance/LEVEL/`: counts where available, relative abundance, filtered tables, full-lineage keys, feature-to-taxon mappings, filter/mass-balance audits, and contracts. Unsupported ranks contain `not_applicable.json`.
- `diversity/LEVEL/`: sample Alpha, three distances, all PCoA axes/eigenvalues, Gamma summaries/tests, Alpha tests, PERMANOVA, dispersion, metadata/sample audits, metric definitions, and R session records.
- `diversity/LEVEL/figures/{pdf,png}/` and `captions.tsv`: Alpha plots, three PCoAs, and three Gamma plots; PNGs use 300 dpi. Descriptive unknown-group labels and sample counts follow the underlying tables.
- `state/*.json`, `execution.jsonl`, `logs/`, and `effective_config.json`: hashes, resolved parameters, code/tool fingerprints, argv, resources, timing, PIDs, exit codes, and validation results.

Taxonomic keys use the full observed lineage through the requested rank. Identical terminal names from different lineages remain distinct. Empty ranks retain empty labels, classification status, and available ancestor lineage. Original features remain traceable; unrelated studies' OTU/ASV IDs are not treated as shared identities. Missing taxonomy is recorded in contracts/mappings. Structural missingness, nonnumeric, negative, and nonfinite matrix values are rejected. Taxids explicitly absent from complete shotgun reports can become true zeros after joining; missing reports cannot.

Aggregation must conserve totals. Default filtering requires prevalence >=0.10 and mean relative abundance >=1e-4. `relative.tsv` uses all original features, including unclassified entries, as its denominator. `filtered_original_denominator.tsv` retains that denominator, so columns can sum below 1. `filtered_tss.tsv` renormalizes retained taxa to a column sum of 1. Both are retained and checked. Zero libraries are excluded from computable diversity and recorded in contracts/definitions.

Default diversity input is `filtered_counts`; `counts`, `relative`, or `filtered_tss` can be selected explicitly. Each rank is analyzed separately.

| Metric | Definition and input |
| --- | --- |
| Observed/richness | Taxa strictly above zero in the selected filtering state; rank richness |
| Shannon | `-sum(p * ln(p))`, with natural logs and probabilities from the selected sample column |
| Simpson | `1 - sum(p^2)` |
| Chao1 | `S + f1*(f1-1)/(2*(f2+1))`; integer counts only; NA for relative or noninteger estimated counts |
| Bray–Curtis | `sum(abs(x-y))/sum(x+y)` on the selected matrix |
| Jaccard | `1 - intersection/union` on positive-value presence/absence |
| Aitchison | Euclidean distance after sample-centered `ln(x+0.5)`; the pseudocount uses the selected matrix's units |
| PCoA | Full eigendecomposition, including exported negative eigenvalues; figure percentages use the sum of positive eigenvalues; no silent distance correction |
| Gamma | Sum columns within each group, then compute Observed, Shannon, and 1-D; depth-weighted for counts and equal-sample-weighted for TSS. `All_samples` is descriptive |

The main Python workflow performs no rarefaction and implements no UniFrac, phylogenetic-tree metric, or NMDS. Seed 1 is the historical permutation default. Example 99-permutation runs are new verification settings; the historical default is 999. A separate old microeco script without confirmed executed outputs was not treated as a validated route.

Two-group Gamma tests preserve group sizes in independent label permutations, use a two-sided absolute difference with a plus-one p value, and apply BH correction. Repeated, blocked, or covariate-adjusted designs retain descriptive Gamma summaries with an explicit NA reason for unsupported tests.

## Statistical designs and historical corrections

The old `analysis_core.betadisper_test` used two PCoA axes and Kruskal testing. The audit copy is retained, but the main workflow uses vegan `betadisper(type="median")` and `permutest`, separately from `adonis2(by="margin")`. All coordinates, eigenvalues, and distance TSVs are exported. These corrected tests do not claim exact numerical reproduction of the old statistical results.

Independent Alpha comparisons use two-group Mann–Whitney (R asymptotic implementation) or multi-group Kruskal, with BH correction. P values can differ from old SciPy automatic exact-test choices. Two-group `within_subject` designs use paired Wilcoxon only when each subject contributes one observation per group. Unsupported repeated or adjusted Alpha designs retain descriptions and an NA reason. PERMANOVA supports configured covariates, complete-case auditing, and design-matrix rank checks. Missing group/block values, fewer than two observations per group, rank deficiency, or absent residual degrees of freedom are reported explicitly.

Known repeated subjects require a valid design. Repeats detected in manifests or common metadata fields trigger preflight errors when unconfigured. `within_subject` requires subject/block alignment and exchangeable group labels within subjects. A group invariant within each subject, such as repeated visits in fixed case/control groups, produces `group_invariant_within_subject`; it cannot use within-subject label permutations. Select appropriate specimens or use a study-specific model. Unprovided subject mappings cannot be inferred from code.

Single-arm and unknown-label datasets still receive Alpha, distance, PCoA, and Gamma descriptions. Each test exports its NA reason and actual complete-case sample size.

## Existing merged-matrix branch

The standalone `integrated_diversity.R` accepts an existing integration root containing the expected Genus matrices/metadata under `FINAL` or `02_merged_cohorts`. It retains unfiltered Alpha, filtered Beta, Dataset adjustment, seed 20261005, and 999 permutations. Its input settings differ from the main counts workflow.

```bash
AREK_POOLED_SCOPE=verified_same_host_site_pipeline \
Rscript legacy/pooled_diversity/integrated_diversity.R /path/to/new_diversity_output /path/to/existing_integration_root
```

Verify host, body site, comparable processing, Genus identities, and inclusion scope before use. Missing-value and repeated-subject guards reject inappropriate inputs. Rounding count-like values into Chao1 counts is disabled. Chao1 requires explicit `AREK_INTEGRATION_INTEGER_COUNTS=1` and an actually integer-count input. The full merged cohort and its source-matrix construction were not rerun during repository preparation.

Legacy shell branches require real `RAW_DATA_ROOT`, `RESULT_ROOT`, and `DATABASE_ROOT` paths plus `AREK_LEGACY_SOURCE_MODE=isolated_research`. They use isolated output directories and do not provide main-workflow resume guarantees. SOLiD additionally requires its separate environment and BFAST references; see [databases.md](databases.md).

## Validation and troubleshooting

Historical source-package evidence is in `validation/VALIDATION_zh.md` and `validation/final_tests/` in the original archive. The public PRJDB10075 fixtures contain the first 2,000 original read pairs per run, recovered from imported QZA. Accessions and hashes are in [examples/provenance.json](../examples/provenance.json). Metadata templates and statistical fault fixtures use software-test labels. See [validation/README.md](../validation/README.md) for current checks.

```bash
# Full verification requires the amplicon environment and classifier.
python validation/validate_release.py --work-root /tmp/arek_recheck_full \
  --classifier "$AREK_CLASSIFIER"

# Matrix/fault/design verification still requires the analysis and R dependencies.
python validation/validate_release.py --work-root /tmp/arek_recheck_matrix --matrix-only
```

Validation covers help/configuration, real pairing discovery, lane/orphan conflicts, no-write dry-runs, dependency/reference failures, first execution, resume, selective Genus recovery, independent dataset continuation, explicit SE, matrix/relative inputs, repeated-subject guards, and independent numerical checks. Archive validation additionally checks extraction CRCs, payload hashes, internal paths, and representative entries.

| Error or state | Action |
| --- | --- |
| `classifier QZA missing` | Supply a verified file with `--classifier` |
| `missing executable/R dependencies` | Activate the appropriate locked environment |
| `orphan/duplicate/multiple lanes` | Confirm run/lane/specimen mappings and supply an explicit manifest |
| `mate ID mismatch/invalid FASTQ record` | Inspect the named file/record, recover valid input, and resume |
| `color-space reads detected` | Use the separate SOLiD diagnostic branch |
| `relative input does not sum to 1` | Establish the actual units and denominator before import |
| `invalid or stale` | Inspect state/logs and resume the affected branch; stage selection retains upstream checks |
| `already running` | Inspect the process IDs and existing dataset lock |
| `not_applicable` | The source taxonomy does not provide named taxa at that rank |
| Insufficient resources | Inspect available memory/disk and revise the pilot/resource settings |

`FILE_MANIFEST.tsv` lists current module paths, bytes, and SHA256, excluding itself. Source provenance retains the original archive identity and original hashes; adapted files are marked accordingly. Full study FASTQ, historical runtime environments, and large references are external. Untested branches remain explicitly identified.
