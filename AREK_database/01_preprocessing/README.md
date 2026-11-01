# 01 · Preprocessing: FASTQ to taxonomic abundance and diversity

The first processing module of the AREK database converts explicitly registered datasets into traceable taxonomic abundance tables, diversity measurements, statistical results, and figures. Its numerical methods are retained from the recovered source package. Documentation, comments, and user-facing messages are maintained in English.

## Supported routes

| Input | `route` | Method |
| --- | --- | --- |
| Paired-end or single-end amplicon FASTQ | `amplicon_vsearch97` | QIIME2, vsearch 97% OTU clustering, and SILVA sklearn classification |
| Shotgun FASTQ | `shotgun_kraken_bracken` | fastp, configurable host removal, Kraken2, and independent Bracken estimates for each rank |
| Existing abundance matrix | `matrix` | Feature matrix with taxonomy, or a matrix with an explicit `native_level` |

The amplicon route uses OTU97 clustering. SOLiD color-space processing is a separate, partial [legacy branch](legacy/README.md). Historical amplicon subset validation and current repository checks are documented separately in [validation](validation/README.md). The full shotgun chain and full merged cohorts were not rerun during repository preparation.

## Quick start

Start in the repository root. The dependency locks target Linux x86_64. Reference databases must be supplied separately.

```bash
cd AREK_database/01_preprocessing

# Restore the environment for the existing-matrix route.
conda create -n arek_analysis --file envs/analysis.linux-64.explicit.txt
conda activate arek_analysis

python run_pipeline.py --help
python run_pipeline.py check --config config/matrix_example.json --dataset REAL_MATRIX
python run_pipeline.py dry-run --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py run --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py validate --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
```

For amplicon reads, restore and activate the environment defined by `envs/amplicon.linux-64.explicit.txt`, then supply a verified classifier:

```bash
export AREK_CLASSIFIER=/path/to/silva-138-99-nb-classifier.qza
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
python run_pipeline.py run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_pe
```

Use `resume` to reuse valid checkpoints. Stage ranges require valid upstream checkpoints:

```bash
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_pe \
  --from-stage levels --to-stage plot
```

Relative FASTQ paths are resolved against the sample table; relative configuration paths are resolved against the configuration file. See the [input/output contract](docs/input_output_contract.md).

## Processing steps

| Step | Operation | CLI command or stage | Main outputs |
| --- | --- | --- | --- |
| 00 | Discover and confirm sample mappings | `scan` | Candidate sample tables and pairing conflicts |
| 01 | Check inputs and inspect the execution plan | `check`, `dry-run` | Input, dependency, and database checks; planned commands |
| 02 | Process reads | `preprocess` | Amplicon import, merging, dereplication, clustering, and chimera filtering; or shotgun cleaning, host removal, and classification reports |
| 03 | Produce feature tables | `features` | OTU tables, independent Bracken rank tables, or imported matrices |
| 04 | Assign taxonomy | `taxonomy` | Amplicon taxonomy; existing shotgun taxonomic identities |
| 05 | Build taxonomic abundance tables | `levels` | Counts, relative abundance, filtered tables, mappings, and mass-balance audits |
| 06 | Compute diversity and statistics | `diversity` | Alpha, distances, PCoA, Gamma, PERMANOVA, PERMDISP, and sample audits |
| 07 | Draw figures | `plot` | PDF, PNG, and captions |
| 08 | Verify outputs | `validate` | Content, parameter, source, and tool fingerprint checks |

The six stage selectors are `preprocess`, `features`, `taxonomy`, `levels`, `diversity`, and `plot`. Each route expands them into its own tool calls. See the [step reference](docs/steps.md).

## Directory layout

| Path | Contents |
| --- | --- |
| `run_pipeline.py` | Unified command-line entry point |
| `workflow/` | Configuration, dataset scheduling, routes, checkpoints, and input checks |
| `scripts/preprocess/` | Host-removal helper for the main workflow |
| `scripts/abundance/` | Bracken merging, rank aggregation, and abundance filtering |
| `scripts/diversity/` | Alpha/Beta/Gamma metrics, statistical tests, and plots |
| `config/` | Defaults, route examples, and metric contract |
| `envs/` | Route-specific Linux conda locks and version records |
| `examples/` | Sample tables, metadata template, public read subsets, and matrices |
| `validation/` | Executable checks that write evidence to a supplied work directory |
| `legacy/` | Historical amplicon, SOLiD, and merged-matrix branches |
| `docs/` | Steps, parameters, databases, contracts, and source records |

## Documentation

- [User guide](docs/user_guide.md): installation, commands, metrics, statistical designs, and troubleshooting.
- [Processing steps](docs/steps.md) and [workflow/checkpoints](docs/workflow.md).
- [Parameters](docs/parameters.md) and [reference databases](docs/databases.md).
- [Validation](validation/README.md), [legacy branches](legacy/README.md), and [source provenance](docs/provenance/README.md).

This module processes each dataset separately. Cross-study integration, drinking/disease label inference, differential abundance, functional prediction, subtyping, prediction models, and web deployment require separate analysis modules. Full study data and large references remain outside the repository; `examples/` is for verification.
