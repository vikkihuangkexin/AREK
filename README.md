# AREK

Research code for the **AREK database**, organized by database construction and analysis stage.

| Stage | Module | Scope |
| --- | --- | --- |
| 01 | [Preprocessing](AREK_database/01_preprocessing/README.md) | FASTQ discovery and processing, taxonomic abundance tables, diversity, statistics, and figures |

The preprocessing module supports paired-end and single-end amplicon reads, shotgun metagenomic reads, and existing abundance matrices. It provides explicit sample mapping, configuration checks, execution plans, content-based checkpoints, and output validation.

## Getting started

1. Open the [preprocessing quick start](AREK_database/01_preprocessing/README.md).
2. Select the appropriate route and restore its Linux environment.
3. Supply the required inputs and external reference databases.
4. Run `check`, inspect `dry-run`, execute `run`, and verify the results with `validate`.

| Resource | Description |
| --- | --- |
| [Database code organization](AREK_database/README.md) | Module layout and interfaces |
| [User guide](AREK_database/01_preprocessing/docs/user_guide.md) | Installation, commands, input requirements, metrics, and limitations |
| [Processing steps](AREK_database/01_preprocessing/docs/steps.md) | Steps, source modules, dependencies, and outputs |
| [Parameters](AREK_database/01_preprocessing/docs/parameters.md) | Defaults, parameter evidence, and method differences |
| [Reference databases](AREK_database/01_preprocessing/docs/databases.md) | Required resources and reproducibility boundaries |
| [Validation](AREK_database/01_preprocessing/validation/README.md) | Reproducible checks and their scope |

Large reference databases and full study datasets are supplied externally. Small public read subsets and matrix fixtures are included for software verification. Historical logs, result snapshots, and cohort inventories remain in the original source archive.

Later analysis modules will be added alongside `01_preprocessing/` when their code is available.

## License

[MIT](LICENSE). The bundled third-party BFAST source retains its own license.
