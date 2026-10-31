# AREK

Research code for the **AREK database**. Database construction and analysis code is maintained under [`AREK_database/`](AREK_database/README.md), with one numbered directory per major stage.

| Stage | Module | Scope |
| --- | --- | --- |
| 01 | [Preprocessing](AREK_database/01_preprocessing/README.md) | FASTQ discovery and processing, taxonomic abundance tables, diversity and figures |

The first module supports amplicon OTU97, shotgun Kraken2/Bracken, and existing-matrix inputs. Its canonical workflow, standalone historical branches, configuration, environments and validation programs are documented separately. Later modules will be added alongside `01_preprocessing/` when their code is available.

中文使用入口：[数据库代码目录](AREK_database/README.md) · [预处理快速开始](AREK_database/01_preprocessing/README.md) · [步骤与源码对应](AREK_database/01_preprocessing/docs/steps_zh.md)。

Large reference databases and full study datasets are supplied externally. Small public read fixtures are included for verification. Historical logs, result snapshots and cohort inventories remain in the original source archive and are excluded from this repository.

License: [MIT](LICENSE). Third-party BFAST source retains its own bundled license.
