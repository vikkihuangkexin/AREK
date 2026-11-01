# AREK database code

This directory contains the code for constructing and analyzing the AREK database. Each numbered module represents a major processing stage and contains its workflow, algorithms, resources, documentation, and validation tools.

## Available modules

| Stage | Directory | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| 01 | [01_preprocessing](01_preprocessing/README.md) | Explicitly mapped FASTQ reads or feature/taxon matrices | Taxonomic abundance, diversity, statistics, figures, and execution records | Available |

Future modules can use names such as `02_<module>/` and `03_<module>/`. Functional prediction, differential abundance, disease modeling, and web services are outside the current module's scope.

## Module layout

| Location | Responsibility |
| --- | --- |
| `01_preprocessing/` | One complete database processing stage |
| `run_pipeline.py`, `README.md` | User entry point and quick start |
| `workflow/` | Configuration, input checks, external commands, logs, and checkpoints |
| `scripts/` | Computational modules used by the main workflow |
| `config/`, `envs/`, `examples/` | Parameters, dependency locks, and small verification inputs |
| `docs/`, `validation/` | Data contracts, source records, documentation, and executable checks |
| `legacy/` | Standalone historical branches with explicitly documented limits |

Specify output locations on the command line. Keep production datasets and large reference databases outside the repository. Downstream modules should consume verified outputs with explicit units, taxonomic ranks, sample identities, and provenance. OTU identifiers from different studies do not establish shared feature identity.

Database analysis code and future web application code can be maintained separately.
