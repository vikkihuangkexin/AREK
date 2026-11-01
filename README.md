# AREK

Code for **AREK database construction** and **alcohol-related microbiome research**, maintained as separate areas of this repository.

| Area | Directory | Purpose |
| --- | --- | --- |
| Database construction | [AREK_database](AREK_database/README.md) | Reusable processing stages for building the AREK database |
| Alcohol research | [Alcohol_research](Alcohol_research/README.md) | Frozen study comparisons, clinical interactions, source robustness, inferred functions, patient patterns, and figures |

## Database construction

The available [preprocessing module](AREK_database/01_preprocessing/README.md) covers FASTQ discovery and processing, taxonomic abundance, diversity, statistics, and figures. It supports amplicon OTU97, shotgun Kraken2/Bracken, and existing-matrix inputs with explicit sample mapping, input checks, content checkpoints, and output validation.

Start with the [database user guide](AREK_database/01_preprocessing/docs/user_guide.md), [processing steps](AREK_database/01_preprocessing/docs/steps.md), [parameters](AREK_database/01_preprocessing/docs/parameters.md), [reference databases](AREK_database/01_preprocessing/docs/databases.md), and [validation notes](AREK_database/01_preprocessing/validation/README.md).

## Alcohol microbiome research

The [Alcohol research package](Alcohol_research/README.md) contains 22 command-line modules with shared Python, R, and shell implementations. Its main scientific entry points consume the frozen analysis profiles and accepted clinical evidence. Dataset definitions, feature families, statistical support gates, multiplicity rules, and original result authorities are retained.

| Analysis | Modules | Main purpose |
| --- | --- | --- |
| Inputs and metadata | 01–02, 06–08 | Validate inputs, replay accepted metadata, export supported abundance ranks, align profiles, and retain comparison membership |
| Diversity and group comparisons | 09–13 | Alpha/Beta diversity, Alcohol-versus-reference M1, disease comparisons, and clinical-complete subset M1 |
| Clinical and source analyses | 14–15 | Clinical main effects/interactions and within-Study/leave-source-out robustness |
| Inferred functions | 16–18 | FAPROTAX annotation replay, functional comparisons, coverage sensitivities, and taxon contributions |
| Patterns and reporting | 19–22 | PAM patient patterns, candidate overlap, figures, and existing HTML presentation reconstruction |
| Optional historical adapters | 03–05 | Isolated source-specific FASTQ processing; original raw-read execution remains unvalidated |

The research package starts at verified frozen profiles for the main comparisons. Its historical raw adapters are separate from the database preprocessing module and are not automatically run by the main analysis driver. Patient patterns use unsupervised PAM; supervised Alcohol/AUD/ALD multiclass prediction is not implemented in this release.

Start with the [research quick start](Alcohol_research/README.md), [22-module guide](Alcohol_research/docs/MODULE_GUIDE.md), [input requirements](Alcohol_research/docs/INPUT_REQUIREMENTS.md), [analysis definitions](Alcohol_research/docs/ANALYSIS_DEFINITIONS.md), and [validation boundaries](Alcohol_research/validation/README.md).

## Data and environments

Production cohorts, full metadata, raw reads, large reference databases, detailed run logs, and presentation assets are supplied outside the repository. Database preprocessing includes small public verification fixtures. Research code includes configuration, source provenance, tests, and compact historical validation summaries; its full-data demo remains in the original Biomni archive.

Each area has its own dependency records. Preserve units, taxonomic identities, sample mappings, host/body-site boundaries, and source provenance when moving data between them; the research package does not directly substitute arbitrary database outputs for its frozen inputs.

## License

[MIT](LICENSE). Bundled third-party BFAST source in the database area retains its own license.
