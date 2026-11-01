# Alcohol-related microbiome research

This directory contains the modular code recovered by Biomni for the completed, frozen alcohol-related microbiome analyses. It is separate from `AREK_database/`. The 22 entry points call implementations included in `src/ac/`; the internal Python namespace remains `ac` after the repository directory was renamed.

The release preserves the scientific definitions and numerical model implementations. It does not add samples, reselect features, adjust thresholds, or replace the frozen result authorities. Full cohorts, metadata, large inputs, detailed execution logs, and presentation assets are external.

## Modules at a glance

| Modules | Task | Main outputs |
| --- | --- | --- |
| 01–02 | Input audit and accepted clinical metadata replay | Schema/hash checks, harmonized values, provenance, unresolved conflicts |
| 03–05 | Optional historical FASTQ adapters | Source-specific QC/profiling; raw execution remains blocked without original tools and references |
| 06–08 | Supported rank export, profile alignment, and frozen memberships | Abundance tables, canonical mappings/CLR checks, comparison counts |
| 09–10 | Alpha and Beta diversity | Diversity measures, M1 Alpha tests, restricted PERMANOVA and dispersion |
| 11–13 | Alcohol and disease M1 comparisons, complete-case subset M1 | Effects, robust uncertainty, fixed-family q-values, candidates, memberships |
| 14 | Clinical main effects and exposure-by-clinical interactions | Per-variable/per-feature support, coefficients, robust Wald tests, multiplicity families |
| 15 | Within-Study and leave-source-out robustness | Source/component fits, direction checks, robustness classifications |
| 16–18 | Inferred functions and taxon contributions | Accepted annotation states, original-scale contributions, functional M1 and coverage sensitivities |
| 19 | Unsupervised patient patterns | PAM K=2–6, stability/size/source gates, assignments, display PCA |
| 20–22 | Candidate overlap, figures, and existing presentation | Overlap tables, eight principal SVG/data pairs, reconstructed HTML |

Module 19 does not train a supervised multiclass classifier. Module 22 preserves the external presentation narrative and does not translate it or create a standalone PPTX/PDF. The optional raw adapters are not replacements for the separate database preprocessing module.

## Quick start

Biomni's reference validation used Python 3.11.13 and the exact package versions in `requirements-python.lock.txt`. From the repository root:

```bash
cd Alcohol_research
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-python.lock.txt

# Point to the contained package directory after extracting the original M1 archive.
M1=/path/to/Alcohol_Independent_Comparisons_M1
OUT=/path/to/writable/alcohol_outputs

python scripts/11_alcohol_vs_health_M1.py --help
python run.py --modules 11 --input-root "$M1" --output-root "$OUT" --dry-run
python run.py --modules 11 --input-root "$M1" --output-root "$OUT"
```

On Windows use `.venv\Scripts\python.exe`; Bash source adapters need Linux or a compatible environment. The driver orders `01,07,08,11`. These modules independently read frozen inputs; it does not replace later inputs with newly generated upstream outputs or launch raw-read profiling.

The primary full-source gut comparison is 1,117 Alcohol versus 2,258 reference specimens. The original M1 reference has 3,201 feature rows, 360 selected family positions, 182 estimable selected positions, and 39 exploratory candidates at q<0.20 (15 positive, 24 negative). Those are frozen-profile regression results, not prediction accuracy. The full-data verification input remains under `AC/demo/` in the original Biomni code archive and can be supplied as an external `--input-root` for modules 01, 07, 08, and 11.

## External roots

| Option or alias | Required content |
| --- | --- |
| `--input-root` / M1 | Extracted `Alcohol_Independent_Comparisons_M1`: frozen profiles, memberships, families, completed results, and evidence-supplemented inputs |
| `--clinical-root` / CLIN | Extracted `Clinical_Interactions_and_Paper_PPT`: accepted exact-ID evidence, clinical model inputs/results, and presentation source |
| `--source-root` / V3 | Source compositions and canonical crosswalk for profile/Alpha reconstruction |
| `--source-root` / V2 | Original genus/function inputs and accepted FAPROTAX annotation resources for module 16 |
| `--database-root` / DB | Original FAPROTAX mapper and reference assets for optional mapper re-execution |

Use the exact schemas in [INPUT_REQUIREMENTS.md](docs/INPUT_REQUIREMENTS.md). The package is a frozen-study replay, not a general interface for arbitrary new cohorts. Database output requires an explicitly verified transformation before it can satisfy these inputs.

```bash
CLIN=/path/to/Clinical_Interactions_and_Paper_PPT
python scripts/12_disease_comparisons_M1.py \
  --input-root "$M1" --output-root "$OUT" --comparison 04_ALD_vs_Health
python scripts/14_clinical_interactions.py \
  --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT" \
  --variable all --panel all
python scripts/19_patient_patterns.py \
  --input-root "$M1" --output-root "$OUT" \
  --comparison 04_ALD_vs_Health --panel gut__common_all
```

Beta/PAM require the separately recorded R environment. Modules 03–05 require original reads, matching references, and Bash tool environments. Official FAPROTAX mapper execution has additional unpinned optional dependencies and was not validated by Biomni.

## Execution and validation

- Existing nodes are refused unless `--resume` matches arguments/input/code/configuration hashes and intact outputs, or explicit `--overwrite` is requested.
- Output roots must be outside authoritative input roots. Module-level `--dry-run` checks dependencies without writing; driver-level `--dry-run` only prints ordered commands.
- `execution.json` records successful execution and hashes. Its generic `RECOMPUTED` label does not mean every component was recalculated: inspect `detail` and the documented reused inputs.
- Valid `not_estimable`, `BLOCKED`, and `NOT_RUN` states retain their meaning. An unavailable analysis is not a validation pass.

[validation/README.md](validation/README.md) distinguishes supplied Biomni evidence from checks performed after this repository reorganization. Full clinical, R Beta/PAM, functional, and raw-read chains are not certified by the packaging checks.

```bash
# Data-dependent interface checks use an external primary verification input.
python validation/check_interfaces.py \
  --input-root /path/to/original/AC/demo --output-root /path/to/new/interface_checks

# Independent statistical tests require the external input and test dependencies.
AC_TEST_INPUT_ROOT=/path/to/original/AC/demo python -m pytest -q tests
```

## Documentation and layout

| Path | Purpose |
| --- | --- |
| [CODE_INDEX.tsv](CODE_INDEX.tsv) | All 22 entries, implementations, dependencies, and original validation boundaries |
| [docs/MODULE_GUIDE.md](docs/MODULE_GUIDE.md) | Inputs, outputs, commands, and limits for each module |
| [docs/ANALYSIS_DEFINITIONS.md](docs/ANALYSIS_DEFINITIONS.md) | Populations, M1/clinical estimators, support, multiplicity, functions, and PAM |
| [docs/INPUT_REQUIREMENTS.md](docs/INPUT_REQUIREMENTS.md) | External file schemas and environments |
| [docs/KNOWN_GAPS.md](docs/KNOWN_GAPS.md) | Blocked, reused, partial, and untested branches |
| `scripts/`, `src/ac/`, `configs/` | Entry points, included implementations, and frozen configuration |
| [history/README.md](history/README.md) | Original source identities and repository preparation provenance |
| [validation/README.md](validation/README.md) | Original compact summaries and current package checks |
| `FILE_MANIFEST.tsv` | Current research-module paths, sizes, and SHA256; excludes itself |

All repository-authored code, help, comments, filenames, and documentation are English. Source field keys, category coding, original hashes, and clinical definitions remain unchanged.
