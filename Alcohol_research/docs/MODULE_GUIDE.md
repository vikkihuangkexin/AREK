# Alcohol analysis — English module guide

## 1. Purpose and verified starting points

The release contains 22 independent command-line modules and their included shared implementations. It reconstructs the completed frozen analysis without new samples, covariates, thresholds or scientific conclusions. Start at external frozen profiles for modelling or original source compositions for source-level reconstruction. The original primary verification input remains in AC/demo/ in the supplied Biomni archive and is not included in this code-only release. Raw-read execution is implemented but **BLOCKED** in the validation environment.

The authoritative external archives are `Alcohol_M1_Completed_FULL.zip` and `Clinical_Interactions_and_Paper_PPT_FULL.zip`. They are not repackaged. Their original narrative, template and manifests remain immutable.

## 2. Installation, roots and common options

Use Python 3.11.13 and `requirements-python.lock.txt`. R-dependent beta/PAM require R 4.4.2, vegan 2.7.6 and cluster 2.1.8.3. Shell/R optional dependencies and unavailable original pins are listed separately. Native Windows supports the Python profile workflow; source Bash adapters need a compatible environment.

From the repository's `Alcohol_research` directory:

```bash
uv venv --python 3.11.13 .venv
uv pip install --python .venv/bin/python -r requirements-python.lock.txt
# Activate this environment before using the portable "python" commands below.
M1=/path/to/Alcohol_Independent_Comparisons_M1
python scripts/11_alcohol_vs_health_M1.py --input-root "$M1" --output-root run/main
```

On Windows use `.venv\Scripts\python.exe`. The external original AC/demo input only supports audit/profile/membership/main-M1 checks; use complete external roots for other analyses.

Set these example aliases to your actual extracted roots:

```bash
M1=/your/data/Alcohol_Independent_Comparisons_M1
CLIN=/your/data/Clinical_Interactions_and_Paper_PPT
V2=/your/data/Alcohol_4603_genus_function_v2
V3=/your/data/Alcohol_4603_pooled_vs_dataset_v3
DB=/your/data/FAPROTAX-mapper_v2.0
OUT=/your/writable/alcohol_rerun
```

Every entry point requires `--input-root` and `--output-root`; add `--clinical-root` where documented. Every successful node has `execution.json` with argument/input/code/output hashes and component detail. Output node names begin with the module ID, then applicable comparison/scope/panel/variable fields. Some legacy engines create nested paths below that isolated node; the receipt indexes them.

Common flags:

- `--help`: exact accepted options and defaults.
- `--dry-run`: check/list dependencies and target, no writes or deletions.
- `--resume`: reuse only identical argument/input/code/config hashes and intact outputs.
- `--overwrite`: replace only this module node; never silently overwrite.

The dependency driver adds only profile-level prerequisites:

```bash
python run.py --modules 11 --input-root "$M1" --output-root run/driver --dry-run
python run.py --modules 11 --input-root "$M1" --output-root run/driver
```

It executes `01 → 07 → 08 → 11`. Modules remain independent readers of frozen inputs: their outputs are not substituted into later models. Driver dry-run prints commands without certifying input availability. `--continue-on-error` records blocked dependants and permits independent work.

## 3. Exact comparison and scope keys

| Key | Frozen gut case / reference specimens |
|---|---:|
| `01_Alcohol_vs_Health` | 1117 / 2258 |
| `02_AUD_or_ALD_vs_Health` | 729 / 2258 |
| `03_AUD_vs_Health` | 453 / 2258 |
| `04_ALD_vs_Health` | 303 / 2258 |
| `05_ALD_vs_non_ALD` | 303 / 31 |

For group modules, scopes are `M1_all`, `M1_clinical_complete` and `supplemented`. Supplemented means the accepted-evidence age/sex/BMI complete subset, not a covariate-adjusted model. Interaction scopes instead are `frozen`, `exclude_explicit_disease_references`, or `all`.

The primary taxa panel is `gut__full_source_reference`. Beta/PAM default to `gut__common_all` (the shared-7 panel). Other exact panel IDs come from `M1/reproducibility/inputs/panels/`; do not invent panels. Clinical functions use body-level `gut`/`oral` panels.

## 4. Independent modules

### 01 — Audit inputs

- **Entry / implementation:** `scripts/01_audit_inputs.py`; `src/ac/early.py:audit`.
- **Inputs:** frozen metadata, comparison memberships, every panel TSS/CLR/family; optional original ZIPs through `--archives`.
- **Outputs:** `schema_checks.tsv`, `input_hashes.tsv`, `archives.json`, `stage_availability.json`.
- **Dependencies:** no earlier output; reads original profiles.
- **Validation:** schemas and actual hashes recomputed. Optional ZIP CRC/identity checks only run when archives are supplied. Raw-stage availability is not read-processing validation.

```bash
python scripts/01_audit_inputs.py --input-root "$M1" --output-root "$OUT"
```

### 02 — Collect and harmonize metadata

- **Entry / implementation:** `scripts/02_collect_harmonize_metadata.py`; `early.py:harmonize_metadata`.
- **Inputs:** clinical exact-ID evidence (`variable`, `sample_id`, `status`, `value`), frozen clinical values, technical metadata, unresolved conflicts and registry.
- **Outputs:** `normalized_clinical_values.tsv.gz`, `normalized_metadata.tsv.gz`, `exact_ID_provenance.tsv.gz`, `conflicts.tsv`, numerical checks.
- **Dependencies:** CLIN; no live metadata crawl.
- **Validation:** accepted individual values reconstructed and compared; technical metadata reused/hash-checked. Unresolved conflicts remain missing. No group-level fill or new evidence retrieval.

```bash
python scripts/02_collect_harmonize_metadata.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT"
```

### 03 — Preprocess original FASTQ

- **Entry / implementation:** `scripts/03_preprocess_fastq.py`; `raw.py` and the included recovered shell adapter.
- **Inputs:** source JSON, original project reads, original classifier or host/Kraken databases, tools.
- **Outputs:** `source_execution.json` and adapter-specific read-layout/QC or fastp products.
- **Dependencies:** source assets, not a downstream profile. Set every source output root inside `$OUT/03/`.
- **Validation:** Bash syntax and expected fail-closed preflight tested; actual read execution **BLOCKED**. Amplicon stops at layout/QC, shotgun at fastp.

```bash
python scripts/03_preprocess_fastq.py --input-root "$M1" --output-root "$OUT" --source-config source_qc.json
```

### 04 — Original 16S profiling

- **Entry / implementation:** `scripts/04_profile_16s.py`; `raw.py`, `amplicon_v5.sh`.
- **Inputs:** adapted `source_amplicon.example.json`, original reads, SILVA 138-99 classifier and QIIME 2/VSEARCH environment.
- **Outputs:** original adapter's mate-aware manifests/QC, 97% OTU table, taxonomy and source profiles under `$OUT/04/`.
- **Dependencies:** its own source QC, not an automatic handoff from 03. Preserve project-specific mate rules. No invented primer-trimming step.
- **Validation:** syntax/preflight only; raw profiling **BLOCKED**.

```bash
python scripts/04_profile_16s.py --input-root "$M1" --output-root "$OUT" --source-config source_16s.json
```

### 05 — Original shotgun profiling

- **Entry / implementation:** `scripts/05_profile_shotgun.py`; `raw.py`, `shotgun_v2.sh`.
- **Inputs:** adapted shotgun JSON, paired reads, hg38 Bowtie2 index, original Kraken2/Bracken DB, compatible tools.
- **Outputs:** fastp QC, host-depleted read processing, Kraken2/Bracken supported S/G/F/P profiles and source audits under `$OUT/05/`.
- **Dependencies:** source assets and its own preprocessing; output and missing-run roots must both be inside the node.
- **Validation:** syntax/preflight only; raw profiling **BLOCKED**. Original tool versions/database hashes not recovered.

```bash
python scripts/05_profile_shotgun.py --input-root "$M1" --output-root "$OUT" --source-config source_shotgun.json
```

### 06 — Build taxonomic abundance

- **Entry / implementation:** `scripts/06_build_taxonomic_abundance.py`; `abundance.py`.
- **Inputs:** by default frozen genus TSS/family. Optional `--profile` specimen-by-feature table and `--taxonomy` (`feature_id` plus requested rank), `--rank`, `--semantics`.
- **Outputs:** `genus_relative_abundance.tsv.gz` / `taxonomy.tsv` in frozen-export mode; supported-rank table, relative abundance, rank support and source taxonomy in aggregation mode.
- **Dependencies:** source-supported taxonomy only; no species inference or absent-rank reconstruction.
- **Validation:** frozen export exercised and hashed. Optional rank aggregation is implemented, not exhaustively rerun across sources.

```bash
python scripts/06_build_taxonomic_abundance.py --input-root "$M1" --output-root "$OUT"
# Optional original source table; counts/estimated_reads/fraction/percent must be declared correctly:
python scripts/06_build_taxonomic_abundance.py --input-root "$M1" --output-root run/rank --profile source_counts.tsv.gz --taxonomy source_taxonomy.tsv --rank family --semantics counts
```

### 07 — Harmonize profiles

- **Entry / implementation:** `scripts/07_harmonize_profiles.py`; `early.py:profiles`, shared CLR in `stats.py`.
- **Inputs:** frozen TSS/CLR/family; optional V3 source TSS/CLR and canonical crosswalk.
- **Outputs:** TSS, CLR, family, canonical taxa, structural measurability mask, optional rebuilt crosswalk.
- **Dependencies:** driver adds 01. Family is frozen, not reselected.
- **Validation:** source integration actually rebuilt and compared when V3 supplied. Otherwise canonical mapping/TSS reused and CLR recomputed/checked. Structural NA preserved.

```bash
python scripts/07_harmonize_profiles.py --input-root "$M1" --source-root "$V3" --output-root "$OUT"
```

### 08 — Define comparisons

- **Entry / implementation:** `scripts/08_define_comparisons.py`; `early.py:comparisons`.
- **Inputs:** frozen metadata/memberships; optional CLIN values for variable-specific availability.
- **Outputs:** frozen comparison membership, recomputed population counts, original three-field complete membership, optional per-variable IDs and explicit-disease references.
- **Dependencies:** driver adds 01.
- **Validation:** labels reused/hash-checked and counts recomputed; overlapping disease definitions retained.

```bash
python scripts/08_define_comparisons.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT"
```

### 09 — Alpha diversity

- **Entry / implementation:** `scripts/09_alpha_diversity.py`; `groups.py:alpha`, shared M1 engine.
- **Inputs:** panel abundance, stored alpha, metadata and memberships; V3 source compositions for fresh full-source measures.
- **Outputs:** `alpha_per_specimen.tsv.gz`, M1 results, membership, candidates and numerical comparison.
- **Dependencies:** driver adds 07/08. Full-source alpha is not the canonical pooled genus projection.
- **Validation:** primary gut full-source measures recomputed from source compositions and M1 reproduced. Without V3, measures explicitly reused/hash-checked and M1 freshly fitted.

```bash
python scripts/09_alpha_diversity.py --input-root "$M1" --source-root "$V3" --output-root "$OUT"
```

### 10 — Beta diversity

- **Entry / implementation:** `scripts/10_beta_diversity.py`; `beta.py`, `beta_kernel.py`, `beta_completion.R`; `supplemented.py` for supplemented scope.
- **Inputs:** TSS/CLR, metadata, memberships, frozen original node and original deterministic permutation definitions.
- **Outputs:** `completed_tests.tsv`, permutation F tables, R observed/dispersion results, receipt/session information and supporting matrices.
- **Dependencies:** driver adds 07/08; R/vegan. Does not recycle historical P values.
- **Validation:** primary gut shared-7, 9,999 permutations replayed; actual vegan marginal observed F/R² and first 99 permutations checked. Supplemented 16S pilot correctly gated at 47/2 working units (NA, zero permutations). Other nodes NOT_RUN.

```bash
python scripts/10_beta_diversity.py --input-root "$M1" --output-root "$OUT" --comparison 01_Alcohol_vs_Health --panel gut__common_all
python scripts/10_beta_diversity.py --input-root "$M1" --output-root run/beta_subset --scope supplemented --panel gut__common_16S
```

### 11 — Full Alcohol versus Health M1

- **Entry / implementation:** `scripts/11_alcohol_vs_health_M1.py`; `groups.py`, `m1.py`, `stats.py`.
- **Inputs:** frozen CLR, family, metadata, memberships and original result table.
- **Outputs:** `all_results.tsv.gz`, `candidates.tsv`, `membership.tsv.gz`, `numerical_comparison.tsv`.
- **Dependencies:** driver adds 07/08. Accepts only primary comparison and `M1_all`; use 12/13 otherwise.
- **Validation:** real-data full-reference replay: 3,201 rows; selected family 360, estimable 182, candidates 39 (15+/24−). Demo reproduces the same full-precision table.

```bash
python scripts/11_alcohol_vs_health_M1.py --input-root "$M1" --output-root "$OUT"
```

### 12 — Disease comparisons M1

- **Entry / implementation:** `scripts/12_disease_comparisons_M1.py`; shared group M1 engine.
- **Inputs / outputs:** same schemas as 11, but exact disease membership and its original selected-feature row scope.
- **Dependencies:** 07/08; select an explicit disease comparison.
- **Validation:** all four primary gut disease nodes replayed; candidate counts 33 / 14 / 22 / 15 in comparison-key order.

```bash
python scripts/12_disease_comparisons_M1.py --input-root "$M1" --output-root "$OUT" --comparison 02_AUD_or_ALD_vs_Health
python scripts/12_disease_comparisons_M1.py --input-root "$M1" --output-root "$OUT" --comparison 03_AUD_vs_Health
python scripts/12_disease_comparisons_M1.py --input-root "$M1" --output-root "$OUT" --comparison 04_ALD_vs_Health
python scripts/12_disease_comparisons_M1.py --input-root "$M1" --output-root "$OUT" --comparison 05_ALD_vs_non_ALD
```

### 13 — Clinical-complete subset M1

- **Entry / implementation:** `scripts/13_clinical_subset_M1.py`; shared group M1 engine.
- **Inputs:** original or supplemented three-field complete memberships, frozen profiles/family.
- **Outputs:** group M1 results/candidates/membership and full-precision checks.
- **Dependencies:** 07/08. These are subset M1 fits, not age/sex/BMI-adjusted or interaction fits.
- **Validation:** original and supplemented primary gut nodes replayed. Original primary subset retains full-reference rows; supplemented fits selected rows only. Supplemented candidates include Haemophilus and Veillonella, identified by their canonical keys in results.

```bash
python scripts/13_clinical_subset_M1.py --input-root "$M1" --output-root "$OUT" --scope M1_clinical_complete
python scripts/13_clinical_subset_M1.py --input-root "$M1" --output-root "$OUT" --scope supplemented
```

### 14 — Clinical main effects and interactions

- **Entry / implementation:** `scripts/14_clinical_interactions.py`; `clinical.py`, `interaction.py`, shared covariance.
- **Inputs:** exact-ID clinical values, technical metadata, profiles/diversity/functions and frozen registry.
- **Outputs:** combined overall tests/coefficients, node summary, membership; per-node coefficients, tests, category effects, support IDs/design manifests and within-Study fits.
- **Dependencies:** driver adds 01; CLIN. `--variable` accepts one exact registry key or `all`.
- **Validation:** full 406 nodes, 24,398 overall positions, 9,796 estimable, 120,801 coefficient rows compared. Full run has six primary variable–genus pairs, zero primary age candidates. A single-variable run withholds joint q-values when the frozen across-variable family is incomplete.

```bash
# All original nodes:
python scripts/14_clinical_interactions.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT" --variable all --panel all
# Independent variable:
python scripts/14_clinical_interactions.py --input-root "$M1" --clinical-root "$CLIN" --output-root run/sex --variable sex --body gut --scope frozen --outcome taxa --panel gut__full_source_reference --coverage base
```

### 15 — Source robustness

- **Entry / implementation:** `scripts/15_source_robustness.py`; `sources.py` and original group M1.
- **Inputs:** original fixed-family results, profile/functions, metadata and working-unit components.
- **Outputs:** `all_fixed_family_source_fits.tsv.gz`, `robustness_classification.tsv`, receipt.
- **Dependencies:** driver adds 01; taxa or `--outcome function`.
- **Validation:** primary gut taxa/functions replayed. Five priority taxa retain original canonical identities; no primary function meets source-direction robustness. Ordinary screen hits are not discarded.

```bash
python scripts/15_source_robustness.py --input-root "$M1" --output-root "$OUT"
python scripts/15_source_robustness.py --input-root "$M1" --output-root "$OUT" --outcome function
```

### 16 — Functional annotation and source contributions

- **Entry / implementation:** `scripts/16_functional_annotation.py`; `functions.py:annotation`.
- **Inputs:** V2 official mapper input/placements/states, genus reference audit, accepted 1/0/NA states, original source TSS and expected source contributions.
- **Outputs:** accepted genus functions, placement/mapping/unmatched audits, per-source contribution/known-mass/global-coverage tables and source checks.
- **Dependencies:** V2 required. Optional `--database-root "$DB" --rerun-mapper` invokes the original official mapper interface; not validated here.
- **Validation:** cached official states reused/hash-checked; source contributions freshly reconstructed and compared. Mapper execution NOT_RUN; missing packages/original pins documented.

```bash
python scripts/16_functional_annotation.py --input-root "$M1" --source-root "$V2" --output-root "$OUT"
# Optional, only after restoring the original mapper environment:
python scripts/16_functional_annotation.py --input-root "$M1" --source-root "$V2" --database-root "$DB" --rerun-mapper --output-root run/official_mapper
```

### 17 — Functional comparisons and clinical function models

- **Entry / implementation:** `scripts/17_functional_comparisons.py`; group M1 or shared clinical engine.
- **Inputs:** original-scale trait contributions, fixed trait family, metadata/memberships; global/trait coverage tables for sensitivities; CLIN for `--variable`.
- **Outputs:** group M1 results/candidates/checks or clinical function coefficients/tests/supports.
- **Dependencies:** driver adds 01. Trait contributions are inferred annotation, not measured function.
- **Validation:** primary base/global/trait-coverage fits replayed; nine base candidates, eight in both sensitivity screens. Independent gut sex/function base run also exercised (one node, 89 positions, 59 estimable, no candidates); joint q withheld for incomplete variable family.

```bash
python scripts/17_functional_comparisons.py --input-root "$M1" --output-root "$OUT" --coverage base
python scripts/17_functional_comparisons.py --input-root "$M1" --output-root "$OUT" --coverage global_coverage
python scripts/17_functional_comparisons.py --input-root "$M1" --output-root "$OUT" --coverage trait_coverage
python scripts/17_functional_comparisons.py --input-root "$M1" --clinical-root "$CLIN" --output-root run/function_sex --variable sex --panel gut --body gut --clinical-scope frozen
```

### 18 — Taxon-to-function contributions

- **Entry / implementation:** `scripts/18_taxon_function_contributions.py`; `functions.py:contributions`.
- **Inputs:** candidate taxon/trait mapping, accepted states, source-original TSS, original per-ID expected contributions, memberships.
- **Outputs:** `contributions_per_ID.tsv.gz`, `group_summary.tsv`.
- **Dependencies:** driver adds 01.
- **Validation:** candidate contributions freshly summed and compared at full precision before group/scope filtering. Preserve original TSS scale and overlapping traits; no cross-trait closure.

```bash
python scripts/18_taxon_function_contributions.py --input-root "$M1" --output-root "$OUT"
```

### 19 — Patient patterns

- **Entry / implementation:** `scripts/19_patient_patterns.py`; `patterns.py`, `pattern_completion.R`, supplemented adapter.
- **Inputs:** case CLR, working-cluster metadata, original node definitions and pinned original seeds.
- **Outputs:** learning metadata/CLR, PCA scores/loadings, K-specific labels/metrics, bootstrap Jaccards and assignment stability, source contingency tables, completion receipt.
- **Dependencies:** 07/08; R/cluster. Full CLR drives PAM, not the PCA display.
- **Validation:** ALD shared-7 accepted K=2, 298 units, sizes 219/79, ARI=1 versus original; Alcohol 1,004 units and no accepted K. Supplemented ALD has 18 units and correctly returns not-estimable without bootstrapping.

```bash
python scripts/19_patient_patterns.py --input-root "$M1" --output-root "$OUT" --comparison 04_ALD_vs_Health
python scripts/19_patient_patterns.py --input-root "$M1" --output-root "$OUT" --comparison 01_Alcohol_vs_Health
python scripts/19_patient_patterns.py --input-root "$M1" --output-root "$OUT" --comparison 04_ALD_vs_Health --scope supplemented
```

### 20 — Candidate overlaps

- **Entry / implementation:** `scripts/20_candidate_overlap.py`; `overlap.py`.
- **Inputs:** original full/clinical-complete/supplemented M1 result tables, optional clinical tests, canonical IDs.
- **Outputs:** canonical overlaps, counts, optional interaction overlap, `disease_overlap.svg` and plot data.
- **Dependencies:** driver adds 01; original model outputs are consumed, not refitted.
- **Validation:** canonical set comparisons recomputed. Same-direction primary/disease overlaps 28 / 10 / 15 / 13; display genus names never merge homonyms.

```bash
python scripts/20_candidate_overlap.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT"
```

### 21 — Figures

- **Entry / implementation:** `scripts/21_make_figures.py`; `figures.py`.
- **Inputs:** completed result tables, beta/PAM summaries, optional clinical tests.
- **Outputs:** eight principal editable-text SVGs when full authorities are supplied, corresponding `.data.tsv.gz` tables and `figure_index.tsv`.
- **Dependencies:** driver adds 01; plotting consumes result values, not a new model.
- **Validation:** plots/data index regenerated from actual values. No simulated points, new statistics, narrative overlays or generic schematic substitutes.

```bash
python scripts/21_make_figures.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT"
```

### 22 — Presentation

- **Entry / implementation:** `scripts/22_build_presentation.py`; `presentation.py`.
- **Inputs:** CLIN's original packed `.deck.html`, source index, local figures/fonts and assets. These remain external immutable inputs.
- **Outputs:** isolated rebuilt `Alcohol_Microbiome_Main_Results.deck.html`, `asset_provenance.tsv`, identity receipt.
- **Dependencies:** driver adds 01; no narrative rewrite or newly authored slide conclusions.
- **Validation:** 16 slides, original narrative/slide IDs/runtime and 565 editable IDs retained, assets re-embedded. Standalone PPTX/PDF not generated; use the presentation editor for any later export.

```bash
python scripts/22_build_presentation.py --input-root "$M1" --clinical-root "$CLIN" --output-root "$OUT"
```

## 5. Statistical and data-integrity rules

- Preserve exact specimen and canonical feature IDs, structural NA, frozen feature families and original scope.
- M1 retains the identifiable original source nuisance span. HC3/CR1, working-unit gates, residual df and cluster thresholds are unchanged.
- Fixed-family BH counts non-estimable positions in its denominator while leaving their displayed q missing. Within-variable, joint-across-variable and full-reference q columns are distinct.
- Primary/original primary complete-case taxa tables retain 3,201 reference rows; disease and supplemented outputs use selected rows only. Do not add unselected disease tests.
- Clinical subset M1, clinical-main and interaction fits share support within each variable/feature. Do not interpret a subset fit as clinical adjustment.
- Full-source alpha needs original source compositions for fresh measures. FAPROTAX is annotation inference; contributions stay on original nonexclusive TSS scale.
- A gate failure is a legitimate missing result, not a negative fitted effect. Source robustness, dispersion and patient-pattern gates do not imply causality or external replication.

Full definitions and input schemas are in `docs/ANALYSIS_DEFINITIONS.md` and `docs/INPUT_REQUIREMENTS.md` inside the code archive.

## 6. Validation and limitations

Recorded checks include all 22 syntax/import/help interfaces, Python tests against independent BH/statsmodels CR1, shell/R parsing, overwrite/resume/dry-run/preflight safeguards, 73 full-precision/identity checks and the full 406-node clinical replay. Actual 9,999-permutation beta and ALD/Alcohol PAM replay were run. Specific supplemented nodes correctly failed their original estimability gates.

Validation **does not mean all 22 modules were fully recomputed from raw reads**. Raw modules 03–05 are BLOCKED. Official mapper invocation and unlisted beta/PAM/rank branches remain NOT_RUN. Frozen labels/families, official annotation states and presentation narrative/runtime are intentionally reused. The coverage table explicitly separates fresh calculations from reused assets.

The GitHub increment contains a research FILE_MANIFEST.tsv, root README/ignore updates, and Windows-safe paths. Full data are external. Archive validation checks independent extraction and hashes, then can replay primary M1 when an external verification input is supplied. Its receipt remains outside the ZIP to avoid a circular self-hash.
