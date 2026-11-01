# Frozen analysis definitions

This release changes code organization and English documentation, not scientific definitions. IDs, missingness, families, models, gates and original narrative remain frozen.

## Populations and scope

- 4,603 specimens from 28 Studies: gut 3,375 (1,117 Alcohol / 2,258 reference), oral 1,228 (958 / 270). Gut contains 475 16S and 2,900 shotgun specimens; oral is 16S.
- Gut comparison populations: AUD-or-ALD 729 / 2,258; AUD 453 / 2,258; ALD 303 / 2,258; ALD versus non-ALD 303 / 31. AUD/ALD overlap is 27, not an independent sixth group.
- `M1_all`: frozen full eligible population. `M1_clinical_complete`: original age/sex/BMI complete-case subset. `supplemented`: accepted evidence-supplemented three-field complete subset, 385 / 2,172 for primary gut comparison. These subset M1 models are **not clinically adjusted and are not interaction models**.
- Clinical interactions use per-variable, per-feature support, not the three-field complete-case subset. Gut availability for age and sex is 425 / 2,226 and BMI 386 / 2,172; oral age 774 / 270, sex 815 / 270, BMI 754 / 262. Estimable feature support can be smaller.
- Exclusion sensitivity removes only 92 explicitly disease-labelled gut reference specimens (retaining 1,117 / 2,166). It neither changes the frozen main model nor certifies every remaining reference as healthy.

## Taxonomy, abundance and missingness

Family-anchored IDs such as `FG::Bacteria::<family>::<genus>` are distinct from ambiguous source-specific `SRC::...` IDs. Homonymous genera are not collapsed by display name. A confident mapping requires the original genus resolution and acceptable family anchor. Supported ranks are aggregated only from source taxonomy.

TSS is within the original source composition for full-source panels. Structural unmeasurability stays `NA`, never a biological zero. Zeros in supported compositions are replaced with the original `1e-6` pseudocount, closed and log-centred. Common panels are closed within that panel; full-source CLR retains its original source coordinate system.

Full-source alpha is Shannon, Gini-Simpson and observed genera from the **individual source-specific composition**, closed per specimen; it is not computed from the canonical pooled projection. Without `V3`, module 09 explicitly reuses/hash-checks stored measures and freshly fits their M1 models. Common-panel alpha uses the closed common panel.

## Group M1

The exposure effect is Alcohol/case versus reference in an OLS model with the identifiable span of Study, analysis set, assay, material and technical batch (batch nested in Study). Ordered rank reduction preserves the original nuisance span; no alternative covariate selection is introduced.

Per-feature support and frozen working clusters are used. Require at least five working units per arm, residual degrees of freedom at least 10, identifiable exposure; clustered CR1 requires at least 20 clusters. HC3 applies when the original independent-unit rule permits; otherwise CR1 uses the original finite-sample correction and reference degrees of freedom. Unsupported fits remain missing with an explicit reason.

BH uses the fixed selected family, counting non-estimable positions as `p=1` for the denominator but displaying their q-values as missing. Exploratory screen: `q<0.20`. Full-reference BH is a separate column/family.

**Original row-scope distinction:** primary full M1 and original primary clinical-complete tables retain all 3,201 reference rows with 360 selected family positions. Disease and supplemented clinical-complete taxa tables fit/export selected rows only (360 in the gut primary panel). Do not “standardize” these into new tests or silently change the full-reference denominator.

## Clinical models and multiplicity

For each variable/feature, fit subset M1, clinical-main-effect and exposure-by-clinical-interaction models on the **same feature support**. Age is centred at 50 and scaled per 10 years; BMI at 25 per 5 kg/m²; sex reference is female. Original category encodings and reference levels remain in `configs/clinical_registry.json`.

Continuous variables require shared-range support of at least five working units and five distinct values per arm. This gate does not truncate the fitted data to the overlap range. Categorical cells require five working units. Multicategory interaction is a whole-block robust Wald F test, not a selected single coefficient.

Within-variable fixed-family BH and joint-across-variable BH are separate. A partial `--variable sex` run intentionally withholds joint q-values unless the full frozen variable family is present for that scope/body/outcome/panel/coverage. The full run reproduces 406 nodes, 24,398 overall positions and 9,796 estimable tests. The six primary findings are variable–genus pairs, not six independently replicated taxa. Primary age interactions have zero candidates.

Accepted exact-ID evidence is replayed without group-level imputation. The 52 perturbed HepatAlc age/sex records are not genuine individual fills; 235 bilirubin conflicts remain missing. CTP is a clinical score; the source's creatine-labelled field means creatinine. Smoking definitions remain separate. The original 1,099-column source-key dictionary is retained for audit; review-only fields are not promoted into models.

## Source robustness

Use the original within-Study and leave-one-Study/component-out fits. Shared working units define connected source components. Direction/estimability classifications supplement, not delete, the ordinary primary screen. A candidate supported by one estimable Study is not multi-Study replication.

## Functions

FAPROTAX represents inferred annotation, not a measured pathway or functional assay. Accept original genus-level placements only; preserve 1 / 0 / unknown states. Contributions are sums of original TSS abundance times accepted state. Traits may overlap; never close across traits. A trait with no known mass is missing, not zero. Global-coverage and trait-known-mass sensitivities preserve their original covariates and families.

The primary function screen has nine candidates; eight pass the intersection of both coverage sensitivities; zero meet the original source-direction robustness rule. Module 16 replays official cached mapper states and freshly reconstructs source contributions unless `--rerun-mapper` is explicitly requested.

## Beta diversity

Use one deterministic representative per working unit, original nuisance span, 9,999 restricted permutations in frozen Study/analysis-set/assay/material blocks, Bray-Curtis and Aitchison. Corrected observed F/R² are **vegan marginal** quantities; permutation kernels are cross-checked against the first 99 vegan permutations. The commuting-kernel route uses inverse permutation indexing. General residual kernels are used otherwise.

Dispersion uses spatial median and bias adjustment, with Lingoes correction for Bray. Significant dispersion prevents interpreting a location test alone as an unqualified location shift. Gate failures remain NA with zero permutations, not a nominal “completed” 9,999-permutation result.

## Patient patterns

PAM uses case working-unit representatives, full CLR geometry and pinned original seeds, not seeds recomputed from translated code hashes. K=2–6, 200 bootstrap attempts; require at least 190 valid bootstraps, minimum cluster size 20 and fraction 0.05, silhouette at least 0.25, minimum mean Jaccard at least 0.75, and original source-NMI exclusion at 0.8. Select the smallest passing K within 0.01 of the maximum passing silhouette. PCA is display geometry, not the clustering input.

ALD shared-7: 298 learning units, accepted K=2 sizes 219 / 79; reproduced partition ARI=1. Alcohol shared-7: 1,004 units, no accepted K. Supplemented ALD shared-7: 18 learning units, below the 40-unit gate; no bootstrap fit was performed.

## Presentation boundary

Module 21 plots existing/recomputed result values into eight principal SVGs with matching data tables; it does not refit statistics. Module 22 reconstructs the existing 16-slide self-contained HTML, preserving the original Chinese narrative, runtime, slide IDs, 565 editable IDs and external assets. It is not a new English narrative or a standalone PPTX/PDF export.
