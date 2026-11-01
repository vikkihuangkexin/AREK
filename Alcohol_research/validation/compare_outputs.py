"""Compare recomputed outputs against full-precision authoritative tables."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE/"src"))
from ac.io import read, save, dump
from ac.stats import bh

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", required=True)
    ap.add_argument("--clinical-root", required=True)
    ap.add_argument("--validation-root", required=True)
    a = ap.parse_args()
    old, clinical, root = map(Path, [a.input_root, a.clinical_root, a.validation_root])
    checks = []
    def check(name, value, **detail):
        if not bool(value):
            raise AssertionError(name+" "+str(detail))
        checks.append(dict(check=name, status="PASS", **detail))
    def compare(name, current, frozen, columns, keys, rtol=2e-7, atol=2e-9):
        x = current.set_index(keys).sort_index()
        y = frozen.set_index(keys).reindex(x.index)
        check(name+"_row_identity", not x.index.duplicated().any() and len(x) == len(frozen) and set(x.index)==set(frozen.set_index(keys).index))
        for col in columns:
            v, w = pd.to_numeric(x[col]).to_numpy(), pd.to_numeric(y[col]).to_numpy()
            q = col.startswith("q") or col == "BH_q"
            rt, at = (1e-10, 1e-12) if q else (rtol, atol)
            np.testing.assert_allclose(v, w, equal_nan=True, rtol=rt, atol=at, err_msg=name+":"+col)
            d = abs(v-w)
            checks.append(dict(check=name+":"+col, status="PASS", rows=len(x), rtol=rt, atol=at,
                max_absolute_error=float(np.nanmax(d)) if np.isfinite(d).any() else None))
        return x,y
    full = next((root/"clinical_full").glob("*/all_interaction_tests.tsv.gz"))
    tests = read(full)
    frozen_tests = read(clinical/"interaction_models/all_interaction_tests.tsv.gz")
    keys = ["scope", "body", "variable", "outcome", "panel", "coverage", "feature_id"]
    x,y = compare("clinical_all_tests", tests, frozen_tests,
        ["P_interaction", "Wald_F", "q_interaction", "q_interaction_joint_variables",
         "fixed_family_size", "joint_fixed_family_size", "n_Alcohol", "n_Health",
         "working_units_Alcohol", "working_units_Health", "n_Study"], keys)
    for col in ["sample_hash", "working_unit_hash", "reason"]:
        check("clinical_identity_"+col, x[col].fillna("").equals(y[col].fillna("")))
    co = read(full.parent/"all_model_coefficients.tsv.gz")
    expected = read(clinical/"interaction_models/all_model_coefficients.tsv.gz")
    compare("clinical_all_coefficients", co, expected,
        ["beta", "SE", "CI_low", "CI_high", "raw_P", "q_Alcohol"], keys+["model_id", "term"])
    check("clinical_nodes", len(read(full.parent/"node_summary.tsv")) == 406, expected=406)
    check("clinical_test_positions", len(tests) == 24398, expected=24398)
    check("clinical_estimable", tests.P_interaction.notna().sum() == 9796, expected=9796)
    primary = tests[tests.scope.eq("frozen")&tests.body.eq("gut")&tests.outcome.eq("taxa")
        &tests.panel.eq("gut__full_source_reference")&tests.coverage.eq("base")]
    check("primary_six_variable_genus_pairs", primary.q_interaction.lt(.2).sum() == 6, expected=6)
    check("primary_age_no_candidates", not primary[primary.variable.eq("age")].q_interaction.lt(.2).any())
    joint_keys = ["scope", "body", "outcome", "panel", "coverage"]
    maximum = 0.
    for _, g in tests.groupby(joint_keys):
        d = abs(bh(g.P_interaction)-g.q_interaction_joint_variables)
        if np.isfinite(d).any():
            maximum = max(maximum, float(np.nanmax(d)))
    check("clinical_independent_joint_BH", maximum < 1e-12, max_absolute_error=maximum)
    single = read(next((root/"clinical_single").glob("*/all_interaction_tests.tsv.gz")))
    check("single_variable_joint_BH_not_mislabelled", single.q_interaction_joint_variables.isna().all()
        and single.joint_BH_status.eq("incomplete_variable_family").all())
    for label, number in [("main",39),("disease_1",33),("disease_2",14),("disease_3",22),("disease_4",15)]:
        receipt = json.loads(next((root/label).glob("*/execution.json")).read_text())["detail"]
        check(label+"_candidate_count", receipt["candidates"] == number, expected=number, observed=receipt["candidates"])
    source = read(next((root/"source_main").rglob("robustness_classification.tsv")))
    priorities = source[source.source_robust&source.ordinary_primary_hit]
    check("five_source_direction_priorities", len(priorities) == 5, observed=priorities.taxon_key.tolist())
    source_expected = read(old/"01_Alcohol_vs_Health/gut/M1_all/source_robustness/taxa/gut__full_source_reference/robustness_classification.tsv")
    compare("source_classification", source, source_expected,
        ["primary_beta","primary_raw_p","primary_BH_q","n_estimable_two_arm_Studies",
         "n_independent_study_components","n_valid_LOSO"], ["taxon_key"])
    for col in ["source_robust", "classification"]:
        check("source_classification_"+col,
            source.set_index("taxon_key")[col].sort_index().equals(source_expected.set_index("taxon_key")[col].sort_index()))
    fs = read(next((root/"source_function").rglob("robustness_classification.tsv")))
    check("function_source_robust_zero", (fs.source_robust&fs.ordinary_primary_hit).sum() == 0)
    candidates = []
    for label in ["function_base","function_global","function_trait"]:
        d = read(next((root/label).glob("*/all_results.tsv.gz")))
        candidates.append(set(d.loc[d.BH_q.lt(.2),"taxon_key"]))
    check("function_candidates_nine", len(candidates[0]) == 9, expected=9)
    check("function_retained_both_eight", len(set.intersection(*candidates)) == 8, expected=8)
    beta = read(next((root/"beta_shared").rglob("completed_tests.tsv")))
    frozen = read(old/"beta_completed_summary.tsv")
    frozen = frozen[frozen.comparison.eq("01_Alcohol_vs_Health")&frozen.scope.eq("M1_all")&frozen.panel_id.eq("gut__common_all")]
    compare("corrected_beta", beta, frozen, ["F","R2","P","dispersion_F","dispersion_P"], ["metric"], rtol=1e-7, atol=1e-7)
    check("beta_all_permutations", beta.B_actual.eq(9999).all())
    for label, comparison in [("PAM_ALD","04_ALD_vs_Health"),("PAM_Alcohol","01_Alcohol_vs_Health")]:
        new = next((root/label).rglob("PAM_completion"))
        original = old/comparison/"gut/M1_all/patterns/gut__common_all/PAM_completion"
        rec = json.loads((new/"completion_receipt.json").read_text())
        frozen = json.loads((original/"completion_receipt.json").read_text())
        check(label+"_status", rec["status"] == frozen["status"], observed=rec["status"])
        check(label+"_seed", rec["dependency"]["seed"] == frozen["dependency"]["seed"], observed=rec["dependency"]["seed"])
        compare(label+"_all_K_metrics", read(new/"all_K_metrics.tsv"), read(original/"all_K_metrics.tsv"),
            ["silhouette","min_mean_Jaccard","algorithm_valid_attempts","NMI_Study"], ["K"], rtol=1e-10, atol=1e-12)
        if label == "PAM_ALD":
            v = read(new/"original_K2.tsv").set_index("sample_id").cluster
            w = read(original/"original_K2.tsv").set_index("sample_id").cluster.reindex(v.index)
            ari = adjusted_rand_score(v,w)
            check("ALD_partition", ari == 1., adjusted_rand_index=ari,
                sizes=sorted(v.value_counts().tolist()), working_units=len(v))
    overlap = read(next((root/"overlaps").glob("*/overlap_counts.tsv")))
    z = overlap[overlap.scope.eq("M1_all")&overlap.body.eq("gut")&overlap.panel.eq("gut__full_source_reference")]
    expected = dict(zip(["02_AUD_or_ALD_vs_Health","03_AUD_vs_Health","04_ALD_vs_Health","05_ALD_vs_non_ALD"],[28,10,15,13]))
    for comparison, n in expected.items():
        observed = int(z.set_index("comparison").loc[comparison,"same_direction"])
        check("same_direction_overlap_"+comparison, observed == n, expected=n, observed=observed)
    save(pd.DataFrame(checks), root/"numerical_checks.tsv")
    dump(dict(status="PASS", checks=len(checks), clinical_coefficient_rows=len(co),
        clinical_test_positions=len(tests)), root/"numerical_receipt.json")
    print("NUMERICAL_VALIDATION_PASS", len(checks), "checks", len(co), "clinical coefficient rows")

if __name__ == "__main__":
    main()
