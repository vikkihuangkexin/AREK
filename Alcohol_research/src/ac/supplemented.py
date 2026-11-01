"""Original supplemented-complete beta and case-pattern adapters."""
from .inputs import *
from .stats import nuisance, independent
from .beta_kernel import permutations, kernel
from scipy.spatial.distance import pdist, squareform
import subprocess

SCOPE = "M1_clinical_complete_evidence_supplemented"

def beta(args, out):
    from .beta import fast
    data = Inputs(args.input_root)
    t, y, _ = data.panel(args.panel)
    m = representatives(data.rows(args.comparison, y.index, "supplemented"), "20261030::")
    a = m.analysis_group.eq("Alcohol").to_numpy(float)
    N, _ = nuisance(m, "M1")
    N, _ = independent(N)
    reason = ("single_arm_or_empty" if m.analysis_group.nunique() < 2 else
        "fewer_than_5_units_in_arm_or_df_below_10" if min(m.analysis_group.value_counts()) < 5 or
        len(m)-N.shape[1]-1 < 10 else
        "exposure_in_required_nuisance_span" if np.linalg.norm(a-N@(np.linalg.pinv(N)@a)) < 1e-9 else "none")
    seed = int(hashlib.sha256((args.comparison+"::"+args.panel+"::evidence_clinical_Beta_v1").encode()).hexdigest()[:8], 16)
    if reason == "none":
        p, blocks, _ = permutations(m, 9999, seed)
        if not any(np.unique(a[b]).size == 2 for b in blocks):
            reason = "no_exchangeable_two_arm_block"
    rows = [dict(comparison=args.comparison, panel_id=args.panel, scope=SCOPE, metric=metric,
        n=len(m), n_case=int(a.sum()), n_reference=len(m)-int(a.sum()), reason=reason,
        F=np.nan, R2=np.nan, P=np.nan, dispersion_P=np.nan, B_actual=0)
        for metric in ("BrayCurtis", "Aitchison")]
    if reason == "none":
        save(t.loc[m.index], out/"TSS.tsv.gz", True)
        save(y.loc[m.index], out/"CLR.tsv.gz", True)
        save(m, out/"metadata.tsv", True)
        save(pd.DataFrame(N, index=m.index), out/"nuisance.tsv", True)
        save(pd.DataFrame(p+1), out/"permutations_1based.tsv.gz")
        for row in rows:
            metric = row["metric"]
            z = t.loc[m.index] if metric == "BrayCurtis" else y.loc[m.index]
            G = kernel(squareform(pdist(z, "braycurtis" if metric == "BrayCurtis" else "euclidean")))
            result, ff = fast(G, N, a, np.argsort(p, axis=1), blocks)
            row.update(result)
            save(pd.DataFrame(dict(F=ff)), out/(metric+"_location_Fperm.tsv.gz"))
        subprocess.run(["Rscript", str(Path(__file__).with_name("beta_completion.R")), str(out)], check=True)
        r = read(out/"R_results.tsv").set_index("metric")
        for row in rows:
            metric = row["metric"]
            ff = read(out/(metric+"_location_Fperm.tsv.gz")).F.to_numpy()
            rf = read(out/(metric+"_vegan_location_Fperm99.tsv.gz")).F.to_numpy()
            np.testing.assert_allclose(ff[:99], rf, atol=1e-7, rtol=1e-7)
            row["F"] = float(r.loc[metric, "vegan_location_F"])
            row["R2"] = float(r.loc[metric, "vegan_location_R2"])
            row["P"] = float(((ff >= row["F"]-np.sqrt(np.finfo(float).eps)).sum()+1)/10000)
            row.update(r.loc[metric].to_dict())
            row["seed"] = seed
            row["permutation_hash"] = hashlib.sha256(p.tobytes()).hexdigest()
        (out/"permutations_1based.tsv.gz").unlink()
    save(pd.DataFrame(rows), out/"completed_tests.tsv")
    reference = data.root/"evidence_supplemented/analyses"/args.comparison/args.panel.split("__")[0]/SCOPE/"diversity"/args.panel/"completion/completed_tests.tsv"
    if reference.exists():
        expected = read(reference).set_index("metric").reindex([r["metric"] for r in rows])
        for key in ["F", "R2", "P", "dispersion_P"]:
            np.testing.assert_allclose([r[key] for r in rows], expected[key], equal_nan=True, rtol=1e-7, atol=1e-7)
    return dict(status="PASS", reason=reason, seed=seed, results=rows, model="subset_M1_not_clinical_adjusted")

def patterns(args, out):
    from . import patterns as implementation
    data = Inputs(args.input_root)
    implementation.rows = lambda comparison, ids, scope, cases_only=False: data.rows(
        comparison, ids, "supplemented", cases_only)
    return implementation.run_node(args.comparison, SCOPE, args.panel, {})
