"""Canonical direction-aware overlaps, retaining panel and scope identities."""
from .inputs import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def candidate_overlap(args, out):
    data = Inputs(args.input_root)
    records, counts = [], []
    primary_file = data.root/"01_Alcohol_vs_Health"/"gut/M1_all/taxa"/args.panel/"all_results.tsv.gz"
    primary = read(primary_file).set_index("taxon_key")
    primary = primary[truth(primary.selected)&primary.BH_q.lt(.2)]
    paths = sorted(data.root.glob("0*/*/M1*/taxa/*/all_results.tsv.gz"))
    paths += sorted((data.root/"evidence_supplemented/analyses").glob("0*/*/M1*/taxa/*/all_results.tsv.gz"))
    for path in paths:
        parts = path.relative_to(data.root).parts
        if parts[0] == "evidence_supplemented":
            parts = parts[2:]
        comparison, body, scope, _, panel = parts[:5]
        result = read(path)
        result = result[truth(result.selected)&result.BH_q.lt(.2)]
        common = result[result.taxon_key.isin(primary.index)].copy()
        common["primary_beta"] = primary.beta.reindex(common.taxon_key).to_numpy()
        common["same_direction"] = np.sign(common.beta).eq(np.sign(common.primary_beta))
        for r in common.to_dict("records"):
            records.append({**r, "comparison": comparison, "body": body, "scope": scope, "panel": panel,
                "primary_panel": args.panel, "primary_scope": "M1_all"})
        counts.append(dict(comparison=comparison, body=body, scope=scope, panel=panel,
            candidates=len(result), canonical_overlap=len(common), same_direction=int(common.same_direction.sum())))
    save(pd.DataFrame(records), out/"canonical_overlap.tsv.gz")
    summary = pd.DataFrame(counts)
    save(summary, out/"overlap_counts.tsv")
    if args.clinical_root:
        tests = read(Path(args.clinical_root)/"interaction_models/all_interaction_tests.tsv.gz")
        hits = tests[tests.q_interaction.lt(.2)].copy()
        hits["overlaps_primary_canonical_ID"] = hits.feature_id.isin(primary.index)
        save(hits, out/"interaction_overlap.tsv")
    matplotlib.rcParams.update({"font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"], "svg.fonttype": "none"})
    selected = summary[summary.body.eq("gut")&summary.scope.eq("M1_all")&summary.panel.eq(args.panel)
        &summary.comparison.isin(COMPARISONS[1:5])].copy()
    selected["label"] = ["AUD or ALD", "AUD", "ALD", "ALD vs non-ALD"][:len(selected)]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.bar(selected.label, selected.same_direction, color="#0279EE")
    ax.set(xlabel="Comparison", ylabel="Same-direction overlap")
    ax.tick_params(axis="x", labelsize=9)
    fig.tight_layout()
    fig.savefig(out/"disease_overlap.svg")
    plt.close(fig)
    save(selected, out/"disease_overlap.data.tsv")
    return dict(primary_candidates=len(primary), comparisons=len(summary),
        matching="canonical IDs, never genus-name strings", source_results="REUSED_HASH_CHECKED")
