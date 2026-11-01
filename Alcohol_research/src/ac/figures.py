"""Principal English paper plots from frozen real data and full-precision results."""
from .inputs import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

def make(args, out):
    matplotlib.rcParams.update({"font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
        "svg.fonttype": "none"})
    data = Inputs(args.input_root)
    c = Path(args.clinical_root) if args.clinical_root else None
    index = []
    def export(fig, name, frame, source):
        fig.tight_layout()
        fig.savefig(out/(name+".svg"), bbox_inches="tight")
        plt.close(fig)
        save(frame, out/(name+".data.tsv.gz"))
        index.append(dict(figure=name, svg=name+".svg", data=name+".data.tsv.gz",
            source=source, source_sha256=sha(source)))
    source = data.root/"01_Alcohol_vs_Health/gut/M1_all/taxa"/args.panel/"all_results.tsv.gz"
    result = read(source)
    result = result[truth(result.selected)&result.BH_q.lt(.2)].sort_values(["beta", "taxon_key"])
    fig, ax = plt.subplots(figsize=(7, 10))
    x = np.arange(len(result))
    colors = np.where(result.beta.gt(0), "#D55E00", "#0072B2")
    for i, (_, row) in enumerate(result.iterrows()):
        ax.plot([row.CI_low, row.CI_high], [i, i], color=colors[i])
        ax.scatter(row.beta, i, color=colors[i], s=18)
    ax.set_yticks(x, result.genus_name)
    ax.axvline(0, color="0.6", lw=.7)
    ax.set(xlabel="Alcohol coefficient (CLR)", ylabel="Genus")
    export(fig, "main_candidates", result, source)
    source = data.i/"metadata_4603.tsv.gz"
    m = data.metadata()
    m["body"] = m.analysis_set.str.split("__").str[0]
    counts = m.groupby(["body", "Study", "group_broad"]).size().rename("specimens").reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={"width_ratios": [2, 1]})
    for ax, body in zip(axes, ["gut", "oral"]):
        z = counts[counts.body.eq(body)].pivot(index="Study", columns="group_broad", values="specimens").fillna(0)
        z.reindex(columns=["Health", "Alcohol"]).plot.barh(stacked=True, ax=ax,
            color=["#999999", "#D55E00"], width=.75)
        ax.set(xlabel="Specimens", ylabel="Study", title=body.capitalize())
        ax.legend(title="Frozen group")
    export(fig, "source_composition", counts, source)
    source = data.root/"beta_completed_summary.tsv"
    if source.exists():
        b = read(source)
        b = b[b.comparison.eq("01_Alcohol_vs_Health")&b.scope.eq("M1_all")]
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.barplot(data=b, x="panel_id", y="R2", hue="metric", ax=ax,
            palette=["#0072B2", "#E69F00"])
        ax.set(xlabel="Panel", ylabel="Marginal R squared")
        ax.tick_params(axis="x", labelrotation=65, labelsize=8)
        export(fig, "beta_effects", b, source)
    source = data.root/"01_Alcohol_vs_Health/gut/M1_all/functions/all_results.tsv.gz"
    if source.exists():
        f = read(source)
        f = f[truth(f.selected)&f.BH_q.lt(.2)].sort_values("beta")
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.errorbar(f.beta, np.arange(len(f)), xerr=[f.beta-f.CI_low, f.CI_high-f.beta],
            fmt="o", color="#0072B2", capsize=2)
        ax.set_yticks(np.arange(len(f)), f.taxon_key.str.replace("_", " "))
        ax.axvline(0, color="0.6", lw=.7)
        ax.set(xlabel="Alcohol coefficient (contribution fraction)", ylabel="Annotated trait")
        export(fig, "functional_candidates", f, source)
    if c:
        source = c/"interaction_models/all_interaction_tests.tsv.gz"
        tests = read(source)
        d = tests[tests.scope.eq("frozen")&tests.body.eq("gut")&tests.outcome.eq("taxa")
            &tests.panel.eq(args.panel)&tests.coverage.eq("base")&tests.q_interaction.lt(.2)]
        cofile = c/"interaction_models/all_model_coefficients.tsv.gz"
        co = read(cofile)
        co = co[co.model_id.eq("interaction")&co.term.str.startswith("Alcohol:")]
        keys = ["variable", "body", "scope", "outcome", "panel", "coverage", "feature_id"]
        f = co.merge(d[keys+["q_interaction", "q_interaction_joint_variables"]], on=keys)
        fig, ax = plt.subplots(figsize=(7, 4))
        f = f.sort_values(["variable", "feature_id", "term"])
        ax.errorbar(f.beta, np.arange(len(f)), xerr=[f.beta-f.CI_low, f.CI_high-f.beta],
            fmt="o", color="#0072B2", capsize=2)
        ax.set_yticks(np.arange(len(f)), f.variable+" x "+f.feature_id.str.split("::").str[-1])
        ax.axvline(0, color="0.6", lw=.7)
        ax.set(xlabel="Interaction coefficient", ylabel="Variable and genus")
        export(fig, "primary_interactions", f, cofile)
    source = data.root/"01_Alcohol_vs_Health/gut/M1_all/diversity"/args.panel/"alpha_M1.tsv"
    if source.exists():
        alpha = read(source)
        fig, axes = plt.subplots(1, 3, figsize=(9, 3))
        for ax, (_, row) in zip(axes, alpha.iterrows()):
            ax.errorbar([row.beta], [0], xerr=[[row.beta-row.CI_low], [row.CI_high-row.beta]],
                fmt="o", color="#0072B2", capsize=3)
            ax.axvline(0, color="0.6", lw=.7)
            ax.set(xlabel="Alcohol coefficient", title=row.metric, yticks=[])
        export(fig, "alpha_M1", alpha, source)
    disease = []
    for comparison in COMPARISONS[1:5]:
        source = data.root/comparison/"gut/M1_all/taxa"/args.panel/"all_results.tsv.gz"
        if source.exists():
            table = read(source)
            disease.append(dict(comparison=comparison, candidates=int(
                (truth(table.selected)&table.BH_q.lt(.2)).sum()), source_sha256=sha(source)))
    if disease:
        d = pd.DataFrame(disease)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.bar(["AUD or ALD", "AUD", "ALD", "ALD vs non-ALD"][:len(d)], d.candidates, color="#0072B2")
        ax.set(xlabel="Comparison", ylabel="Candidates")
        ax.tick_params(axis="x", labelsize=9)
        export(fig, "disease_candidates", d, source)
    source = data.root/"04_ALD_vs_Health/gut/M1_all/patterns/gut__common_all/PAM_completion/case_PCA_scores.tsv.gz"
    if source.exists():
        pca = read(source, index_col=0)
        assignment = read(source.parent/"original_K2.tsv").set_index("sample_id")
        pca = pca.join(assignment)
        fig, ax = plt.subplots(figsize=(5, 4))
        sns.scatterplot(data=pca, x="PC1", y="PC2", hue="cluster", palette=["#0072B2", "#D55E00"],
            s=18, ax=ax, linewidth=0)
        ax.set(xlabel="PC1", ylabel="PC2")
        ax.legend(title="Exploratory partition")
        export(fig, "ALD_case_patterns", pca.reset_index(), source)
    save(pd.DataFrame(index), out/"figure_index.tsv")
    return dict(figures=len(index), plotted_values="real frozen results; statistics not refitted by plotting",
        upstream_validation="REUSED_HASH_CHECKED; SVG rendering RECOMPUTED")
