"""Completed group M1 analyses with frozen per-feature support and families."""
from .inputs import *
from .stats import bh
from .m1 import fit_matrix, apply_families

def model_inputs(args, outcome="taxa"):
    data = Inputs(args.input_root)
    body = args.panel.split("__")[0]
    folder = data.i/("panels" if outcome == "taxa" else "functions" if outcome == "function" else "diversity")
    folder = folder/(body if outcome == "function" else args.panel)
    matrix = folder/("CLR_unadjusted.tsv.gz" if outcome == "taxa" else
        "contribution.tsv.gz" if outcome == "function" else "alpha_per_specimen.tsv.gz")
    files = [data.i/"metadata_4603.tsv.gz", data.i/"comparison_membership_4603.tsv.gz", matrix]
    if outcome != "alpha":
        files.append(folder/"family.tsv")
    if args.scope == "supplemented":
        files.append(data.root/"evidence_supplemented/metadata_evidence_supplemented.tsv.gz")
    if getattr(args, "coverage", "base") != "base":
        files.append(folder/("global_coverage.tsv.gz" if args.coverage == "global_coverage" else "known_mass.tsv.gz"))
    return files

def run_m1(args, out, outcome="taxa"):
    data = Inputs(args.input_root)
    files = model_inputs(args, outcome)
    y = read(files[2], index_col=0)
    if outcome == "alpha":
        y = y[["Shannon", "Simpson", "Observed_genera"]]
        family = pd.DataFrame(dict(taxon_key=y.columns, selected=True))
    else:
        family = read(files[3])
        # The original main and original clinical-complete tables retain the full
        # reference coordinates. Disease and supplemented tables use selected rows.
        if outcome == "taxa" and not (args.comparison == "01_Alcohol_vs_Health" and
                args.scope in {"M1_all", "M1_clinical_complete"}):
            family = family[truth(family.selected)]
        y = y.reindex(columns=family.taxon_key)
    m = data.rows(args.comparison, y.index, args.scope)
    extra = None
    coverage = getattr(args, "coverage", "base")
    if coverage == "global_coverage":
        extra = read(files[-1], index_col=0)[["original_abundance_coverage"]].reindex(m.index)
        m = m[extra.notna().all(axis=1)]
        extra = extra.loc[m.index]
    if coverage == "trait_coverage":
        known = read(files[-1], index_col=0)
        results = []
        for key in y:
            ids = m.index[known.reindex(m.index)[key].notna()]
            results.append(fit_matrix(m.loc[ids], y.loc[ids, [key]], "M1",
                extra=known.loc[ids, [key]].rename(columns={key: "trait_known_mass"})))
        result = pd.concat(results, ignore_index=True)
    else:
        result = fit_matrix(m, y.loc[m.index], "M1", extra=extra)
    result = apply_families(result, family)
    save(result, out/"all_results.tsv.gz")
    save(m, out/"membership.tsv.gz", True)
    save(result[result.BH_q.lt(.2)], out/"candidates.tsv")
    comparison = []
    if args.scope == "supplemented":
        reference = data.root/"evidence_supplemented/analyses"/args.comparison/args.panel.split("__")[0]/"M1_clinical_complete_evidence_supplemented"
    else:
        reference = data.root/args.comparison/args.panel.split("__")[0]/args.scope
    reference = reference/("taxa"/Path(args.panel) if outcome == "taxa" else Path("functions") if outcome == "function" else Path("diversity")/args.panel)
    expected_file = reference/("alpha_M1.tsv" if outcome == "alpha" else "all_results.tsv.gz")
    if coverage != "base":
        expected_file = reference/"coverage_sensitivity_results.tsv.gz"
    if expected_file.exists():
        expected = read(expected_file)
        if outcome == "alpha":
            expected = expected.rename(columns={"metric": "taxon_key"})
        if coverage != "base":
            labels = ["M1_coverage", "M1_global_coverage"] if coverage == "global_coverage" else ["M1_trait_coverage"]
            expected = expected[expected.model.isin(labels)]
        expected = expected.set_index("taxon_key").reindex(result.taxon_key)
        for column in ("beta", "SE", "CI_low", "CI_high", "raw_p", "BH_q", "full_reference_BH_q"):
            if column not in expected:
                continue
            actual = result[column].to_numpy(float)
            frozen = expected[column].to_numpy(float)
            rtol, atol = (1e-10, 1e-12) if column.endswith("BH_q") else (2e-7, 2e-9)
            np.testing.assert_allclose(actual, frozen, rtol=rtol, atol=atol, equal_nan=True,
                err_msg=str(expected_file)+":"+column)
            error = abs(actual-frozen)
            comparison.append(dict(column=column, status="PASS", rtol=rtol, atol=atol,
                max_absolute_error=float(np.nanmax(error)) if np.isfinite(error).any() else None))
        save(pd.DataFrame(comparison), out/"numerical_comparison.tsv")
    return dict(outcome=outcome, comparison=args.comparison, scope=args.scope,
        panel=args.panel, coverage=coverage, specimens=m.analysis_group.value_counts().to_dict(),
        fixed_family=int(result.selected.sum()), full_reference_features=len(result),
        estimable=int(result.loc[result.selected, "raw_p"].notna().sum()),
        candidates=int(result.BH_q.lt(.2).sum()),
        positive=int((result.BH_q.lt(.2)&result.beta.gt(0)).sum()),
        negative=int((result.BH_q.lt(.2)&result.beta.lt(0)).sum()),
        numerical_reference=str(expected_file) if expected_file.exists() else "NOT_AVAILABLE")

def alpha(args, out):
    """Recompute Shannon, Gini-Simpson and observed genera before M1."""
    data = Inputs(args.input_root)
    t, _, _ = data.panel(args.panel)
    from scipy.stats import entropy
    def measures(table):
        p = table.div(table.sum(axis=1), axis=0).to_numpy()
        return pd.DataFrame(dict(Shannon=entropy(p, axis=1), Simpson=1-(p*p).sum(axis=1),
            Observed_genera=(p>0).sum(axis=1)), index=table.index)
    stored = read(data.i/"diversity"/args.panel/"alpha_per_specimen.tsv.gz", index_col=0)
    mode = "RECOMPUTED"
    if "full_source" in args.panel:
        source = getattr(args, "source_root", None)
        if source:
            tables = sorted((Path(source)/"provenance/source_inputs").glob("*/TSS.tsv.gz"))
            if not tables:
                raise FileNotFoundError("BLOCKED: original source-specific TSS compositions are absent.")
            values = pd.concat([measures(read(p, index_col=0)) for p in tables]).loc[t.index]
        else:
            values = stored[["Shannon", "Simpson", "Observed_genera"]].copy()
            mode = "REUSED_HASH_CHECKED"
    else:
        values = measures(t)
    np.testing.assert_allclose(values, stored.loc[values.index, values.columns],
        atol=1e-10, rtol=1e-10, equal_nan=True)
    save(values, out/"alpha_per_specimen.tsv.gz", True)
    result = run_m1(args, out, "alpha")
    result["alpha_measures"] = mode
    result["source_definition"] = "Source-specific full compositions, not the canonical pooled genus projection."
    return result
