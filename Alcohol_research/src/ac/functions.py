"""Accepted FAPROTAX mapping and original-scale taxon-to-trait contributions."""
from .inputs import *
import subprocess

def annotation_inputs(args):
    v = Path(args.source_root)
    base = v/"functional/faprotax"
    files = [base/"genus_reference_audit.tsv.gz", base/"official_mapper_input.tsv",
        base/"official/genus_function.tsv", base/"official/placements.tsv"]
    if args.rerun_mapper:
        db = Path(args.database_root)
        files += [db/"faprotax.R"] + sorted((db/"FAPROTAX_ref_data").glob("*"))
    for body in ("gut", "oral"):
        for layer in sorted((v/"inputs_v1"/body/"layers").iterdir()):
            if layer.is_dir():
                files.append(layer/"TSS.tsv.gz")
    return [x for x in files if not x.is_dir()]

def annotation(args, out):
    v = Path(args.source_root)
    base = v/"functional/faprotax"
    audit = read(base/"genus_reference_audit.tsv.gz", keep_default_na=False)
    query = read(base/"official_mapper_input.tsv", keep_default_na=False)
    source = base/"official"
    if args.rerun_mapper:
        db = Path(args.database_root)
        source = out/"official"
        source.mkdir()
        save(query, source/"input.tsv")
        command = ["Rscript", str(db/"faprotax.R"), "-i", str(source/"input.tsv"), "-k", "OTU",
            "-t", "column:taxonomy", "-d", str(db/"FAPROTAX_ref_data"), "-n", "none",
            "-o", str(source/"function_table.tsv"), "-e", str(source/"genus_function.tsv"),
            "-q", str(source/"probabilities.tsv"), "--out_nearest_refs", str(source/"placements.tsv"),
            "-j", str(source/"intermediates"), "--placement_failure", "omit_otus", "--Nthreads", "1",
            "-l", str(source/"log.txt"), "-r", str(source/"report.txt")]
        subprocess.run(command, check=True)
        dump(command, out/"mapper_command.json")
    placements = read(source/"placements.tsv", comment="#", keep_default_na=False)
    placements = query.merge(placements, on="OTU", validate="one_to_one")
    placements["genus_placement_pass"] = placements.taxonomy.eq(placements.parent_node_taxonomy)
    states = read(source/"genus_function.tsv", comment="#", index_col=0)
    frozen = read(base/"accepted_genus_function_1_0_NA.tsv.gz", index_col=0)
    accepted = placements.loc[placements.genus_placement_pass, "OTU"]
    np.testing.assert_allclose(states.loc[accepted].reindex(index=frozen.index, columns=frozen.columns),
        frozen, equal_nan=True, rtol=0, atol=0)
    save(states.loc[accepted], out/"accepted_genus_function.tsv.gz", True)
    save(placements, out/"placement_audit.tsv")
    audit = audit.merge(placements[["OTU", "genus_placement_pass"]],
        left_on="annotation_id", right_on="OTU", how="left", validate="many_to_one")
    audit["genus_placement_pass"] = audit.genus_placement_pass.eq(True)
    save(audit, out/"mapping_and_unmatched.tsv.gz")
    summaries = []
    for body in ("gut", "oral"):
        for layer in sorted((v/"inputs_v1"/body/"layers").iterdir()):
            if not layer.is_dir():
                continue
            t = read(layer/"TSS.tsv.gz", index_col=0)
            a = audit[audit.analysis_set.eq(layer.name)].set_index("feature_id").reindex(t.columns)
            z = pd.DataFrame(np.nan, index=t.columns, columns=states.columns)
            good = a.genus_placement_pass.eq(True)
            z.loc[good] = states.reindex(a.loc[good, "annotation_id"]).to_numpy()
            mass = t.to_numpy()@z.notna().to_numpy(float)
            values = t.to_numpy()@z.eq(1).to_numpy(float)
            values[mass <= 0] = np.nan
            folder = out/body/layer.name
            save(pd.DataFrame(values, index=t.index, columns=z.columns), folder/"contribution.tsv.gz", True)
            save(pd.DataFrame(mass, index=t.index, columns=z.columns), folder/"known_mass.tsv.gz", True)
            coverage = pd.DataFrame(dict(original_abundance_coverage=(t*z.notna().any(axis=1)).sum(axis=1)))
            save(coverage, folder/"global_coverage.tsv.gz", True)
            expected = v/body/"faprotax"/layer.name/"annotated_abundance_contribution.tsv.gz"
            np.testing.assert_allclose(values, read(expected, index_col=0).reindex(index=t.index, columns=z.columns),
                rtol=1e-10, atol=1e-12, equal_nan=True)
            summaries.append(dict(analysis_set=layer.name, specimens=len(t),
                accepted_genus_placements=int(good.sum()), comparison="PASS"))
    save(pd.DataFrame(summaries), out/"source_checks.tsv")
    return dict(layers=len(summaries), accepted_genera=len(accepted), traits=states.shape[1],
        mapper="RECOMPUTED" if args.rerun_mapper else "REUSED_HASH_CHECKED",
        source_contributions="RECOMPUTED", original_TSS_scale=True, cross_trait_reclosure=False)

def contributions(args, out):
    data = Inputs(args.input_root)
    mapping = read(data.i/"candidate_taxon_trait_mapping.tsv.gz")
    pieces = []
    for aset, group in mapping.groupby("analysis_set", sort=True):
        t = read(data.i/"candidate_function_original_TSS"/(aset+".tsv.gz"), index_col=0)
        for (taxon, trait), g in group.groupby(["taxon_key", "trait"], sort=True):
            values = pd.Series(0., index=t.index)
            for row in g.itertuples():
                values = values+t[row.original_feature_id]*float(row.accepted_state)
            pieces.append(pd.DataFrame(dict(sample_id=values.index, analysis_set=aset,
                taxon_key=taxon, trait=trait, annotated_taxon_contribution=values.to_numpy())))
    long = pd.concat(pieces, ignore_index=True)
    expected = read(data.i/"candidate_taxon_trait_contribution_all_supported_samples.tsv.gz")
    keys = ["sample_id", "analysis_set", "taxon_key", "trait"]
    a = long.set_index(keys).annotated_taxon_contribution.sort_index()
    b = expected.set_index(keys).annotated_taxon_contribution.reindex(a.index)
    np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12, equal_nan=True)
    labels = data.membership()
    long["comparison_group"] = labels.loc[long.sample_id, args.comparison].to_numpy()
    long = long[long.comparison_group.isin(["case", "reference"])]
    if args.scope != "M1_all":
        m = data.rows(args.comparison, labels.index, args.scope)
        long = long[long.sample_id.isin(m.index)]
    save(long, out/"contributions_per_ID.tsv.gz")
    summary = long.groupby(["taxon_key", "trait", "comparison_group"]).annotated_taxon_contribution.agg(
        n_specimens="count", mean_contribution="mean", median_contribution="median").reset_index()
    save(summary, out/"group_summary.tsv")
    return dict(rows=len(long), groups=len(summary), contribution_scale="original_TSS_nonexclusive_unreclosed")
