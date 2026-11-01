"""Input auditing, evidence replay, frozen harmonization, and memberships."""
from .inputs import *
from .stats import clr
import shutil
import zipfile
import re

def audit(args, out):
    data = Inputs(args.input_root)
    m, membership = data.metadata(), data.membership()
    assert m.index.equals(membership.index) and m.index.is_unique
    records = []
    for folder in sorted((data.i/"panels").iterdir()):
        t, y, family = data.panel(folder.name)
        assert t.index.equals(y.index) and t.columns.equals(y.columns)
        assert set(t.index) <= set(m.index)
        assert not family.taxon_key.duplicated().any()
        assert set(family.taxon_key) == set(t.columns)
        assert np.array_equal(t.isna(), y.isna())
        records.append(dict(panel=folder.name, specimens=len(t), profile_features=t.shape[1],
            selected=int(truth(family.selected).sum()), status="PASS"))
    save(pd.DataFrame(records), out/"schema_checks.tsv")
    save(pd.DataFrame(file_manifest(data.i)), out/"input_hashes.tsv")
    archives = []
    for path in getattr(args, "archives", []) or []:
        p = Path(path)
        with zipfile.ZipFile(p) as z:
            bad = z.testzip()
            if bad:
                raise ValueError("Archive CRC failure: "+bad)
        identity = sha(p)
        pins = json.loads((Path(__file__).parents[2]/"configs/frozen_definitions.json").read_text())["archive_hashes"]
        if p.name in pins and identity != pins[p.name]:
            raise ValueError("Authoritative archive identity mismatch: "+p.name)
        archives.append(dict(filename=p.name, bytes=p.stat().st_size, sha256=identity))
    dump(archives, out/"archives.json")
    available = {tool: shutil.which(tool) is not None for tool in
        ["qiime", "fastp", "bowtie2", "samtools", "kraken2", "bracken", "Rscript"]}
    dump(dict(profile_stage="AVAILABLE", raw_reads="NOT_VERIFIED_AVAILABLE",
        raw_database_assets="NOT_VERIFIED_AVAILABLE", executables=available,
        raw_execution="BLOCKED_without_original_reads_and_database_assets"), out/"stage_availability.json")
    return dict(specimens=len(m), Studies=m.Study.nunique(), panels=len(records), raw_FASTQ_validation="NOT_RUN")

def harmonize_metadata(args, out):
    """Replay accepted exact-ID evidence and unresolved conflict decisions."""
    c = Path(args.clinical_root)
    evidence = read(c/"evidence/individual_value_evidence.tsv.gz", dtype=str, keep_default_na=False)
    expected = read(c/"harmonised_clinical_values.tsv.gz", index_col=0, low_memory=False)
    m = read(c/"reproducibility/inputs/metadata.tsv.gz", index_col=0)
    assert expected.index.equals(m.index) and m.index.is_unique
    values = pd.DataFrame(index=m.index)
    registry = json.loads((Path(__file__).parents[2]/"configs/clinical_registry.json").read_text())
    continuous = {x["variable"] for x in registry if x["kind"] == "continuous"}
    continuous |= {"bilirubin", "INR", "prothrombin_time", "Child_Pugh", "MELD",
        "heart_rate", "weight", "height", "waist", "hip", "BSS"}
    conflicts = []
    for (variable, sid), group in evidence.groupby(["variable", "sample_id"], sort=False):
        if sid not in m.index:
            raise ValueError("Evidence ID is outside the frozen specimen population: "+sid)
        accepted = group.status.str.startswith("accepted")
        if not accepted.all():
            conflicts.append(dict(variable=variable, sample_id=sid,
                reason="unresolved_source_conflict_set_missing"))
            continue
        raw = group.value
        equal = np.allclose(pd.to_numeric(raw), float(raw.iloc[0]), rtol=1e-6, atol=1e-6) if variable in continuous else raw.nunique() == 1
        if not equal:
            raise ValueError("Accepted evidence contains unresolved conflicting values.")
        if variable not in values:
            values[variable] = pd.Series(index=m.index, dtype=object)
        values.at[sid, variable] = raw.iloc[0]
    checks = []
    for variable in values:
        actual = values[variable]
        frozen = expected[variable]
        assert actual.notna().equals(frozen.notna()), variable
        if variable in continuous or pd.api.types.is_numeric_dtype(frozen):
            values[variable] = pd.to_numeric(actual)
            np.testing.assert_allclose(values[variable].to_numpy(float), pd.to_numeric(frozen).to_numpy(float),
                equal_nan=True, rtol=1e-12, atol=1e-12)
        else:
            assert actual.fillna("").astype(str).equals(frozen.fillna("").astype(str)), variable
        checks.append(dict(variable=variable, nonmissing=int(actual.notna().sum()), status="PASS"))
    save(values, out/"normalized_clinical_values.tsv.gz", True)
    save(m, out/"normalized_metadata.tsv.gz", True)
    save(evidence, out/"exact_ID_provenance.tsv.gz")
    save(pd.DataFrame(conflicts), out/"conflicts.tsv")
    save(pd.DataFrame(checks), out/"numerical_comparison.tsv")
    return dict(specimens=len(m), reconstructed_evidence_variables=len(values.columns),
        conflicts=len(conflicts), review_only_fields="Remain in external frozen inventory; not promoted into models.",
        technical_metadata="REUSED_HASH_CHECKED; clinical values recomputed from accepted evidence")

def profiles(args, out):
    data = Inputs(args.input_root)
    t, y, family = data.panel(args.panel)
    source = getattr(args, "source_root", None)
    canonical_mode = "REUSED_HASH_CHECKED"
    if source:
        root = Path(source)
        crosswalk = read(root/"provenance/taxon_crosswalk.tsv.gz", keep_default_na=False)
        def canonical(row):
            parts = dict(re.findall(r"(?:^|;)([dpcofg])__([^;]+)", row.taxonomy_lineage))
            anchor = parts.get("f", "")
            confident = row.resolution_status == "genus" and anchor and not re.search(
                "unclassified|uncultured|incertae|Family_", anchor, re.I)
            return "FG::"+parts.get("d", "")+"::"+anchor+"::"+row.genus_name if confident else (
                "SRC::"+row.analysis_set+"::"+row.original_feature_id)
        assert crosswalk.apply(canonical, axis=1).equals(crosswalk.taxon_key)
        metadata = data.metadata().loc[t.index]
        ts, ys = [], []
        for name in sorted(metadata.analysis_set.unique()):
            folder = root/"provenance/source_inputs"/name
            st = read(folder/"TSS.tsv.gz", index_col=0)
            sy = read(folder/"CLR_unadjusted.tsv.gz", index_col=0)
            mapping = crosswalk[crosswalk.analysis_set.eq(name)&
                ~truth(crosswalk.excluded_from_genus_inference)].set_index("original_feature_id")
            labels = mapping.loc[st.columns, "taxon_key"]
            if labels.duplicated().any():
                raise ValueError("Ambiguous within-source canonical mapping.")
            st.columns = labels.to_list()
            sy.columns = labels.to_list()
            ts.append(st)
            ys.append(sy)
        reconstructed_t = pd.concat(ts, axis=0, sort=True).reindex(index=t.index, columns=t.columns)
        if "full_source" in args.panel:
            reconstructed = pd.concat(ys, axis=0, sort=True).reindex(index=y.index, columns=y.columns)
        else:
            reconstructed_t = reconstructed_t.div(reconstructed_t.sum(axis=1), axis=0)
            reconstructed = clr(reconstructed_t)
        np.testing.assert_allclose(reconstructed_t, t, equal_nan=True, rtol=1e-10, atol=1e-12)
        t = reconstructed_t
        save(crosswalk, out/"source_canonical_crosswalk.tsv.gz")
        canonical_mode = "RECOMPUTED"
    else:
        reconstructed = clr(t)
    assert np.array_equal(reconstructed.isna(), y.isna())
    np.testing.assert_allclose(reconstructed, y, equal_nan=True, rtol=1e-10, atol=1e-10)
    assert family.taxon_key.is_unique
    save(t, out/"TSS.tsv.gz", True)
    save(reconstructed, out/"CLR_unadjusted.tsv.gz", True)
    save(family, out/"family.tsv")
    save(t.notna().astype("uint8"), out/"measurability_mask.tsv.gz", True)
    save(family[["taxon_key", "genus_name", "taxon_status"]], out/"canonical_taxa.tsv")
    return dict(panel=args.panel, profile_features=len(family), selected=int(truth(family.selected).sum()),
        structural_missing=int(t.isna().sum().sum()), selection="frozen_not_reselected",
        canonical_mapping=canonical_mode,
        earliest_stage="source-specific frozen TSS/CLR plus crosswalk" if source else "canonical specimen-by-feature TSS")

def comparisons(args, out):
    data = Inputs(args.input_root)
    m, labels = data.metadata(), data.membership()
    assert m.index.equals(labels.index)
    body = m.analysis_set.str.split("__").str[0]
    records = []
    for comparison in COMPARISONS:
        for site in sorted(body.unique()):
            v = labels.loc[body.eq(site), comparison]
            records.append(dict(comparison=comparison, body=site,
                case=int(v.eq("case").sum()), reference=int(v.eq("reference").sum()),
                excluded=int(v.eq("excluded").sum())))
    save(labels, out/"comparison_membership.tsv.gz", True)
    save(pd.DataFrame(records), out/"population_counts.tsv")
    save(labels[truth(labels.clinical_complete_age_sex_BMI)], out/"original_clinical_complete.tsv.gz", True)
    if args.clinical_root:
        c = Path(args.clinical_root)
        v = read(c/"harmonised_clinical_values.tsv.gz", index_col=0, low_memory=False)
        cm = read(c/"reproducibility/inputs/metadata.tsv.gz", index_col=0)
        rows = []
        for var in ("age", "sex", "BMI"):
            for site in ("gut", "oral"):
                ids = cm.index[cm.body.eq(site)&v[var].notna()]
                for sid in ids:
                    rows.append(dict(variable=var, body=site, sample_id=sid,
                        group=cm.at[sid, "analysis_group"]))
        save(pd.DataFrame(rows), out/"individual_variable_membership.tsv.gz")
        save(cm[truth(cm.explicit_disease_reference)], out/"explicit_disease_references.tsv.gz", True)
    return dict(populations=records, definitions="Frozen labels preserved, including overlapping AUD and ALD.")
