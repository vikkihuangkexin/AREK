"""Frozen profile and membership readers; no label inference or reselection."""
from .io import *
import importlib.metadata

COMPARISONS = [
    "01_Alcohol_vs_Health", "02_AUD_or_ALD_vs_Health", "03_AUD_vs_Health",
    "04_ALD_vs_Health", "05_ALD_vs_non_ALD", "06_AUD_vs_non_AUD",
    "07_Cirrhosis_vs_non_cirrhotic_ALD", "08_ALD_vs_non_alcoholic_liver_disease"]

class Inputs:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.i = self.root/"reproducibility/inputs"
        if not self.i.is_dir():
            raise FileNotFoundError("Use the extracted Alcohol_Independent_Comparisons_M1 root or an external original verification-input root.")

    def metadata(self):
        m = read(self.i/"metadata_4603.tsv.gz", index_col=0, keep_default_na=False)
        require_columns(m, ["Study", "analysis_set", "assay_analysis", "material_analysis",
            "technical_batch", "working_cluster", "force_cluster_source", "group_broad"])
        for column in ("age_model", "BMI_model"):
            m[column] = pd.to_numeric(m[column], errors="coerce")
        m["sex_model"] = m.sex_model.replace({"": np.nan, "NA": np.nan, "unknown": np.nan})
        m["force_cluster_source"] = truth(m.force_cluster_source)
        m["sample_id"] = m.index
        return m

    def membership(self):
        return read(self.i/"comparison_membership_4603.tsv.gz", index_col=0)

    def panel(self, name):
        folder = self.i/"panels"/name
        return read(folder/"TSS.tsv.gz", index_col=0), read(folder/"CLR_unadjusted.tsv.gz", index_col=0), read(folder/"family.tsv")

    def rows(self, comparison, ids, scope="M1_all", cases_only=False):
        m = self.metadata().loc[ids].copy()
        labels = self.membership().loc[ids]
        m["analysis_group"] = labels[comparison].map({"case": "Alcohol", "reference": "Health"})
        m = m[m.analysis_group.eq("Alcohol") if cases_only else m.analysis_group.notna()]
        if scope == "M1_clinical_complete":
            m = m[truth(labels.loc[m.index, "clinical_complete_age_sex_BMI"])]
        elif scope == "supplemented":
            evidence = read(self.root/"evidence_supplemented/metadata_evidence_supplemented.tsv.gz", index_col=0)
            m = m[evidence.loc[m.index, ["age_model", "sex_model", "BMI_model"]].notna().all(axis=1)]
        elif scope != "M1_all":
            raise ValueError("Unsupported frozen M1 scope: " + scope)
        return m

def representatives(m, salt):
    m = m.copy()
    m["_selection"] = [hashlib.sha256((salt+s).encode()).hexdigest() for s in m.index]
    return m.sort_values("_selection").drop_duplicates("working_cluster").sort_index().drop(columns="_selection")

def versions():
    return {p: importlib.metadata.version(p)
        for p in ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]}
