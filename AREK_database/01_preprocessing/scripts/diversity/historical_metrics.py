# Extracted in-scope historical implementation; see docs/source_script_map.tsv.
from __future__ import annotations
from typing import Dict
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
def clr_transform(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    x = counts.astype(float) + pseudocount
    logx = np.log(x)
    gm = logx.mean(axis=0)
    clr = logx.sub(gm, axis=1)
    return clr

def alpha_metrics(counts: pd.DataFrame) -> pd.DataFrame:
    # counts taxa x samples
    rows = []
    for sample in counts.columns:
        x = counts[sample].astype(float).values
        total = x.sum()
        if total <= 0:
            observed = shannon = simpson = chao1 = np.nan
        else:
            pos = x[x > 0]
            p = pos / total
            observed = len(pos)
            shannon = -np.sum(p * np.log(p))
            simpson = 1.0 - np.sum(p ** 2)
            f1 = np.sum(x == 1)
            f2 = np.sum(x == 2)
            chao1 = observed + (f1 * (f1 - 1)) / (2 * (f2 + 1))
        rows.append({
            "analysis_sample_id": sample,
            "Observed": observed,
            "Shannon": shannon,
            "Simpson": simpson,
            "Chao1": chao1,
        })
    return pd.DataFrame(rows).set_index("analysis_sample_id")

def compare_groups_numeric(values: pd.Series, groups: pd.Series) -> Dict[str, float | str]:
    df = pd.DataFrame({"value": values, "group": groups}).dropna()
    levels = list(df["group"].astype(str).unique())
    if len(levels) < 2:
        return {"method": "not_tested", "statistic": np.nan, "p_value": np.nan, "effect_size": np.nan}
    if len(levels) == 2:
        a = df.loc[df["group"].astype(str) == levels[0], "value"].astype(float)
        b = df.loc[df["group"].astype(str) == levels[1], "value"].astype(float)
        if len(a) < 2 or len(b) < 2:
            return {"method": "not_enough_samples", "statistic": np.nan, "p_value": np.nan, "effect_size": np.nan}
        try:
            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except Exception:
            stat, p = np.nan, np.nan
        effect = float(np.nanmedian(a) - np.nanmedian(b))
        return {"method": "Mann-Whitney_U", "statistic": stat, "p_value": p, "effect_size": effect, "contrast": f"{levels[0]}_vs_{levels[1]}"}
    else:
        arrays = [df.loc[df["group"].astype(str) == lev, "value"].astype(float).values for lev in levels]
        arrays = [x for x in arrays if len(x) >= 2]
        if len(arrays) < 2:
            return {"method": "not_enough_samples", "statistic": np.nan, "p_value": np.nan, "effect_size": np.nan}
        try:
            stat, p = stats.kruskal(*arrays)
        except Exception:
            stat, p = np.nan, np.nan
        return {"method": "Kruskal-Wallis", "statistic": stat, "p_value": p, "effect_size": np.nan, "contrast": "overall"}

def bray_curtis_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    # samples x taxa
    X = counts.T.astype(float).values
    D = squareform(pdist(X, metric="braycurtis"))
    return pd.DataFrame(D, index=counts.columns, columns=counts.columns)

def jaccard_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    X = (counts.T.astype(float).values > 0).astype(int)
    D = squareform(pdist(X, metric="jaccard"))
    return pd.DataFrame(D, index=counts.columns, columns=counts.columns)

def aitchison_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    clr = clr_transform(counts).T.values
    D = squareform(pdist(clr, metric="euclidean"))
    return pd.DataFrame(D, index=counts.columns, columns=counts.columns)

def pcoa_from_distance(D: pd.DataFrame) -> pd.DataFrame:
    # Classical MDS / PCoA
    X = D.values.astype(float)
    n = X.shape[0]
    if n < 2:
        return pd.DataFrame(index=D.index, columns=["PCoA1", "PCoA2", "var1", "var2"])
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(X ** 2).dot(J)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
    total = np.sum(np.maximum(eigvals, 0))
    var = np.maximum(eigvals[:2], 0) / total if total > 0 else [np.nan, np.nan]
    return pd.DataFrame({"PCoA1": coords[:, 0], "PCoA2": coords[:, 1], "var1": var[0], "var2": var[1]}, index=D.index)

def permanova(D: pd.DataFrame, groups: pd.Series, n_perm: int = 999, seed: int = 1) -> Dict[str, float]:
    # One-way PERMANOVA pseudo-F on distance matrix.
    rng = np.random.default_rng(seed)
    common = [s for s in D.index if s in groups.index and pd.notna(groups.loc[s])]
    if len(common) < 4:
        return {"R2": np.nan, "F_statistic": np.nan, "p_value": np.nan, "n_samples": len(common)}
    D2 = D.loc[common, common].values.astype(float) ** 2
    g = groups.loc[common].astype(str).values
    levels = np.unique(g)
    if len(levels) < 2 or min([(g == lev).sum() for lev in levels]) < 2:
        return {"R2": np.nan, "F_statistic": np.nan, "p_value": np.nan, "n_samples": len(common)}
    n = len(g)
    ss_total = D2.sum() / n / 2.0

    def calc_F(labels):
        ss_within = 0.0
        for lev in np.unique(labels):
            idx = np.where(labels == lev)[0]
            if len(idx) > 0:
                ss_within += D2[np.ix_(idx, idx)].sum() / len(idx) / 2.0
        ss_between = ss_total - ss_within
        df_between = len(np.unique(labels)) - 1
        df_within = n - len(np.unique(labels))
        if df_between <= 0 or df_within <= 0 or ss_within <= 0:
            return np.nan, ss_between
        F = (ss_between / df_between) / (ss_within / df_within)
        return F, ss_between

    F_obs, ssb = calc_F(g)
    if not np.isfinite(F_obs):
        return {"R2": np.nan, "F_statistic": np.nan, "p_value": np.nan, "n_samples": n}
    perm_F = []
    for _ in range(n_perm):
        pg = rng.permutation(g)
        Fp, _ = calc_F(pg)
        if np.isfinite(Fp):
            perm_F.append(Fp)
    p = (np.sum(np.array(perm_F) >= F_obs) + 1) / (len(perm_F) + 1) if perm_F else np.nan
    R2 = ssb / ss_total if ss_total > 0 else np.nan
    return {"R2": R2, "F_statistic": F_obs, "p_value": p, "n_samples": n}

def betadisper_test(D: pd.DataFrame, groups: pd.Series) -> Dict[str, float]:
    coords = pcoa_from_distance(D)
    common = [s for s in coords.index if s in groups.index and pd.notna(groups.loc[s])]
    if len(common) < 4:
        return {"statistic": np.nan, "p_value": np.nan}
    df = coords.loc[common, ["PCoA1", "PCoA2"]].copy()
    df["group"] = groups.loc[common].astype(str)
    dists = []
    for lev, sub in df.groupby("group"):
        center = sub[["PCoA1", "PCoA2"]].mean().values
        dist = np.sqrt(((sub[["PCoA1", "PCoA2"]].values - center) ** 2).sum(axis=1))
        dists.append(dist)
    dists2 = [x for x in dists if len(x) >= 2]
    if len(dists2) < 2:
        return {"statistic": np.nan, "p_value": np.nan}
    try:
        stat, p = stats.kruskal(*dists2)
    except Exception:
        stat, p = np.nan, np.nan
    return {"statistic": stat, "p_value": p}
