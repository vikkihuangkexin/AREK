# Extracted in-scope historical implementation; see docs/source_script_map.tsv.
from __future__ import annotations
from typing import Dict
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
def gamma_metric(counts: pd.DataFrame) -> Dict[str,float]:
    pooled = counts.sum(axis=1).astype(float)
    pos = pooled[pooled>0]
    total = pos.sum()
    p = pos/total if total>0 else pd.Series(dtype=float)
    return {"Observed": float(len(pos)), "Shannon": float(-(p*np.log(p)).sum()) if len(p) else np.nan, "Simpson": float(1-(p**2).sum()) if len(p) else np.nan}

def gamma_pairwise(counts: pd.DataFrame, groups: pd.Series, a: str, b: str, n_perm: int, seed=1) -> pd.DataFrame:
    ids = groups[groups.astype(str).isin([a,b])].dropna().index.intersection(counts.columns)
    g = groups.loc[ids].astype(str)
    if (g==a).sum()<2 or (g==b).sum()<2: return pd.DataFrame()
    obs_a = gamma_metric(counts.loc[:,g[g==a].index]); obs_b = gamma_metric(counts.loc[:,g[g==b].index])
    rng=np.random.default_rng(seed); perm_d={m:[] for m in obs_a}
    lab=g.values.copy()
    for _ in range(n_perm):
        pg=rng.permutation(lab)
        for m in obs_a:
            ga=gamma_metric(counts.loc[:,ids[pg==a]])[m]; gb=gamma_metric(counts.loc[:,ids[pg==b]])[m]
            perm_d[m].append(ga-gb)
    rows=[]
    for m in obs_a:
        d=obs_a[m]-obs_b[m]; arr=np.abs(np.array(perm_d[m])); p=(np.sum(arr>=abs(d))+1)/(len(arr)+1)
        rows.append({"metric":m,"group_a":a,"group_b":b,"gamma_a":obs_a[m],"gamma_b":obs_b[m],"effect_size":d,"p_value":p})
    return pd.DataFrame(rows)
