"""Shared fixed-family multiplicity, nuisance bases, and robust covariance."""
from .io import *
from scipy import stats, linalg

NUISANCE = ["Study", "analysis_set", "assay_analysis", "material_analysis", "technical_batch"]
NA = {"", "NA", "nan", "None", "unknown", "Unknown"}
COVS = ["age_model", "sex_model", "BMI_model"]
SEED = 20261030

def present(value):
    return pd.notna(value) and str(value) not in NA

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if not n:
        return p
    values = np.where(np.isfinite(p), p, 1)
    order = np.argsort(values, kind="stable")
    q = np.minimum.accumulate((values[order]*n/np.arange(1, n+1))[::-1])[::-1].clip(0, 1)
    result = np.empty(n)
    result[order] = q
    result[~np.isfinite(p)] = np.nan
    return result

def independent(x):
    keep = []
    q = np.empty((len(x), 0))
    for j in range(x.shape[1]):
        residual = x[:, j].copy()
        if q.shape[1]:
            residual -= q@(q.T@residual)
            residual -= q@(q.T@residual)
        norm = np.linalg.norm(residual)
        if norm > 1e-9:
            keep.append(j)
            q = np.column_stack([q, residual/norm])
    return x[:, keep], keep

def nuisance(m, model="M1", extra=None, clinical=False):
    """Preserve the two archived missing-category encodings without changing spans."""
    columns, names = [np.ones(len(m))], ["Intercept"]
    if model not in {"M0", "M0_CC"}:
        for field in NUISANCE:
            if field == "technical_batch":
                available = m[field].notna().any() if clinical else m[field].map(present).any()
                if not available:
                    continue
            values = m[field].fillna("NA").astype(str) if clinical else m[field].astype(str)
            if field == "technical_batch":
                values = m.Study.astype(str) + ":" + values
            for level in sorted(values.unique())[1:]:
                columns.append(values.eq(level).to_numpy(float))
                names.append(field+"="+level)
    if extra is not None:
        for field in extra:
            values = np.asarray(extra[field], float)
            columns.append(values-values.mean())
            names.append(field)
    return np.column_stack(columns), names

def robust_covariance(pinv_interest, residuals, leverage, codes, clusters, n, rank, clustered):
    """HC3 or working-cluster CR1 for arbitrary coefficients of interest."""
    if clustered:
        order = np.argsort(codes, kind="stable")
        cuts = np.r_[0, np.flatnonzero(np.diff(codes[order]))+1]
        scores = np.add.reduceat(
            (pinv_interest.T[:, :, None]*residuals[:, None, :])[order], cuts, axis=0)
        return np.einsum("gim,gjm->mij", scores, scores, optimize=True)*clusters/(clusters-1)*(n-1)/(n-rank)
    return np.einsum("in,nm,jn->mij", pinv_interest,
        (residuals/(1-leverage[:, None]))**2, pinv_interest, optimize=True)

def clr(t, close_first=False):
    values = t.div(t.sum(axis=1), axis=0) if close_first else t.copy()
    values = values.mask(values.eq(0), 1e-6)
    values = values.div(values.sum(axis=1), axis=0)
    logs = np.log(values)
    return logs.sub(logs.mean(axis=1), axis=0)
