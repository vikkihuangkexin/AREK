"""Independent checks using only genuine frozen project data."""
from pathlib import Path
import os
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import pytest
from ac.io import read
from ac.inputs import Inputs
from ac.stats import bh
from ac.m1 import fit

@pytest.fixture
def input_root():
    configured = os.environ.get("AC_TEST_INPUT_ROOT")
    if not configured:
        pytest.skip("Set AC_TEST_INPUT_ROOT to an external original primary verification input.")
    return Path(configured).resolve()

def test_fixed_family_BH_real_pvalues(input_root):
    multipletests = pytest.importorskip("statsmodels.stats.multitest").multipletests
    root = input_root
    result = read(root/"01_Alcohol_vs_Health/gut/M1_all/taxa/gut__full_source_reference/all_results.tsv.gz")
    result = result[result.selected.astype(str).str.lower().eq("true")]
    p = result.raw_p.to_numpy()
    reference = multipletests(np.where(np.isfinite(p), p, 1), method="fdr_bh")[1]
    reference[~np.isfinite(p)] = np.nan
    np.testing.assert_allclose(bh(p), reference, equal_nan=True, rtol=1e-12, atol=1e-12)

def test_real_project_CR1_against_statsmodels(input_root):
    sm = pytest.importorskip("statsmodels.api")
    from ac.m1 import make_design
    data = Inputs(input_root)
    _, y, family = data.panel("gut__full_source_reference")
    key = next(x for x in family.taxon_key if x.endswith("::Haemophilus"))
    m = data.rows("01_Alcohol_vs_Health", y.index[y[key].notna()])
    values = y.loc[m.index, [key]].to_numpy()
    result, state = fit(m, values)
    _, x, _, _, _, _ = make_design(m)
    clusters = m.working_cluster.astype("category").cat.codes
    reference = sm.OLS(values[:, 0], x).fit(cov_type="cluster",
        cov_kwds={"groups": clusters, "use_correction": True}, use_t=True)
    np.testing.assert_allclose(result.beta.iloc[0], reference.params[-1], rtol=2e-7, atol=2e-9)
    np.testing.assert_allclose(result.SE.iloc[0], reference.bse[-1], rtol=2e-7, atol=2e-9)
    np.testing.assert_allclose(result.raw_p.iloc[0], reference.pvalues[-1], rtol=2e-7, atol=2e-9)
