"""Replay saved beta locations and complete actual vegan median dispersion."""
import argparse, subprocess
from scipy.spatial.distance import pdist,squareform
from .io import *
from .inputs import Inputs, representatives, versions
ROOT=Path(os.environ["AC_INPUT_ROOT"])
OUTPUT=Path(os.environ["AC_OUTPUT_ROOT"])
I=ROOT/"reproducibility/inputs"
_inputs=Inputs(ROOT)
metadata=_inputs.metadata
membership=_inputs.membership
panel=_inputs.panel
rows=_inputs.rows

from .stats import nuisance,independent
from .beta_kernel import permutations,kernel,test_kernel
def fast(G,N,a,perm,blocks):
    N,_=independent(N);Q=np.linalg.qr(N,mode="reduced")[0]
    q=a-Q@(Q.T@a);q/=np.linalg.norm(q);df=len(a)-N.shape[1]-1
    total=np.trace(G);null=total-np.sum(Q*(G@Q));ss=q@G@q
    if not all(np.allclose(N[b],N[b[0]],atol=1e-12) for b in blocks):
        return test_kernel(G,N,a,perm,blocks,force_general=True)
    active=np.concatenate([b for b in blocks if np.max(abs(q[b]))>1e-10])
    sub=G[np.ix_(active,active)];ff=[]
    for off in range(0,len(perm),128):
        v=q[perm[off:off+128][:,active]].T;s=np.sum(v*(sub@v),axis=0)
        ff.extend(s/(null-s)*df)
    ff=np.asarray(ff);F=ss/(null-ss)*df
    ref,rf=test_kernel(G,N,a,perm[:99],blocks)
    assert np.allclose(rf,ff[:99],atol=1e-7,rtol=1e-7)
    return dict(F=float(F),R2=float(ss/total),P=float(((ff>=F-1e-12).sum()+1)/(len(ff)+1)),
                B_actual=len(ff),full_kernel_max_error=float(np.max(abs(rf-ff[:99])))),ff
def run_node(compare,scope,pid,force_dispersion=False):
    start=time.time();body=pid.split("__")[0]
    dest=ROOT/compare/body/scope/"diversity"/pid
    old=read(dest/"beta_tests.tsv");valid=old[old.P.notna()]
    if valid.empty:return []
    out=OUTPUT/compare/body/scope/"diversity"/pid/"completion";receipt=out/"completion_receipt.json"
    if receipt.exists() and json.loads(receipt.read_text()).get("validation_version")=="vegan_marginal_observed_v2":
        return read(out/"completed_tests.tsv").to_dict("records")
    t,y,_=panel(pid)
    spec=dest/"permutation_specification.json"
    if spec.exists():
        recipe=json.loads(spec.read_text())
        m=read(dest/"formal_metadata.tsv.gz",index_col=0,keep_default_na=False)
        mm=rows(compare,m.index,scope);assert mm.index.equals(m.index)
        assert mm.analysis_group.equals(m.analysis_group)
        m=mm
        seed=int(recipe["seed"])
    else:
        # Original v3 samples/recipe must be restored before validating reuse.
        p=dest/"original_Biomni_permutation_files"/"permutation_specification.json"
        if not p.exists():raise RuntimeError("original_recipe_requires_restoration:"+str(dest))
        recipe=json.loads(p.read_text());seed=int(recipe["seed"])
        m=read(dest/"original_Biomni_permutation_files"/"metadata.tsv",index_col=0,keep_default_na=False)
        m=rows(compare,m.index,scope)
    perm,blocks,key=permutations(m,9999,seed)
    h=hashlib.sha256(perm.tobytes()).hexdigest()
    assert h==recipe.get("permutation_hash",recipe.get("permutation_array_sha256")),pid
    N,_=nuisance(m,"M1");N,_=independent(N);a=m.analysis_group.eq("Alcohol").to_numpy(float)
    assert m.analysis_group.nunique()==2 and min(m.analysis_group.value_counts())>=5
    assert any(np.unique(a[b]).size==2 for b in blocks)
    assert np.linalg.norm(a-N@(np.linalg.pinv(N)@a))>1e-9
    out.mkdir(parents=True,exist_ok=True)
    save(t.loc[m.index],out/"TSS.tsv.gz",True);save(y.loc[m.index],out/"CLR.tsv.gz",True)
    save(m,out/"metadata.tsv",True)
    save(pd.DataFrame(N,index=m.index,columns=["N"+str(j) for j in range(N.shape[1])]),out/"nuisance.tsv",True)
    save(pd.DataFrame(perm+1),out/"permutations_1based.tsv.gz")
    records=[];errors={}
    for metric in ["BrayCurtis","Aitchison"]:
        vals=t.loc[m.index] if metric=="BrayCurtis" else y.loc[m.index]
        G=kernel(squareform(pdist(vals,"braycurtis" if metric=="BrayCurtis" else "euclidean")))
        original,ff_old=fast(G,N,a,perm,blocks)
        expected=valid[valid.metric.eq(metric)].iloc[0]
        for c in ["F","R2","P"]:assert np.isclose(original[c],expected[c],atol=1e-7,rtol=1e-7),(compare,pid,c,original[c],expected[c])
        stored=dest/(metric+"_permutation_F.tsv.gz")
        if not stored.exists():stored=dest/"original_Biomni_permutation_files"/(metric+"_permutation_F.tsv.gz")
        ef=read(stored).iloc[:,0].to_numpy(float);assert len(ef)==9999 and np.allclose(ef,ff_old,atol=1e-7,rtol=1e-7)
        errors[metric]=float(np.max(abs(ef-ff_old)))
        commutes=all(np.allclose(N[b],N[b[0]],atol=1e-12) for b in blocks)
        r,ff=fast(G,N,a,np.argsort(perm,axis=1),blocks) if commutes else test_kernel(G,N,a,perm,blocks,force_general=True)
        for c in ["F","R2"]:assert np.isclose(r[c],expected[c],atol=1e-7,rtol=1e-7)
        save(pd.DataFrame(dict(F=ff)),out/(metric+"_location_Fperm.tsv.gz"))
        records.append(dict(comparison=compare,scope=scope,body=body,panel_id=pid,metric=metric,
            n=len(m),n_case=int(a.sum()),n_reference=int(len(m)-a.sum()),seed=seed,
            permutation_hash=h,reason="none",location_status="vegan_aligned_permutation_direction",
            original_location_P=float(expected.P),location_P_change=r["P"]-float(expected.P),
            permutation_convention="vegan_residual_kernel_P; equivalent_commuting_test_vector_inverse_P",**r))
    subprocess.run(["Rscript",str(Path(__file__).with_name("beta_completion.R")),str(out)],check=True)
    rr=read(out/"R_results.tsv").set_index("metric")
    for row in records:
        metric=row["metric"];rf=read(out/(metric+"_vegan_location_Fperm99.tsv.gz")).F.to_numpy()
        pf=read(out/(metric+"_location_Fperm.tsv.gz")).F.to_numpy()
        assert np.allclose(rf,pf[:99],atol=1e-7,rtol=1e-7),(compare,pid,metric,"vegan_location_mismatch")
        # vegan's marginal observed numerator is a difference of eigentruncated
        # residual inertias, not the direct one-dimensional quadratic form.
        # Retain both and use the actual vegan marginal observed statistic.
        vf=float(rr.loc[metric,"vegan_location_F"])
        assert np.isclose(vf,row["F"],atol=1e-6,rtol=1e-4),("vegan_observed_F_large_discrepancy",compare,pid,metric,vf,row["F"])
        row["original_location_F"]=row["F"];row["original_location_R2"]=row["R2"]
        row["observed_vegan_minus_direct_kernel_F"]=vf-row["F"]
        row["F"]=vf;row["R2"]=float(rr.loc[metric,"vegan_location_R2"])
        row["P"]=float(((pf>=vf-np.sqrt(np.finfo(float).eps)).sum()+1)/(len(pf)+1))
        row["location_P_change"]=row["P"]-row["original_location_P"]
        row["vegan_first99_F_max_error"]=float(np.max(abs(rf-pf[:99])))
        row.update(rr.loc[metric].to_dict())
        df=read(out/(metric+"_dispersion_Fperm.tsv.gz")).F.to_numpy()
        pv=((df>=float(row["dispersion_F"])-np.sqrt(np.finfo(float).eps)).sum()+1)/10000
        assert np.isclose(pv,row["dispersion_P"],atol=1e-4),("dispersion_p",pv,row["dispersion_P"])
        row["dispersion_status"]="completed_actual_vegan"
        previous=valid[valid.metric.eq(metric)].iloc[0]
        if pd.notna(previous.get("dispersion_P",np.nan)):
            assert np.isclose(row["dispersion_P"],previous.dispersion_P,atol=1e-7),("existing_dispersion",compare,pid,metric)
            row["dispersion_status"]="replayed_and_verified_existing_vegan"
    save(pd.DataFrame(records),out/"completed_tests.tsv")
    # Exact recipe/hash and F vectors are retained; huge text permutation matrix is reproducible.
    (out/"permutations_1based.tsv.gz").unlink()
    dump(dict(status="PASS",validation_version="vegan_marginal_observed_v2",comparison=compare,scope=scope,panel=pid,actual_B=9999,
              sample_order=m.index.tolist(),permutation_specification=recipe,
              full_saved_Fperm_max_errors=errors,vegan_location_validation_B=99,
              elapsed_seconds=time.time()-start,versions=versions()),receipt)
    print("BETA_COMPLETED",compare,scope,pid,len(m),"seconds",round(time.time()-start,2),flush=True)
    return records
