"""Explicit-design HC3/CR1 M1 regression and frozen multiplicity families."""
from .stats import *

def make_design(m, model="M1", extra=None):
    n,names=nuisance(m,model,extra)
    ni,keep=independent(n)
    a=m.analysis_group.eq("Alcohol").to_numpy(float)
    resid=a-ni@(np.linalg.pinv(ni)@a) if len(m) else a
    identifiable=np.dot(resid,resid)>1e-9
    x=np.column_stack([ni,a])
    aliases=pd.DataFrame(dict(candidate_column=names,retained=[j in keep for j in range(len(names))],
        reason=["independent" if j in keep else "aliased_or_constant_in_full_nuisance_space" for j in range(len(names))]))
    return ni,x,[names[j] for j in keep]+["Alcohol"],aliases,resid,identifiable

def fit(m,y,model="M1",extra=None,design_dir=None):
    y=np.asarray(y,float); n,p=y.shape
    a=m.analysis_group.eq("Alcohol").to_numpy(float)
    codes,unique=pd.factorize(m.working_cluster.astype(str),sort=True); G=len(unique)
    cluster=bool(m.get("force_cluster_source",pd.Series(False,index=m.index)).any() or G<n)
    units={g:m.loc[m.analysis_group.eq(g),"working_cluster"].nunique() for g in ["Alcohol","Health"]}
    info=dict(n_samples=n,n_Alcohol=int(a.sum()),n_Health=int(n-a.sum()),
        working_units_Alcohol=units["Alcohol"],working_units_Health=units["Health"],n_clusters=G,
        n_studies=m.Study.nunique(),sample_hash=hash_ids(m.index),
        covariance="working_cluster_CR1" if cluster else "HC3",inference_eligible=False,
        design_rank=0,residual_df=0,t_df=np.nan,design_hash="NA")
    out=pd.DataFrame({f:np.full(p,np.nan) for f in ["beta","SE","CI_low","CI_high","raw_p"]})
    def finish(reason):
        for k,v in info.items(): out[k]=v
        out["reason"]=reason
        return out, None
    if n<2 or m.analysis_group.nunique()<2: return finish("single_arm_or_empty")
    ni,x,names,aliases,aresid,identifiable=make_design(m,model,extra)
    info["residual_exposure_SS"]=float(aresid@aresid)
    if not identifiable: return finish("Alcohol_in_full_required_nuisance_span")
    k=x.shape[1];df=n-k;info.update(design_rank=k,residual_df=df,t_df=G-1 if cluster else df)
    info["design_hash"]=hashlib.sha256(x.tobytes()+hash_ids(m.index).encode()).hexdigest()
    pinv=np.linalg.pinv(x);coef=pinv@y; e=y-x@coef;w=pinv[-1]
    out["beta"]=coef[-1]
    if design_dir is not None:
        d=Path(design_dir)/info["design_hash"]
        if not (d/"design.tsv.gz").exists():
            save(pd.DataFrame(x,index=m.index,columns=names),d/"design.tsv.gz",True)
            save(aliases,d/"aliases.tsv")
            save(pd.DataFrame(dict(sample_id=m.index,working_cluster=m.working_cluster.values,
                residualized_Alcohol=aresid,residualized_Alcohol_SS=aresid**2)),
                d/"residual_exposure_information.tsv.gz")
    if df<=0: return finish("no_residual_df")
    h=np.sum(x*pinv.T,axis=1)
    if cluster and G<2: return finish("insufficient_clusters")
    if not cluster and (h>=1-1e-10).any(): return finish("unit_leverage_HC3_undefined")
    var=robust_covariance(pinv[-1:],e,h,codes,G,n,k,cluster)[:,0,0]
    se=np.sqrt(np.maximum(var,0));out["SE"]=se
    critical=stats.t.ppf(.975,info["t_df"])
    out["CI_low"]=out.beta-critical*se;out["CI_high"]=out.beta+critical*se
    eligible=min(units.values())>=(10 if model=="M2" else 5) and df>=(20 if model=="M2" else 10) and (not cluster or G>=20)
    good=eligible & np.isfinite(se) & (se>1e-12) & (np.ptp(y,axis=0)>1e-12)
    out["raw_p"]=np.where(good,2*stats.t.sf(np.abs(out.beta/se),info["t_df"]),np.nan)
    for key,value in info.items(): out[key]=value
    out["inference_eligible"]=good
    out["reason"]=np.where(good,"none",
        np.where(np.ptp(y,axis=0)<=1e-12,"constant_response",
                 "working_unit_df_cluster_threshold_or_undefined_variance"))
    state=dict(N=ni,X=x,pinv=pinv,h=h,w=w,codes=codes,G=G,cluster=cluster,
               residual_exposure=aresid)
    return out,state

def fit_matrix(m,y,model="M1",extra=None,B=0,context="",out_dir=None):
    """Group columns ONLY by identical actual response support; preserve all observations."""
    if B: raise ValueError("Wild bootstrap is outside these completed ordinary-P modules.")
    result=[]; manifests=[]
    patterns={}
    for col in y.columns:
        mask=y[col].notna().to_numpy()
        patterns.setdefault(np.packbits(mask).tobytes(),(mask,[]))[1].append(col)
    for ix,(mask,cols) in enumerate(patterns.values()):
        mm=m.loc[y.index[mask]];yy=y.loc[mm.index,cols]
        ex=extra.loc[mm.index] if extra is not None else None
        rr,state=fit(mm,yy.to_numpy(),model,ex,
                     Path(out_dir)/"designs" if out_dir else None)
        rr["taxon_key"]=cols;rr["empirical_p"]=np.nan
        rr["empirical_reason"]="not_scheduled" if not B else "M1_not_inferable"
        result.append(rr)
    r=pd.concat(result,ignore_index=True).set_index("taxon_key").reindex(y.columns).rename_axis("taxon_key").reset_index()
    if out_dir:
        dump(manifests,Path(out_dir)/"resampling_manifest.json")
    return r

def apply_families(r,family):
    r=r.merge(family,on="taxon_key",how="left",validate="many_to_one")
    r["selected"]=r.selected.astype(str).str.lower().eq("true")
    r["full_reference_BH_q"]=bh(r.raw_p);r["full_reference_empirical_BH_q"]=bh(r.empirical_p)
    r["BH_q"]=np.nan;r["empirical_BH_q"]=np.nan
    mask=r.selected
    r.loc[mask,"BH_q"]=bh(r.loc[mask,"raw_p"])
    r.loc[mask,"empirical_BH_q"]=bh(r.loc[mask,"empirical_p"])
    r["primary_family_size"]=int(mask.sum());r["full_reference_family_size"]=len(r)
    r["primary_family_valid_tests"]=int(r.loc[mask,"raw_p"].notna().sum())
    r["status"]=np.where(r.raw_p.notna(),"inferable","not_estimable")
    return r
