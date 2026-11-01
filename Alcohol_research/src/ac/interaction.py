"""Same-support nested OLS with explicit HC3/CR1 interaction contrasts."""
from .stats import *

def encode(v,config):
    if config["kind"]=="continuous":
        c=((pd.to_numeric(v)-config["center"])/config["scale"]).to_numpy(float)[:,None]
        names=[config["variable"]]
    else:
        lev=[x for x in config["levels"] if x!=config["reference"]]
        c=np.column_stack([v.astype(str).eq(x).to_numpy(float) for x in lev]) if lev else np.empty((len(v),0))
        names=[config["variable"]+"="+x for x in lev]
    return c,names

def support(m,v,config):
    if len(m)<2 or m.analysis_group.nunique()<2:return "single_arm_or_empty",{}
    units={a:m.loc[m.analysis_group.eq(a),"working_cluster"].nunique() for a in ["Alcohol","Health"]}
    if min(units.values())<5:return "fewer_than_5_working_units_in_arm",{}
    extra={}
    if config["kind"]=="continuous":
        a=pd.to_numeric(v[m.analysis_group.eq("Alcohol")]);h=pd.to_numeric(v[m.analysis_group.eq("Health")])
        lo=max(a.min(),h.min());hi=min(a.max(),h.max());extra.update(common_low=lo,common_high=hi)
        if not hi>lo:return "no_continuous_common_range",extra
        for arm in ["Alcohol","Health"]:
            vi=pd.to_numeric(v[m.analysis_group.eq(arm)]);ids=vi[vi.between(lo,hi)].index
            if m.loc[ids,"working_cluster"].nunique()<5 or vi.loc[ids].nunique()<5:return "common_range_below_5_units_or_5_distinct_values",extra
    else:
        if len(config["levels"])<2:return "constant_covariate",{}
        for lev in config["levels"]:
            for arm in ["Alcohol","Health"]:
                ids=m.index[m.analysis_group.eq(arm)&v.astype(str).eq(lev)]
                if m.loc[ids,"working_cluster"].nunique()<5:return "category_by_group_cell_below_5_working_units",{}
    return "none",extra

def fit_one_stage(m,y,c,c_names,stage,extra=None):
    n=len(m);pout=y.shape[1]
    rawN,names=nuisance(m,clinical=True);N,keep=independent(rawN);names=[names[j] for j in keep]
    if extra is not None:
        ee=np.asarray(extra,float)
        if ee.ndim==1:ee=ee[:,None]
        N,nkeep=independent(np.column_stack([N,ee-ee.mean(axis=0)]))
        # Numeric names are provenance only, full required space remains preserved.
        names=[f"nuisance_{j}" for j in range(N.shape[1])]
    a=m.analysis_group.eq("Alcohol").to_numpy(float)
    blocks=[N,a[:,None]];term_names=["Alcohol"]
    if stage!="subset_M1":blocks.append(c);term_names+=c_names
    if stage=="interaction":blocks.append(a[:,None]*c);term_names+=["Alcohol:"+s for s in c_names]
    X=np.column_stack(blocks);kN=N.shape[1];interest=np.arange(kN,X.shape[1])
    rank=np.linalg.matrix_rank(X,tol=1e-9);df=n-rank
    pinv=np.linalg.pinv(X,rcond=1e-12);co=pinv@y;e=y-X@co;h=(X*pinv.T).sum(axis=1)
    codes,unq=pd.factorize(m.working_cluster.astype(str),sort=True);G=len(unq)
    cluster=bool(m.force_cluster_source.any() or G<n);tdf=G-1 if cluster else df
    eligible=df>=10 and (not cluster or G>=20) and m.analysis_group.nunique()==2 and min(
        m.loc[m.analysis_group.eq(arm),"working_cluster"].nunique() for arm in ["Alcohol","Health"])>=5
    reason="none" if eligible else "original_unit_df_cluster_threshold"
    if not cluster and np.any(h>=1-1e-10):eligible=False;reason="unit_leverage_HC3_undefined"
    B=pinv[interest]
    cov=np.full((pout,len(interest),len(interest)),np.nan)
    if df>0 and G>1 and (cluster or np.all(h<1-1e-10)):
        cov=robust_covariance(B,e,h,codes,G,n,rank,cluster)
    projected=pinv@X
    estim=np.max(np.abs(projected[interest]-np.eye(X.shape[1])[interest]),axis=1)<1e-7
    ss_c=np.linalg.norm(c-N@(np.linalg.pinv(N)@c)) if c.size else 0
    gate_c="none" if ss_c>1e-9 and c.shape[1]>0 else "clinical_variable_in_required_nuisance_span"
    inf=dict(n_samples=n,n_Alcohol=int(a.sum()),n_Health=int(n-a.sum()),working_units_Alcohol=m.loc[m.analysis_group.eq("Alcohol"),"working_cluster"].nunique(),
        working_units_Health=m.loc[m.analysis_group.eq("Health"),"working_cluster"].nunique(),n_Study=m.Study.nunique(),
        n_two_arm_Study=sum(g.analysis_group.nunique()==2 for _,g in m.groupby("Study")),n_clusters=G,design_rank=rank,
        residual_df=df,t_df=tdf,covariance="CR1" if cluster else "HC3",sample_hash=hash_ids(m.index),
        working_unit_hash=hash_ids(m.working_cluster),design_hash=hashlib.sha256(X.tobytes()+hash_ids(m.index).encode()).hexdigest())
    return dict(beta=co[interest].T,cov=cov,estimable=estim,terms=term_names,info=inf,eligible=eligible,
                reason=reason,clinical_gate=gate_c,X=X,pinv=pinv,full_co=co,e=e,h=h,N=N,interest=interest)

def contrasts_for(m,v,config,state,j,reason):
    rows=[];terms=state["terms"];b=state["beta"][j];V=state["cov"][j]
    if state["info"]["t_df"]<=0:return rows
    aidx=terms.index("Alcohol");ii=[i for i,t in enumerate(terms) if t.startswith("Alcohol:")]
    if config["kind"]=="continuous":
        ar=pd.to_numeric(v[m.analysis_group.eq("Alcohol")]);hr=pd.to_numeric(v[m.analysis_group.eq("Health")])
        lo=max(ar.min(),hr.min());hi=min(ar.max(),hr.max())
        grid=np.linspace(lo,hi,100) if np.isfinite(lo+hi) and hi>lo else []
        levels=[(float(x),(x-config["center"])/config["scale"]) for x in grid]
    else:levels=[(lev,lev) for lev in config["levels"]]
    for val,encoded in levels:
        L=np.zeros(len(terms));L[aidx]=1
        if config["kind"]=="continuous":
            if ii:L[ii[0]]=encoded
        elif val!=config["reference"]:
            name="Alcohol:"+config["variable"]+"="+str(val)
            if name in terms:L[terms.index(name)]=1
        est=bool(np.all(state["estimable"][np.flatnonzero(L)]))
        variance=float(L@V@L);se=np.sqrt(max(variance,0)) if np.isfinite(variance) else np.nan
        effect=float(L@b);critical=stats.t.ppf(.975,state["info"]["t_df"])
        good=reason=="none" and est and np.isfinite(se) and se>1e-12
        rows.append(dict(clinical_value=val,effect=effect,SE=se if good else np.nan,
            CI_low=effect-critical*se if good else np.nan,CI_high=effect+critical*se if good else np.nan,
            contrast_P=2*stats.t.sf(abs(effect/se),state["info"]["t_df"]) if good else np.nan,
            estimable=good,contrast_reason="none" if good else reason if reason!="none" else "contrast_not_estimable",
            **state["info"]))
    return rows

def fit_group(m,y,v,config,extra=None,make_curves=False):
    """All columns in y must have the same nonmissing support."""
    y=np.asarray(y,float);c,cn=encode(v,config);n,p=y.shape
    clinical_reason,common=support(m,v,config)
    coefs=[];overall=[];curve=[];states={}
    for stage in ["subset_M1","clinical_main","interaction"]:
        state=fit_one_stage(m,y,c,cn,stage,extra);states[stage]=state
        tdf=state["info"]["t_df"];critical=stats.t.ppf(.975,tdf) if tdf>0 else np.nan
        for j in range(p):
            for it,term in enumerate(state["terms"]):
                b=state["beta"][j,it];var=state["cov"][j,it,it];se=np.sqrt(max(var,0)) if np.isfinite(var) else np.nan
                base_good=state["eligible"] and state["estimable"][it] and np.isfinite(se) and se>1e-12 and np.ptp(y[:,j])>1e-12
                term_reason=state["reason"] if not state["eligible"] else "coefficient_not_estimable" if not state["estimable"][it] else "constant_response_or_undefined_variance" if not base_good else "none"
                if stage=="interaction" and term.startswith("Alcohol:") and (clinical_reason!="none" or state["clinical_gate"]!="none"):
                    base_good=False;term_reason=clinical_reason if clinical_reason!="none" else state["clinical_gate"]
                coefs.append(dict(feature_index=j,model_id=stage,term=term,beta=b,SE=se if base_good else np.nan,
                    CI_low=b-critical*se if base_good else np.nan,CI_high=b+critical*se if base_good else np.nan,
                    raw_P=2*stats.t.sf(abs(b/se),tdf) if base_good else np.nan,reason=term_reason,**state["info"]))
    state=states["interaction"];ii=[i for i,t in enumerate(state["terms"]) if t.startswith("Alcohol:")];q=len(ii)
    for j in range(p):
        reason=clinical_reason if clinical_reason!="none" else state["clinical_gate"] if state["clinical_gate"]!="none" else state["reason"]
        pp=np.nan;F=np.nan
        if reason=="none" and (not q or not np.all(state["estimable"][ii])):reason="whole_interaction_block_not_estimable"
        V=state["cov"][j][np.ix_(ii,ii)]
        if reason=="none" and (not np.isfinite(V).all() or np.linalg.matrix_rank(V)<q or np.any(np.diag(V)<=1e-24) or np.ptp(y[:,j])<=1e-12):reason="singular_robust_covariance_or_constant_response"
        if reason=="none":
            b=state["beta"][j,ii];F=float(b@np.linalg.pinv(V)@b/q);pp=float(stats.f.sf(F,q,state["info"]["t_df"]))
        overall.append(dict(feature_index=j,P_interaction=pp,Wald_F=F,interaction_df=q,reason=reason,**common,**state["info"]))
        if make_curves:curve.extend(dict(feature_index=j,**z) for z in contrasts_for(m,v,config,state,j,reason))
    return pd.DataFrame(coefs),pd.DataFrame(overall),pd.DataFrame(curve),states
