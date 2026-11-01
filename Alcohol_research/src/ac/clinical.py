"""Batch all frozen panels; support-specific models, source checks and joint BH."""
from .io import *
from .stats import bh
from .interaction import fit_group
import argparse, time, resource

def fit_node(m,y,v,config,context,dest,extra=None,source_checks=True):
    start=time.time();dest.mkdir(parents=True,exist_ok=True)
    coef=[];overall=[];contrasts=[];source=[];mapping=[];designs=[]
    patterns={}
    for feature in y:
        mask=y[feature].notna()
        if extra is not None:mask &= extra.notna().all(axis=1)
        patterns.setdefault(np.packbits(mask.to_numpy()).tobytes(),(mask,[]))[1].append(feature)
    for mask,features in patterns.values():
        mm=m.loc[y.index[mask]];yy=y.loc[mm.index,features];vv=v.loc[mm.index]
        ex=extra.loc[mm.index].to_numpy() if extra is not None else None
        sid=hash_ids(mm.index)
        mapping.extend(dict(feature_id=f,support_hash=sid) for f in features)
        save(mm[["analysis_group","working_cluster","Study"]],dest/"supports"/(sid+".tsv.gz"),True)
        if len(mm)<2:
            for f in features:
                inf=dict(feature_id=f,n_samples=len(mm),n_Alcohol=int(mm.analysis_group.eq("Alcohol").sum()),
                         n_Health=int(mm.analysis_group.eq("Health").sum()),sample_hash=sid,working_unit_hash=hash_ids(mm.working_cluster),
                         working_units_Alcohol=mm.loc[mm.analysis_group.eq("Alcohol"),"working_cluster"].nunique(),
                         working_units_Health=mm.loc[mm.analysis_group.eq("Health"),"working_cluster"].nunique(),n_Study=mm.Study.nunique(),n_two_arm_Study=0,
                         n_clusters=mm.working_cluster.nunique(),design_rank=0,residual_df=0,t_df=np.nan,covariance="not_fitted",
                         design_hash="NA",reason="no_response_support")
                overall.append(dict(**inf,P_interaction=np.nan,Wald_F=np.nan,interaction_df=np.nan))
                for model in ["subset_M1","clinical_main","interaction"]:
                    coef.append(dict(**inf,model_id=model,term="Alcohol",beta=np.nan,SE=np.nan,CI_low=np.nan,CI_high=np.nan,raw_P=np.nan))
            continue
        cr,ov,ct,states=fit_group(mm,yy.to_numpy(),vv,config,ex,make_curves=config["kind"]=="categorical")
        for z,target in [(cr,coef),(ov,overall),(ct,contrasts)]:
            if len(z):
                z["feature_id"]=z.feature_index.map(dict(enumerate(features)));z=z.drop(columns="feature_index")
                target.extend(z.to_dict("records"))
        for stage,s in states.items():
            designs.append(dict(support_hash=sid,model_id=stage,terms=json.dumps(s["terms"]),**s["info"]))
        # Only supported real two-arm sources are fitted, retaining all attempted statuses.
        if source_checks:
            for study,sm in mm.groupby("Study"):
                if sm.analysis_group.nunique()<2:continue
                eligible_features=ov[ov.P_interaction.notna()].feature_index.tolist()
                if not eligible_features:continue
                sv=vv.loc[sm.index];sy=yy.loc[sm.index].iloc[:,eligible_features]
                sx=extra.loc[sm.index].to_numpy() if extra is not None else None
                sc,so,_,_=fit_group(sm,sy.to_numpy(),sv,config,sx)
                ints=sc[(sc.model_id=="interaction")&sc.term.str.startswith("Alcohol:")].copy()
                ints["feature_id"]=ints.feature_index.map({i:features[f] for i,f in enumerate(eligible_features)})
                po=so.set_index("feature_index").P_interaction
                ints["P_interaction"]=ints.feature_index.map(po);ints["source_Study"]=study
                source.extend(ints.drop(columns="feature_index").to_dict("records"))
    C=pd.DataFrame(coef);O=pd.DataFrame(overall);T=pd.DataFrame(contrasts);S=pd.DataFrame(source)
    for d in [C,O,T,S]:
        if len(d):
            for k,val in context.items():d[k]=val
    O["q_interaction"]=bh(O.P_interaction);O["fixed_family_size"]=len(y.columns)
    C["q_Alcohol"]=np.nan
    for stage in C.model_id.unique():
        z=C.model_id.eq(stage)&C.term.eq("Alcohol")
        C.loc[z,"q_Alcohol"]=bh(C.loc[z,"raw_P"])
    save(C,dest/"coefficients.tsv.gz");save(O,dest/"interaction_tests.tsv.gz")
    if len(T):save(T,dest/"category_effects.tsv.gz")
    if len(S):save(S,dest/"within_Study_interactions.tsv.gz")
    save(pd.DataFrame(mapping),dest/"feature_support_index.tsv.gz")
    save(pd.DataFrame(designs),dest/"design_manifest.tsv")
    receipt=dict(context=context,features=len(y.columns),response_support_patterns=len(patterns),
                 eligible_interactions=int(O.P_interaction.notna().sum()),candidates=int(O.q_interaction.lt(.2).sum()),
                 elapsed_seconds=time.time()-start,engine_sha256=sha(Path(__file__).with_name("interaction.py")),
                 sample_population_hash=hash_ids(m.index),clinical_value_hash=hashlib.sha256(v.to_csv().encode()).hexdigest())
    dump(receipt,dest/"receipt.json")
    if context["coverage"]!="trait_coverage":
        print("NODE",context,receipt["eligible_interactions"],receipt["candidates"],round(receipt["elapsed_seconds"],2),flush=True)
    return O,C

def run(a):
    ROOT=Path(a.clinical_root)
    INPUT=Path(a.input_root)
    start=time.time()
    configs=json.loads((Path(__file__).parents[2]/"configs/clinical_registry.json").read_text())
    m=read(ROOT/"reproducibility/inputs/metadata.tsv.gz",index_col=0)
    for col in ["force_cluster_source","explicit_disease_reference"]:m[col]=m[col].astype(str).str.lower().eq("true")
    vv=read(ROOT/"harmonised_clinical_values.tsv.gz",index_col=0,low_memory=False)
    if a.pilot:configs=[x for x in configs if x["variable"] in ["age","sex","BMI"] and x["scope"]=="frozen" and x["body"]=="gut"]
    if a.variable:configs=[x for x in configs if x["variable"]==a.variable]
    if a.body:configs=[x for x in configs if x["body"]==a.body]
    if a.scope:configs=[x for x in configs if x["scope"]==a.scope]
    filtered=bool(a.variable or a.body or a.panel or a.scope or a.outcome or a.coverage)
    if filtered and not a.output_prefix and not a.pilot:
        ap.error("Filtered runs require --output-prefix; never overwrite the full-run aggregate.")
    prefix=Path(a.output_prefix).resolve() if a.output_prefix else ROOT/("validation/pilot" if a.pilot else "interaction_models")
    if filtered and prefix==(ROOT/"interaction_models").resolve():ap.error("Filtered output cannot be the frozen full-run directory.")
    alltests=[];allcoef=[];members=[]
    cache={}
    def load(p):
        if str(p) not in cache:cache[str(p)]=read(p,index_col=0)
        return cache[str(p)]
    for conf in configs:
        var=conf["variable"];body=conf["body"];scope=conf["scope"]
        mm=m[m.body.eq(body)&vv[var].notna()].copy()
        if scope!="frozen":mm=mm[~mm.explicit_disease_reference]
        v=vv.loc[mm.index,var]
        members.extend(dict(variable=var,body=body,scope=scope,sample_id=sid,clinical_value=v.loc[sid],
            analysis_group=row.analysis_group,working_cluster=row.working_cluster,Study=row.Study) for sid,row in mm.iterrows())
        jobs=[]
        for panel in sorted((INPUT/"reproducibility/inputs/panels").glob(body+"__*")):
            if a.panel and panel.name!=a.panel:continue
            if a.pilot and panel.name not in ["gut__common_16S","gut__full_source_reference"]:continue
            fam=read(panel/"family.tsv");fs=fam[fam.selected.astype(str).str.lower().eq("true")].taxon_key.tolist()
            jobs.append(("taxa",panel.name,"base",load(panel/"CLR_unadjusted.tsv.gz").reindex(columns=fs),"frozen CLR",None))
            if not a.pilot:
                jobs.append(("alpha",panel.name,"base",load(INPUT/"reproducibility/inputs/diversity"/panel.name/"alpha_per_specimen.tsv.gz")[["Shannon","Simpson","Observed_genera"]],"original alpha",None))
        if not a.pilot and (not a.panel or a.panel==body):
            fdir=INPUT/"reproducibility/inputs/functions"/body
            family=read(fdir/"family.tsv");keys=family[family.selected.astype(str).str.lower().eq("true")].taxon_key.tolist()
            fy=load(fdir/"contribution.tsv.gz").reindex(columns=keys)
            globalc=load(fdir/"global_coverage.tsv.gz")[["original_abundance_coverage"]]
            known=load(fdir/"known_mass.tsv.gz")
            jobs.extend([("function",body,"base",fy,"original contribution fraction",None),
                         ("function",body,"global_coverage",fy,"original contribution fraction",globalc),
                         ("function",body,"trait_coverage",fy,"original contribution fraction",known)])
        for outcome,panel,coverage,y,transformation,extra in jobs:
            if a.outcome and outcome!=a.outcome:continue
            if a.coverage and coverage!=a.coverage:continue
            y=y.reindex(mm.index)
            context=dict(variable=var,body=body,scope=scope,outcome=outcome,panel=panel,coverage=coverage,transformation=transformation,
                         unit=conf["unit"],center=conf["center"],scale=conf["scale"],reference=conf["reference"])
            dest=prefix/scope/body/var/outcome/panel/coverage
            if coverage=="trait_coverage":
                oo=[];cc=[]
                for feature in y:
                    ex=extra[[feature]].reindex(mm.index)
                    sub=dest/("trait_"+hashlib.sha256(feature.encode()).hexdigest()[:16])
                    o,c=fit_node(mm,y[[feature]],v,conf,context,sub,ex,source_checks=False);oo.append(o);cc.append(c)
                O=pd.concat(oo,ignore_index=True);C=pd.concat(cc,ignore_index=True)
                O["q_interaction"]=bh(O.P_interaction);O["fixed_family_size"]=len(y.columns)
                for stage in C.model_id.unique():
                    z=C.model_id.eq(stage)&C.term.eq("Alcohol")
                    C.loc[z,"q_Alcohol"]=bh(C.loc[z,"raw_P"])
                save(O,dest/"interaction_tests.tsv.gz");save(C,dest/"coefficients.tsv.gz")
            else:
                ex=extra.reindex(mm.index) if extra is not None else None
                O,C=fit_node(mm,y,v,conf,context,dest,ex,source_checks=coverage=="base" and not a.pilot)
            alltests.append(O);allcoef.append(C)
    if not alltests:ap.error("No nodes matched the requested filters.")
    O=pd.concat(alltests,ignore_index=True);C=pd.concat(allcoef,ignore_index=True)
    O["q_interaction_joint_variables"]=np.nan;O["joint_fixed_family_size"]=0
    registry=json.loads((Path(__file__).parents[2]/"configs/clinical_registry.json").read_text())
    for group,idx in O.groupby(["scope","body","outcome","panel","coverage"]).groups.items():
        required={x["variable"] for x in registry if x["scope"]==group[0] and x["body"]==group[1]}
        available=set(O.loc[idx,"variable"])
        complete=required==available
        O.loc[idx,"joint_BH_status"]="complete_frozen_family" if complete else "incomplete_variable_family"
        if complete:
            O.loc[idx,"q_interaction_joint_variables"]=bh(O.loc[idx,"P_interaction"])
            O.loc[idx,"joint_fixed_family_size"]=len(idx)
    save(O,prefix/"all_interaction_tests.tsv.gz");save(C,prefix/"all_model_coefficients.tsv.gz")
    save(pd.DataFrame(members),prefix/"subset_membership.tsv.gz" if a.output_prefix or a.pilot else ROOT/"subset_membership.tsv.gz")
    groups=O.groupby(["scope","body","variable","outcome","panel","coverage"],dropna=False)
    summary=groups.agg(fixed_family_size=("feature_id","size"),n_estimable=("P_interaction","count"),
        n_candidates=("q_interaction",lambda x:int(x.lt(.2).sum())),n_joint_candidates=("q_interaction_joint_variables",lambda x:int(x.lt(.2).sum()))).reset_index()
    summary.loc[summary.n_estimable.eq(0),["n_candidates","n_joint_candidates"]]=np.nan
    save(summary,prefix/"node_summary.tsv")
    dump(dict(status="PASS",nodes=len(summary),feature_tests=len(O),estimable=int(O.P_interaction.notna().sum()),
              total_seconds=time.time()-start,peak_rss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
              pilot=a.pilot),prefix/"run_receipt.json")
    print("RUN_COMPLETE",len(summary),len(O),round(time.time()-start,2),"seconds",flush=True)

