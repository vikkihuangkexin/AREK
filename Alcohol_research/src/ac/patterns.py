"""Case-only fixed-panel PAM; independent cohorts and auditable OOB stability."""
import argparse,itertools,subprocess,shutil
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist,squareform
from sklearn.metrics import silhouette_score,normalized_mutual_info_score
from sklearn.decomposition import PCA
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

def matching(j):
    ri,ci=linear_sum_assignment(-j);best=j[ri,ci].sum()
    for p in itertools.permutations(range(len(j))):
        if np.isclose(sum(j[a,b] for a,b in enumerate(p)),best,rtol=0,atol=1e-12):
            return {b:a for a,b in enumerate(p) if j[a,b]>0}
    raise AssertionError("assignment")
def summarize(out,x,m):
    B=200;n=len(x);d=squareform(pdist(x))
    draw=np.loadtxt(out/"bootstrap_draw_indices.tsv",dtype=int,ndmin=2)-1
    states=read(out/"bootstrap_status.tsv")
    metrics=[];freqs=[];jrecords=[];maps=[]
    for k in range(2,7):
        orig=read(out/f"original_K{k}.tsv").cluster.to_numpy(int)-1
        boot=pd.read_csv(out/f"bootstrap_labels_K{k}.tsv.gz",sep="\t",header=None).to_numpy(float)-1
        js=np.zeros((B,k));counts=np.zeros((n,k),int);attempt=np.zeros(n,int)
        matched=np.zeros(n,int);failed=np.zeros(n,int);unmatched=np.zeros(n,int);valid=0
        for b in range(B):
            ins=np.unique(draw[b]);oob=np.setdiff1d(np.arange(n),ins);attempt[oob]+=1
            state=states[(states.K==k)&(states.attempt==b+1)].status.iloc[0]
            if state!="valid" or not np.isfinite(boot[b]).all():failed[oob]+=1;continue
            valid+=1;lab=boot[b].astype(int);j=np.zeros((k,k))
            for a in range(k):
                for c in range(k):
                    u=orig[ins]==a;v=lab[ins]==c
                    j[a,c]=(u&v).sum()/(u|v).sum() if (u|v).any() else 0.
            mp=matching(j)
            for c in range(k):
                a=mp.get(c)
                maps.append(dict(K=k,attempt=b+1,bootstrap_cluster=c+1,original_cluster=a+1 if a is not None else np.nan,Jaccard=j[a,c] if a is not None else 0))
                if a is not None:js[b,a]=j[a,c]
            for c in range(k):
                who=oob[lab[oob]==c];a=mp.get(c)
                if a is None:unmatched[who]+=1
                else:counts[who,a]+=1;matched[who]+=1
        sil=silhouette_score(d,orig,metric="precomputed");size=np.bincount(orig,minlength=k)
        gate=bool(size.min()>=max(20,int(np.ceil(.05*n))) and sil>=.25 and js.mean(0).min()>=.75 and valid>=190)
        nmis={f:normalized_mutual_info_score(m[f].astype(str),orig) for f in ["Study","analysis_set","assay_analysis","material_analysis"]}
        aligned=max(nmis.values())>=.8
        metrics.append(dict(K=k,n=n,min_cluster_size=int(size.min()),required_min_size=max(20,int(np.ceil(.05*n))),
            silhouette=sil,min_mean_Jaccard=js.mean(0).min(),algorithm_valid_attempts=valid,attempts=B,
            numerical_gates_pass=gate,source_highly_aligned=aligned,all_gates_pass=gate and not aligned,
            **{"NMI_"+f:v for f,v in nmis.items()}))
        for b in range(B):
            for a in range(k):jrecords.append(dict(K=k,attempt=b+1,cluster=a+1,Jaccard=js[b,a]))
        for i,sid in enumerate(x.index):
            assert attempt[i]==matched[i]+failed[i]+unmatched[i]
            for a in range(k):
                freqs.append(dict(sample_id=sid,K=k,original_cluster=int(orig[i]+1),cluster=a+1,
                    count=counts[i,a],OOB_attempts=attempt[i],matched_attempts=matched[i],
                    failed_attempts=failed[i],unmatched_attempts=unmatched[i],
                    assignment_stability_frequency=counts[i,a]/attempt[i] if attempt[i] else np.nan,
                    conditional_matched_frequency=counts[i,a]/matched[i] if matched[i] else np.nan))
        for f in nmis:
            save(pd.crosstab(pd.Series(orig+1,index=m.index,name="cluster"),m[f]),out/f"K{k}_cluster_by_{f}.tsv",True)
    met=pd.DataFrame(metrics);passing=met[met.all_gates_pass]
    chosen=None if passing.empty else int(passing[passing.silhouette.ge(passing.silhouette.max()-.01)].K.min())
    save(met,out/"all_K_metrics.tsv");save(pd.DataFrame(freqs),out/"OOB_assignment_stability_frequency.tsv.gz")
    save(pd.DataFrame(jrecords),out/"all_bootstrap_cluster_Jaccard.tsv.gz")
    save(pd.DataFrame(maps),out/"bootstrap_one_to_one_label_mapping.tsv.gz")
    return dict(status="stable_discrete_partition" if chosen else "continuous_only_no_K_passes",
        chosen_K=chosen,n_learning_working_units=n,bootstrap_attempts=200,
        source_aligned_numerically_passing_K=met.loc[met.numerical_gates_pass&met.source_highly_aligned,"K"].tolist(),
        assignment_output="OOB_stability_frequency_not_clinical_response_probability")
def run_node(c,s,p,cache):
    start=time.time();out=OUTPUT/c/p.split("__")[0]/s/"patterns"/p/"PAM_completion"
    if (out/"completion_receipt.json").exists():return json.loads((out/"completion_receipt.json").read_text())
    _,y,_=panel(p);m=representatives(rows(c,y.index,s,cases_only=True),"20261030::");x=y.loc[m.index]
    out.mkdir(parents=True,exist_ok=True)
    dep=dict(matrix_hash=framehash(x),sample_order=m.index.tolist(),
        study_order=m.Study.tolist(),working_partition=m.working_cluster.astype(str).tolist(),
        algorithm="PAM_Euclidean_unscaled_fixed_CLR",B=200,engine_sha256=sha(Path(__file__).with_name("pattern_completion.R")))
    ident=hashlib.sha256(json.dumps(dep,sort_keys=True).encode()).hexdigest()
    seeds=json.loads((Path(__file__).parents[2]/"configs/pattern_seeds.json").read_text())
    seed=int(seeds["::".join([c,s,p])]);dep["seed"]=seed
    if len(m)<40 or not np.isfinite(x.to_numpy()).all() or x.var().sum()<1e-12:
        r=dict(status="not_estimable",reason="fewer_than_40_case_working_units_for_two_min20_clusters" if len(m)<40 else "nonfinite_or_constant_CLR",
               n_learning_working_units=len(m),bootstrap_attempts=0,comparison=c,scope=s,panel=p,dependency=dep)
        dump(r,out/"completion_receipt.json");return r
    if ident in cache:
        src=cache[ident];shutil.copytree(src,out,dirs_exist_ok=True)
        r=json.loads((out/"completion_receipt.json").read_text())
        r.update(comparison=c,scope=s,panel=p,origin="verified_exact_case_geometry_reuse",
                 reused_from=str(src.relative_to(OUTPUT)))
        dump(r,out/"completion_receipt.json");return r
    save(x,out/"learning_CLR.tsv",True);save(m[["sample_id","Study","analysis_set","assay_analysis","material_analysis","working_cluster"]],out/"learning_metadata.tsv")
    pc=PCA(n_components=min(5,x.shape[1]-1,len(x)-1),svd_solver="full").fit(x)
    for i in range(len(pc.components_)):
        if pc.components_[i,np.argmax(abs(pc.components_[i]))]<0:pc.components_[i]*=-1
    save(pd.DataFrame(pc.transform(x),index=x.index,columns=["PC"+str(i+1) for i in range(pc.n_components_)]),out/"case_PCA_scores.tsv.gz",True)
    save(pd.DataFrame(pc.components_.T,index=x.columns),out/"case_PCA_loadings.tsv",True)
    subprocess.run(["Rscript",str(Path(__file__).with_name("pattern_completion.R")),str(out),str(seed)],check=True)
    r=summarize(out,x,m)
    r.update(comparison=c,scope=s,panel=p,origin="new_independent_case_PAM",dependency=dep,
             elapsed_seconds=time.time()-start,versions=versions())
    dump(r,out/"completion_receipt.json");cache[ident]=out
    print("PATTERN_COMPLETED",c,s,p,len(m),r["status"],"seconds",round(time.time()-start,2),flush=True)
    return r
