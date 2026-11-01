"""Fixed-family support-specific within-Study and LOSO M1, without hit deletion."""
import argparse
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

from . import m1 as models
def nodes():
    for path in sorted(list(ROOT.glob("0*/*/M1*/taxa/*/all_results.tsv.gz"))+
                       list(ROOT.glob("0*/*/M1*/functions/all_results.tsv.gz"))):
        c,b,s,kind=path.relative_to(ROOT).parts[:4]
        pid=path.parent.name if kind=="taxa" else b
        p=read(path);p=p[p.selected.astype(str).str.lower().eq("true")].copy()
        src=I/"panels"/pid/"CLR_unadjusted.tsv.gz" if kind=="taxa" else I/"functions"/b/"contribution.tsv.gz"
        y=read(src,index_col=0).reindex(columns=p.taxon_key)
        yield c,b,s,kind,pid,p,y
def components(studies,sets):
    remaining=set(studies);n=0
    while remaining:
        n+=1;todo=[remaining.pop()]
        while todo:
            st=todo.pop();linked=[t for t in remaining if sets[st]&sets[t]]
            for t in linked:remaining.remove(t);todo.append(t)
    return n
def one(c,b,s,kind,pid,p,y):
    out=OUTPUT/c/b/s/"source_robustness"/kind/pid
    mm=rows(c,y.index,s);y=y.loc[mm.index]
    config=dict(metadata_hash=framehash(mm),response_hash=framehash(y),primary_hash=framehash(p),
                implementation_hash=sha(Path(__file__)),model_hash=sha(Path(models.__file__)))
    if (out/"receipt.json").exists() and json.loads((out/"receipt.json").read_text())["config"]==config:
        return read(out/"robustness_classification.tsv")
    fits=[];sets={st:set(g.working_cluster.astype(str)) for st,g in mm.groupby("Study")}
    for st in sorted(mm.Study.unique()):
        for mode,mask in [("within_Study",mm.Study.eq(st)),("LOSO",~mm.Study.eq(st))]:
            m=mm[mask];r=models.fit_matrix(m,y.loc[m.index],"M1")
            r["BH_q"]=models.bh(r.raw_p);r["branch"]=mode;r["Study_removed_or_retained"]=st
            r["fixed_family_size"]=y.shape[1];fits.append(r)
    f=pd.concat(fits,ignore_index=True) if fits else pd.DataFrame(columns=["taxon_key","branch","raw_p","beta","Study_removed_or_retained"])
    save(f,out/"all_fixed_family_source_fits.tsv.gz");rec=[]
    for _,r in p.iterrows():
        d=f[f.taxon_key.eq(r.taxon_key)];w=d[d.branch.eq("within_Study")&d.raw_p.notna()]
        l=d[d.branch.eq("LOSO")&d.raw_p.notna()];sgn=np.sign(r.beta)
        wa=bool(len(w) and np.isfinite(sgn) and sgn!=0 and np.sign(w.beta).eq(sgn).all())
        la=bool(len(l) and np.isfinite(sgn) and sgn!=0 and np.sign(l.beta).eq(sgn).all())
        n=components(w.Study_removed_or_retained,sets)
        robust=bool(np.isfinite(r.raw_p) and n>=2 and len(l)>=2 and wa and la)
        reason=("source_robust_directional" if robust else "primary_not_estimable" if not np.isfinite(r.raw_p) else
                "fewer_than_two_independent_estimable_two_arm_Studies" if n<2 else
                "within_Study_direction_disagreement" if not wa else
                "insufficient_valid_LOSO" if len(l)<2 else "LOSO_direction_disagreement")
        rec.append(dict(taxon_key=r.taxon_key,primary_beta=r.beta,primary_raw_p=r.raw_p,primary_BH_q=r.BH_q,
                        ordinary_primary_hit=bool(r.BH_q<.2),n_estimable_two_arm_Studies=len(w),
                        n_independent_study_components=n,n_valid_LOSO=len(l),all_within_Study_directions_agree=wa,
                        all_valid_LOSO_directions_agree=la,source_robust=robust,classification=reason,
                        primary_hit_removed=False,comparison=c,body=b,scope=s,module=kind,panel=pid))
    d=pd.DataFrame(rec);save(d,out/"robustness_classification.tsv")
    dump(dict(status="PASS",config=config,fixed_family=len(p),source_fit_rows=len(f),primary_hits_removed=0,
              required_each_source_significant=False,
              independence="connected_components_of_shared_frozen_working_units"),out/"receipt.json")
    print("SOURCE_COMPLETE",c,s,pid,len(mm),len(p),flush=True)
    return d
