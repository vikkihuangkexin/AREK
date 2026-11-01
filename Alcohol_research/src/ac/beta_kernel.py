"""General Freedman-Lane kernel implementation plus guarded commuting-block reduction."""
import argparse, subprocess
from .stats import *
from scipy.spatial.distance import pdist,squareform

def kernel(d):
    a=d*d;return -.5*(a-a.mean(0)[None,:]-a.mean(1)[:,None]+a.mean())

def permutations(m,B,seed,study_only=False):
    fields=["Study"] if study_only else ["Study","analysis_set","assay_analysis","material_analysis"]
    key=m[fields].astype(str).agg("::".join,axis=1)
    blocks=[np.flatnonzero(key.to_numpy()==v) for v in sorted(key.unique())]
    rng=np.random.default_rng(seed);arr=np.tile(np.arange(len(m),dtype=np.int32),(B,1))
    for i in range(B):
        for block in blocks:arr[i,block]=rng.permutation(block)
    return arr,blocks,key

def test_kernel(G,N,a,perms,blocks,force_general=False):
    N,_=independent(N);Q=np.linalg.qr(N,mode="reduced")[0]
    q=a-Q@(Q.T@a);length=np.linalg.norm(q)
    if length<1e-9:raise ValueError("nonidentifiable")
    q/=length;df=len(a)-N.shape[1]-1
    total=float(np.trace(G));null_resid=float(total-np.sum(Q*(G@Q)))
    ss=float(q@G@q);obs=ss/(null_resid-ss)*df
    commutes=all(np.allclose(N[b],N[b[0]],atol=1e-12) for b in blocks)
    if commutes and not force_general:
        ff=[]
        for off in range(0,len(perms),96):
            v=q[perms[off:off+96]].T;sss=np.sum(v*(G@v),0)
            ff.extend((sss/(null_resid-sss)*df).tolist())
        method="validated_commuting_full_nuisance_block_reduction"
    else:
        # Match installed vegan partial-dbRDA: permute the reduced residual kernel
        # AND its condition space, then re-residualize the tested constraint.
        # A fixed original q is wrong when nuisance varies within exchange blocks.
        qg=Q.T@G
        gr=G-Q@qg-(Q@qg).T+Q@(qg@Q)@Q.T
        ff=[]
        for perm in perms:
            qp=Q[perm]
            ar=a-qp@(qp.T@a);ar/=np.linalg.norm(ar)
            star=gr[np.ix_(perm,perm)]
            sr=float(ar@star@ar)
            nr=float(np.trace(gr))
            ff.append(sr/(nr-sr)*df)
        method="vegan_equivalent_general_reduced_residual_condition_space_permutation"
    ff=np.array(ff)
    return dict(F=obs,R2=ss/total,P=(int((ff>=obs-1e-12).sum())+1)/(len(ff)+1),
                residual_df=df,total_inertia=total,residual_inertia=null_resid,Alcohol_inertia=ss,
                B_actual=len(ff),method=method,commuting_guard_passed=commutes),ff

