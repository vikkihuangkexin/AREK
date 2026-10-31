"""Historical metric definitions with strict inputs, exported distances/eigenvalues."""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.spatial.distance import pdist,squareform
from workflow.common import frame,tsv,dump
from scripts.diversity import historical_metrics as h
from scripts.diversity.historical_gamma import gamma_metric,gamma_pairwise

def compute(matrix,out,unit,cfg,metadata=None):
    out=Path(out);x=frame(matrix);zero=x.sum()==0;excluded=list(x.columns[zero]);x=x.loc[:,~zero]
    if not len(x.columns):raise ValueError('no positive samples remain')
    a=h.alpha_metrics(x)
    if unit!='counts' or not np.equal(x.values,np.floor(x.values)).all():a['Chao1']=np.nan
    tsv(a,out/'alpha.tsv')
    for name,dist in [('BrayCurtis',h.bray_curtis_matrix(x)),('Jaccard',h.jaccard_matrix(x)),('Aitchison',h.aitchison_matrix(x))]:
        if not np.isfinite(dist.values).all():raise ValueError(f'{name}: invalid distance')
        tsv(dist,out/(name+'.distance.tsv'))
        n=len(dist);j=np.eye(n)-np.ones((n,n))/n;b=-.5*j@(dist.values**2)@j
        vals,vec=np.linalg.eigh(b);order=np.argsort(vals)[::-1];vals=vals[order];vec=vec[:,order]
        positive=np.maximum(vals,0);total=positive.sum();coords=vec*np.sqrt(positive)
        co=pd.DataFrame(coords,index=dist.index,columns=[f'PCoA{i+1}' for i in range(n)])
        if n<2:co['PCoA2']=0
        tsv(co,out/(name+'.coordinates.tsv'))
        tsv(pd.DataFrame({'eigenvalue':vals,'positive_fraction':positive/total if total>0 else np.zeros(n)}),out/(name+'.eigenvalues.tsv'),False)
    groups=None;meta=pd.DataFrame(index=x.columns)
    if metadata:
        meta=pd.read_csv(metadata,sep='\t',index_col=0,dtype=str,keep_default_na=False).reindex(x.columns).replace('',np.nan)
        if cfg['group_column'] in meta:groups=meta[cfg['group_column']]
    rows=[{'group':'All_samples','n_samples':x.shape[1],**gamma_metric(x)}]
    if groups is not None:
        for label in groups.dropna().unique():
            ids=groups.index[groups==label];rows.append({'group':label,'n_samples':len(ids),**gamma_metric(x[ids])})
    tsv(pd.DataFrame(rows),out/'gamma.tsv',False)
    tests=[];reason='group_missing_or_single_arm'
    if cfg['design']=='descriptive':reason='descriptive_design_no_group_inference'
    if groups is not None and groups.nunique()>=2 and cfg['design']!='descriptive':
        repeated=cfg['subject_column'] in meta and meta[cfg['subject_column']].dropna().duplicated().any()
        if repeated or cfg['covariates'] or cfg['block_column']:reason='gamma permutation NA: repeated/block/adjusted design requires study-specific pooled estimand'
        else:
            reason=''
            for a,b in itertools.combinations(sorted(groups.dropna().unique()),2):
                gt=gamma_pairwise(x,groups,a,b,cfg['permutations'],cfg['seed'])
                if len(gt):tests.append(gt)
    gt=pd.concat(tests,ignore_index=True) if tests else pd.DataFrame([{'metric':'all','p_value':np.nan,'reason':reason or 'insufficient_samples'}])
    if 'reason' not in gt:gt['reason']=''
    if tests:
        from scipy.stats import false_discovery_control
        gt['FDR_p_value']=false_discovery_control(gt.p_value.values)
    tsv(gt,out/'gamma_tests.tsv',False)
    dump(out/'definitions.json',{'input_matrix':str(matrix),'unit':unit,'n_samples':x.shape[1],'sample_order':list(x.columns),'zero_samples_excluded':excluded,'Observed':'number of strictly positive taxa after configured filtering','Shannon':'-sum(p*ln(p)); natural log','Simpson':'1-sum(p^2)','Chao1':'S+f1*(f1-1)/(2*(f2+1)); only integer counts','rarefaction':'none; no depth/iterations','BrayCurtis':'sum(abs(x-y))/sum(x+y), uses selected matrix directly','Jaccard':'1-intersection/union of positive taxa','Aitchison':'Euclidean CLR(x+0.5), natural log; same historical pseudocount','PCoA':'full eigendecomposition; negative eigenvalues exported; positive eigenvalue sum for fractions; no silent distance correction','Gamma':'sum selected sample columns within group, then Observed/Shannon/1-D; counts weighted by depth, filtered_tss equally weighted samples','Gamma_test':reason or 'historical two-sided absolute pooled metric difference, independent label permutations, add-one p, BH correction'})
