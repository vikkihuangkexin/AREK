"""Redraw directly from real metric TSVs. PDF + 300dpi PNG, captions carry n."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from workflow.common import tsv

def plot(folder,metadata=None,group_column='',dataset='',level=''):
    folder=Path(folder);a=pd.read_csv(folder/'alpha.tsv',sep='\t',index_col=0)
    meta=pd.read_csv(metadata,sep='\t',index_col=0,dtype=str,keep_default_na=False) if metadata else pd.DataFrame(index=a.index)
    groups=meta[group_column].reindex(a.index).replace('','Unknown') if group_column in meta else pd.Series('Unknown',index=a.index)
    groups=groups.fillna('Unknown');captions=[]
    def save(fig,stem,caption):
        for fmt in ['pdf','png']:
            p=folder/'figures'/fmt/(stem+'.'+fmt);p.parent.mkdir(parents=True,exist_ok=True);fig.savefig(p,bbox_inches='tight',dpi=300)
        plt.close(fig);captions.append({'figure':stem,'caption':caption})
    for metric in a:
        fig,ax=plt.subplots(figsize=(6,4.8));levels=list(groups.unique());arrays=[a.loc[groups==g,metric].dropna().values for g in levels]
        usable=[(g,v) for g,v in zip(levels,arrays) if len(v)]
        if usable:ax.boxplot([v for g,v in usable],labels=[f'{g}\nn={len(v)}' for g,v in usable])
        ax.set(title=f'{dataset} {level}: {metric}',ylabel=metric)
        save(fig,'alpha_'+metric,f'{dataset} {level}; {metric}; '+ '; '.join(f'{g}: n={len(v)}' for g,v in usable)+'; see alpha.tsv; Unknown kept descriptive')
    for name in ['BrayCurtis','Jaccard','Aitchison']:
        co=pd.read_csv(folder/(name+'.coordinates.tsv'),sep='\t',index_col=0);ev=pd.read_csv(folder/(name+'.eigenvalues.tsv'),sep='\t');fig,ax=plt.subplots(figsize=(6,5))
        for g in groups.unique():
            ids=groups.index[groups==g].intersection(co.index);ax.scatter(co.loc[ids,'PCoA1'],co.loc[ids,'PCoA2'],label=f'{g} (n={len(ids)})')
        fr=list(ev.positive_fraction)+[0,0];ax.set(xlabel=f'PCoA1 ({100*fr[0]:.1f}%)',ylabel=f'PCoA2 ({100*fr[1]:.1f}%)',title=f'{dataset} {level}: {name}');ax.legend(fontsize=8)
        save(fig,'PCoA_'+name,f'{name}; n={len(co)}; coordinates and full eigenvalues in matching TSV; fraction uses positive eigenvalue sum')
    ga=pd.read_csv(folder/'gamma.tsv',sep='\t')
    for metric in ['Observed','Shannon','Simpson']:
        fig,ax=plt.subplots(figsize=(6,4.8));ax.bar(ga.group.astype(str),ga[metric]);ax.set(title=f'{dataset} {level}: gamma {metric}',ylabel=metric);ax.tick_params(axis='x',rotation=25)
        save(fig,'gamma_'+metric,'Pooled selected matrix; '+ '; '.join(f'{r.group}: n={r.n_samples}' for r in ga.itertuples())+'; see gamma.tsv')
    tsv(pd.DataFrame(captions),folder/'captions.tsv',False)
