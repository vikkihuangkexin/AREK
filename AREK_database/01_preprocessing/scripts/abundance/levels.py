"""Lineage-preserving aggregation, explicit denominators and audited filtering."""
import json,re
import numpy as np
import pandas as pd
from workflow.common import frame,taxonomy,tsv,dump
RANKS={'Domain':'d__','Kingdom':'k__','Phylum':'p__','Class':'c__','Order':'o__','Family':'f__','Genus':'g__','Species':'s__'}
ORDER=list(RANKS)
UNKNOWN={'','uncultured','unidentified','metagenome','unknown','na','unassigned'}

def lineage_at(taxon,level):
    parts=[p.strip() for p in str(taxon).split(';') if p.strip()]
    valid={}
    for p in parts:
        for name,prefix in RANKS.items():
            if p.startswith(prefix):valid[name]=p[len(prefix):].strip()
    above=ORDER[:ORDER.index(level)+1]
    lineage=';'.join(RANKS[k]+valid[k] for k in above if k in valid)
    named=valid.get(level,'').lower() not in UNKNOWN
    # Empty rank stays empty; status is separate from the biological taxonomy label.
    key=json.dumps([level,lineage,'assigned' if named else 'unclassified'],ensure_ascii=False,separators=(',',':'))
    return key,lineage,valid.get(level,''),'assigned' if named else 'unclassified'

def aggregate(table,taxfile,level,out,cfg,unit='counts',native_level=None):
    from pathlib import Path
    out=Path(out);(out/'not_applicable.json').unlink(missing_ok=True);x=frame(table);tax=taxonomy(taxfile) if taxfile else pd.Series(x.index,index=x.index)
    extra=tax.index.difference(x.index);missing=x.index.difference(tax.index)
    tax=tax.reindex(x.index).fillna('')
    if native_level:
        if level!=native_level:raise ValueError('native matrix only supports its declared level')
        # Preserve source keys as identity; unknown lineage is never invented.
        keys=list(x.index);details=[(k,str(tax.loc[k]),str(k),'source_native_key') for k in keys]
    else:
        details=[lineage_at(tax.loc[k],level) for k in x.index];keys=[d[0] for d in details]
        if not any(d[3]=='assigned' for d in details):
            dump(out/'not_applicable.json',{'level':level,'reason':'no named rank in source taxonomy','n_features':len(x)});return False
    grouped=x.groupby(keys,sort=False).sum();grouped.index.name='TaxonKey'
    before=x.sum();after=grouped.sum()
    if not np.allclose(before,after,atol=1e-8):raise ValueError('aggregation mass not conserved')
    if unit=='relative':
        # Input declared as 0-1 TSS. Percent or incomplete denominator requires explicit conversion upstream.
        if not np.allclose(before[before>0],1,atol=cfg['relative_tolerance']):raise ValueError('relative input does not sum to 1; declare/convert the true denominator; no automatic normalization')
        rel=grouped.copy()
    else:rel=grouped.div(before.replace(0,np.nan),axis=1).fillna(0)
    nonzero=before[before>0].index
    prevalence=(grouped[nonzero]>0).mean(axis=1) if len(nonzero) else pd.Series(0,index=grouped.index)
    mean=rel[nonzero].mean(axis=1) if len(nonzero) else pd.Series(0,index=grouped.index)
    keep=(prevalence>=cfg['min_prevalence'])&(mean>=cfg['min_mean_relative'])
    kept=grouped.loc[keep,nonzero];retained=kept.sum()
    filtered_tss=kept.div(retained.replace(0,np.nan),axis=1).fillna(0)
    tsv(grouped,out/('counts.tsv' if unit=='counts' else 'relative_input.tsv'))
    tsv(rel,out/'relative.tsv');tsv(rel.loc[keep],out/'filtered_original_denominator.tsv')
    tsv(kept,out/('filtered_counts.tsv' if unit=='counts' else 'filtered_relative_input.tsv'));tsv(filtered_tss,out/'filtered_tss.tsv')
    source_map=pd.DataFrame({'FeatureID':x.index,'TaxonKey':keys,'source_lineage':tax.values,'taxonomy_missing':x.index.isin(missing)})
    tsv(source_map,out/'feature_to_taxon.tsv',False)
    labels=pd.DataFrame(details,columns=['TaxonKey','lineage','rank_label','classification_status']).drop_duplicates('TaxonKey')
    tsv(labels,out/'taxonomy.tsv',False)
    tsv(pd.DataFrame({'prevalence':prevalence,'mean_relative':mean,'keep':keep}),out/'filter_audit.tsv')
    tsv(pd.DataFrame({'input_total':before,'aggregated_total':after,'retained_total':retained.reindex(before.index).fillna(0),'relative_sum':rel.sum(),'filtered_original_denominator_sum':rel.loc[keep].sum(),'filtered_tss_sum':filtered_tss.sum().reindex(before.index).fillna(0),'zero_library':before==0}),out/'mass_balance.tsv')
    dump(out/'contract.json',{'orientation':'taxon x sample','level':level,'input_unit':unit,'relative_unit':'0-1','relative_denominator':'all input features including unclassified','filtered_tss_denominator':'retained taxa after prevalence/mean filtering','input_samples':list(x.columns),'diversity_samples':list(filtered_tss.columns[filtered_tss.sum()>0]),'missing_taxonomy_features':list(missing),'extra_taxonomy_features':list(extra),'zero_libraries_excluded':list(before[before==0].index),'zero_retained_excluded':list(retained[retained==0].index),'structural_missing':'error; never zero-filled','unclassified':'preserve actual partial lineage, separate classification_status; no fake rank names','n_taxa':len(grouped),'n_filtered_taxa':int(keep.sum())})
    if not len(kept) or not (retained>0).any():raise ValueError(f'{level}: no nonzero retained taxa for diversity')
    return True

def verify(out):
    from pathlib import Path
    out=Path(out)
    if (out/'not_applicable.json').exists():return
    c=json.loads((out/'contract.json').read_text());mass=pd.read_csv(out/'mass_balance.tsv',sep='\t',index_col=0)
    if not np.allclose(mass.input_total,mass.aggregated_total):raise ValueError('mass mismatch')
    x=frame(out/'relative.tsv');s=x.sum()
    if not np.allclose(s[s>0],1,atol=1e-6):raise ValueError('TSS sum mismatch')
    f=frame(out/'filtered_tss.tsv');s=f.sum()
    if not np.allclose(s[s>0],1,atol=1e-6):raise ValueError('filtered TSS sum mismatch')
    if list(x.columns)!=c['input_samples']:raise ValueError('sample order changed')
