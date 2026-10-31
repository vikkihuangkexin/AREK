#!/usr/bin/env python3
"""Labeled software fixtures for covariate adjustment, restricted permutations, descriptive-only Gamma."""
from pathlib import Path
import json,sys,subprocess
import pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P));sys.dont_write_bytecode=True
from workflow.common import frame,dump,sha
w=Path(sys.argv[1]).resolve();w.mkdir(exist_ok=True,parents=True);ev=w/'evidence';ev.mkdir(exist_ok=True);rows=[]
ids=list(frame(P/'examples/matrix/counts.tsv').columns)
for name in ['covariate','restricted_block','descriptive_group']:
    meta=pd.DataFrame({'sample_id':ids,'group':['A','A','B','B'],'subject_id':['S1','S2','S3','S4'],'age':[31,47,36,58],'batch':['B1','B2','B1','B2']});m=w/(name+'.metadata.tsv');meta.to_csv(m,sep='\t',index=False)
    cfg={'route':'matrix','matrix':str(P/'examples/matrix/counts.tsv'),'taxonomy':str(P/'examples/matrix/taxonomy.tsv'),'levels':['Genus'],'memory_gb':1,'metadata':str(m),'diversity':{'group_column':'group','subject_column':'subject_id','design':'descriptive' if name=='descriptive_group' else 'independent','covariates':['age'] if name=='covariate' else [],'block_column':'batch' if name=='restricted_block' else '', 'permutations':19}}
    conf=w/(name+'.json');dump(conf,cfg);log=ev/(name+'.log');argv=[sys.executable,str(P/'run_pipeline.py'),'run','--config',str(conf),'--dataset',name,'--out',str(w/'outputs')]
    with log.open('w') as f:r=subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT)
    assert not r.returncode,log.read_text();d=w/'outputs'/name/'diversity/Genus';per=pd.read_csv(d/'permanova.tsv',sep='\t');gamma=pd.read_csv(d/'gamma_tests.tsv',sep='\t')
    if name=='descriptive_group':assert per.p_value.isna().all() and gamma.p_value.isna().all() and gamma.reason.str.contains('descriptive').all()
    else:assert per.p_value.notna().all(),per.to_string()
    rows.append({'fixture':name,'synthetic_metadata':True,'real_count_matrix':True,'exit_code':r.returncode,'argv':argv,'log_sha256':sha(log),'permanova':per.to_dict('records'),'gamma_test_reason':gamma.reason.tolist()})
    dump(ev/'results.json',rows)
    for file in ['alpha_tests.tsv','permanova.tsv','dispersion.tsv','gamma_tests.tsv','metadata.tsv','statistics_session.txt']:
        (ev/(name+'_'+file)).write_bytes((d/file).read_bytes())
dump(ev/'COMPLETE.json',{'status':'passed','n_checks':len(rows)});print('PASS',len(rows),'additional statistical designs')
