#!/usr/bin/env python3
"""Real-matrix contract transitions and stage-range regression checks in isolated output."""
from pathlib import Path
import sys,json,subprocess,shutil,time
import pandas as pd
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P));sys.dont_write_bytecode=True
from workflow.common import sha,dump
w=Path(sys.argv[1]).resolve();w.mkdir(parents=True,exist_ok=True);ev=w/'evidence';ev.mkdir(exist_ok=True)
tax=w/'taxonomy.tsv';original=P/'examples/matrix/taxonomy.tsv';shutil.copyfile(original,tax)
cfg={'route':'matrix','matrix':str(P/'examples/matrix/counts.tsv'),'taxonomy':str(tax),'levels':['Genus'],'memory_gb':1,'diversity':{'design':'descriptive','permutations':19}};conf=w/'config.json';dump(conf,cfg);result=[]
def run(name,command='resume',*extra):
 log=ev/(name+'.log');args=[sys.executable,str(P/'run_pipeline.py'),command,'--config',str(conf),'--dataset','SWITCH','--out',str(w/'outputs'),*extra]
 with log.open('w') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT)
 result.append({'test':name,'argv':args,'exit_code':p.returncode,'log_sha256':sha(log)});dump(ev/'results.json',result);assert p.returncode==0,log.read_text();return log.read_text()
run('initial','run');out=w/'outputs/SWITCH';assert (out/'abundance/Genus/counts.tsv').exists()
original_feature_state=sha(out/'state/features.json')
x=pd.read_csv(tax,sep='\t',keep_default_na=False);x.Taxon='d__Bacteria;p__;g__';x.to_csv(tax,sep='\t',index=False)
run('rank_becomes_unavailable')
assert (out/'abundance/Genus/not_applicable.json').exists();assert not (out/'abundance/Genus/counts.tsv').exists();assert not (out/'diversity/Genus').exists();assert list((out/'history').glob('*/Genus'))
run('validate_unavailable','validate')
shutil.copyfile(original,tax);run('rank_restored');assert not (out/'abundance/Genus/not_applicable.json').exists();assert (out/'diversity/Genus/alpha.tsv').exists()
run('validate_restored','validate')
state=sha(out/'state/features.json');cfg['abundance']={'min_prevalence':.5};dump(conf,cfg)
text=run('filter_configuration_change');assert 'REUSE SWITCH features' in text and 'RUN SWITCH levels.Genus' in text;assert sha(out/'state/features.json')==state
run('stage_range_existing_upstream','validate','--from-stage','levels','--to-stage','diversity')
# Explicit errors for missing/noninteger data remain visible; no constant-data p-value silently NA.
# The main full harness covers format faults and statistical guards; this file focuses on output contracts.
dump(ev/'COMPLETE.json',{'status':'passed','n_checks':len(result),'work_root':str(w)})
print('PASS',len(result),'contract checks')
