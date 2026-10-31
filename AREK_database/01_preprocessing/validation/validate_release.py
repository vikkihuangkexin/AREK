#!/usr/bin/env python3
"""Reproducible release verification. Real biological FASTQ chain and labeled software fault tests."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,shutil,subprocess,sys,time
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(P))
from workflow.common import sha,dump,frame

def main():
    a=argparse.ArgumentParser();a.add_argument('--work-root',required=True);a.add_argument('--classifier');a.add_argument('--matrix-only',action='store_true');a.add_argument('--evidence-dir');args=a.parse_args()
    work=Path(args.work_root).resolve();work.mkdir(parents=True,exist_ok=True);ev=Path(args.evidence_dir).resolve() if args.evidence_dir else work/'evidence';ev.mkdir(parents=True,exist_ok=True);results=[]
    def invoke(name,command,expected=0):
        begin=time.time();log=ev/(name+'.log')
        with log.open('w') as f:proc=subprocess.run(list(map(str,command)),cwd=P,stdout=f,stderr=subprocess.STDOUT)
        result={'test':name,'argv':list(map(str,command)),'exit_code':proc.returncode,'expected_exit_code':expected,'elapsed_seconds':time.time()-begin,'log':str(log),'log_sha256':sha(log)};results.append(result);dump(ev/'verification_results.json',results)
        if proc.returncode!=expected:raise AssertionError(f'{name}: exit {proc.returncode}, expected {expected}; {log}')
        return log.read_text()
    def cli(name,cmd,*rest,expected=0):return invoke(name,[sys.executable,P/'run_pipeline.py',cmd,*rest],expected)
    cli('help','--help')
    for cmd in ['scan','check','dry-run','run','resume','validate','redraw']:cli('help_'+cmd,cmd,'--help')
    for f in (P/'config').glob('*.json'):
        from workflow.cli import load_config
        load_config(f)
    results.append({'test':'all_example_configs_parse','passed':True})
    cli('scan_real','scan','--root',P/'examples/real_fastq','--out',work/'scan_real')
    manifest=pd.read_csv(work/'scan_real/samples.tsv',sep='\t');assert len(manifest)==4 and (manifest.layout=='PE').all()
    # Fault inputs are software tests only; no synthetic biological validation claim.
    faults=work/'faults';faults.mkdir(exist_ok=True)
    for name in ['X_L001_R1.fastq','X_L001_R2.fastq','X_L002_R1.fastq','X_L002_R2.fastq','ORPHAN_R1.fastq']:(faults/name).write_text('@a\nACGT\n+\nIIII\n')
    cli('scan_conflict','scan','--root',faults,'--out',work/'scan_conflicts',expected=2)
    conflicts=pd.read_csv(work/'scan_conflicts/conflicts.tsv',sep='\t');assert conflicts.error.str.contains('multiple lanes').any() and conflicts.error.str.contains('orphan').any()
    matrixcfg=P/'config/matrix_example.json';matrixargs=['--config',matrixcfg,'--dataset','REAL_MATRIX','--out',work/'matrix']
    cli('matrix_check','check',*matrixargs)
    planned=cli('matrix_dry_run','dry-run',*matrixargs);assert not (work/'matrix').exists();assert 'internal Python' in planned
    cli('matrix_run','run',*matrixargs);cli('matrix_validate','validate',*matrixargs)
    cli('redraw_real_matrix','redraw','--folder',work/'matrix/REAL_MATRIX/diversity/Genus','--dataset','REAL_MATRIX','--level','Genus')
    # Independent elementary checks, not calls back into the implementation under test.
    numeric=[]
    for folder in (work/'matrix/REAL_MATRIX/diversity').iterdir():
        if not folder.is_dir():continue
        x=frame(work/'matrix/REAL_MATRIX/abundance'/folder.name/'filtered_counts.tsv');alpha=pd.read_csv(folder/'alpha.tsv',sep='\t',index_col=0)
        for sid in alpha.index:
            v=x[sid].values;pp=v[v>0]/sum(v);f1=sum(v==1);f2=sum(v==2)
            ref=[sum(v>0),-sum(float(p)*__import__('math').log(float(p)) for p in pp),1-sum(float(p)**2 for p in pp),sum(v>0)+f1*(f1-1)/(2*(f2+1))]
            assert np.allclose(alpha.loc[sid,['Observed','Shannon','Simpson','Chao1']],ref,atol=1e-12)
        dist=pd.read_csv(folder/'BrayCurtis.distance.tsv',sep='\t',index_col=0)
        a,b=x.iloc[:,0].values,x.iloc[:,1].values;ref=sum(abs(float(i)-float(j)) for i,j in zip(a,b))/sum(float(i)+float(j) for i,j in zip(a,b));assert abs(dist.iloc[0,1]-ref)<1e-12
        mass=pd.read_csv(work/'matrix/REAL_MATRIX/abundance'/folder.name/'mass_balance.tsv',sep='\t');assert np.allclose(mass.input_total,mass.aggregated_total)
        numeric.append({'level':folder.name,'n_samples':len(alpha),'n_features_filtered':len(x),'alpha_max_abs_tolerance':1e-12,'BrayCurtis_01':ref,'mass_conserved':True})
    dump(ev/'independent_numeric_checks.json',numeric)
    if not args.matrix_only:
        assert args.classifier,'--classifier required for full biological validation';classifier=Path(args.classifier).resolve();assert classifier.is_file()
        cfg=P/'config/real_pe_example.json';base=['--config',cfg,'--dataset','PRJDB10075','--classifier',classifier,'--out',work/'real_pe']
        cli('pe_check','check',*base);cli('pe_dry_run','dry-run',*base);assert not (work/'real_pe').exists()
        start=time.time();cli('pe_first_run','run',*base);elapsed=time.time()-start
        ds=work/'real_pe/PRJDB10075';states={p.name:sha(p) for p in (ds/'state').glob('*.json')}
        resumed=cli('pe_resume_unchanged','resume',*base);assert '\nRUN PRJDB10075' not in resumed
        assert states=={p.name:sha(p) for p in (ds/'state').glob('*.json')}
        broken=ds/'diversity/Genus/BrayCurtis.distance.tsv';broken.unlink()
        repaired=cli('pe_resume_missing_genus_distance','resume',*base)
        ran=[line for line in repaired.splitlines() if line.startswith('RUN ')];assert ran==['RUN PRJDB10075 diversity.Genus','RUN PRJDB10075 statistics.Genus','RUN PRJDB10075 plot.Genus'],ran
        cli('pe_validate','validate',*base)
        # Batch failure must not prevent the already valid independent real dataset from completing.
        bad=work/'broken.fastq';bad.write_text('@broken\nACGT\n+\nI\n')
        sample_rows=list(csv.DictReader((P/'examples/real_samples.tsv').open(),delimiter='\t'))
        for r in sample_rows:
            for k in ['r1','r2']:r[k]=str((P/'examples'/r[k]).resolve())
        badrow=dict(sample_rows[0]);badrow.update(dataset_id='CRA_FAULT_CASE',sample_id='BROKEN',run_id='BROKEN',r1=str(bad))
        sample_rows.insert(0,badrow);manifest=work/'batch_with_fault.tsv'
        with manifest.open('w') as f:w=csv.DictWriter(f,fieldnames=list(sample_rows[0]),delimiter='\t');w.writeheader();w.writerows(sample_rows)
        batchcfg=json.loads(cfg.read_text());batchcfg['samples']=str(manifest);batchcfg['classifier']=str(classifier);f=work/'batch.json';dump(f,batchcfg)
        log=cli('partial_failure_continues','resume','--config',f,'--out',work/'real_pe',expected=2);assert 'invalid FASTQ record 1' in log and 'REUSE PRJDB10075 import' in log
        # Restore canonical config downstream fingerprints after manifest-change test.
        cli('pe_restore_original_manifest','resume',*base);cli('pe_final_validate','validate',*base)
        # Single-end uses the explicitly selected biological R1, separately from the PE dataset.
        cli('se_real_chain','run','--config',P/'config/real_se_example.json','--dataset','PRJDB10075','--classifier',classifier,'--out',work/'real_se')
        cli('se_validate','validate','--config',P/'config/real_se_example.json','--dataset','PRJDB10075','--classifier',classifier,'--out',work/'real_se')
        dump(ev/'real_chain_measurements.json',{'dataset':'PRJDB10075','runs':4,'reads_per_run':2000,'paired_reads':True,'pe_elapsed_seconds':elapsed,'threads':2,'memory_requested_gb':8,'note':'memory request is a preflight estimate, not measured peak RSS; this fixture has no invented group labels','upstream_features':len(frame(ds/'upstream/features/feature-table.tsv')),'supported_levels':[p.name for p in (ds/'diversity').iterdir() if p.is_dir()],'external_classifier_sha256':sha(classifier)})
        # Archive real raw-stage results/logs/checkpoints, excluding optional scratch caches.
        copy=ev/'real_pe_results';shutil.copytree(ds,copy,ignore=shutil.ignore_patterns('tmp','pipeline.lock'),dirs_exist_ok=True)
        shutil.copytree(work/'real_se/PRJDB10075',ev/'real_se_results',ignore=shutil.ignore_patterns('tmp','pipeline.lock'),dirs_exist_ok=True)
    # Explicit synthetic metadata/statistical design checks on the real matrix, not biological inference.
    matrix=Path(P/'examples/matrix/counts.tsv');ids=list(frame(matrix).columns)
    for mode in ['independent','within_subject','invariant_subject','unknown_groups']:
        meta=pd.DataFrame({'sample_id':ids,'group':['A','A','B','B'],'subject_id':['S1','S2','S3','S4']})
        config=json.loads(matrixcfg.read_text());config.update(matrix=str(matrix),taxonomy=str(P/'examples/matrix/taxonomy.tsv'),levels=['Genus'])
        if mode in ['within_subject','invariant_subject']:
            meta.subject_id=['S1','S2','S1','S2'] if mode=='within_subject' else ['S1','S1','S2','S2']
        if mode=='unknown_groups':meta.group=['Unknown']*4
        m=work/(mode+'.metadata.tsv');meta.to_csv(m,sep='\t',index=False);config['metadata']=str(m);config['diversity'].update(group_column='group',subject_column='subject_id',design='within_subject' if mode in ['within_subject','invariant_subject'] else 'independent',block_column='subject_id' if mode in ['within_subject','invariant_subject'] else '',permutations=19)
        f=work/(mode+'.json');dump(f,config);cli('stat_'+mode,'run','--config',f,'--dataset',mode,'--out',work/'stats')
        per=pd.read_csv(work/'stats'/mode/'diversity/Genus/permanova.tsv',sep='\t')
        if mode=='invariant_subject':assert per.reason.str.contains('group_invariant_within_subject').all()
        if mode=='unknown_groups':assert per.p_value.isna().all()
        if mode=='independent':assert per.p_value.notna().all(),per.to_string()
    # Relative-only must never emit a counts file.
    rel=frame(matrix);rel=rel/rel.sum();f=work/'relative.tsv';rel.to_csv(f,sep='\t');cfg=json.loads(matrixcfg.read_text());cfg.update(matrix=str(f),taxonomy=str(P/'examples/matrix/taxonomy.tsv'),matrix_unit='relative',levels=['Genus']);cfg['diversity']['input']='filtered_tss';conf=work/'relative.json';dump(conf,cfg)
    cli('relative_only','run','--config',conf,'--dataset','relative_only','--out',work/'relative')
    assert not list((work/'relative/relative_only').rglob('*counts.tsv'))
    alpha=pd.read_csv(work/'relative/relative_only/diversity/Genus/alpha.tsv',sep='\t');assert alpha.Chao1.isna().all()
    # Specific missing DB error, without starting a compute job.
    cfg=json.loads((P/'config/real_pe_example.json').read_text());cfg['samples']=str(P/'examples/real_samples.tsv');cfg['classifier']=str(work/'missing-classifier.qza');conf=work/'missing_db.json';dump(conf,cfg)
    err=cli('missing_database','check','--config',conf,'--dataset','PRJDB10075',expected=2);assert 'classifier QZA missing' in err
    dump(ev/'verification_results.json',results);dump(ev/'VERIFICATION_COMPLETE.json',{'status':'passed','matrix_only':args.matrix_only,'n_command_tests':len(results),'finished_epoch':time.time(),'package':str(P),'work_root':str(work)})
    print('PASS',len(results),'checks;',ev)
if __name__=='__main__':main()
