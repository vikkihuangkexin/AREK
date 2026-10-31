from pathlib import Path
import argparse,csv,copy,fcntl,json,os,shutil,sys,traceback
import numpy as np,pandas as pd
from workflow.common import ROOT,Runner,dump,frame,taxonomy,tsv,digest
from workflow.scanner import scan,qc_pair
from workflow.routes import amplicon,shotgun
from scripts.abundance import levels as lev
from scripts.diversity.compute import compute
from scripts.diversity.plot import plot
STAGES=['preprocess','features','taxonomy','levels','diversity','plot']

def merge(base,extra):
    for k,v in extra.items():
        if k not in base:raise ValueError(f'unknown config key: {k}')
        if isinstance(base[k],dict):
            if not isinstance(v,dict):raise ValueError(f'config {k} must be object')
            merge(base[k],v)
        else:base[k]=v
    return base

def load_config(path):
    p=Path(path).resolve();cfg=merge(json.loads((ROOT/'config/historical.json').read_text()),json.loads(p.read_text()))
    for k in ['classifier','samples','matrix','taxonomy','metadata']:
        if cfg[k]:cfg[k]=str((p.parent/cfg[k]).resolve())
    for k in ['host_index','kraken_db']:
        if cfg['shotgun'][k]:cfg['shotgun'][k]=str((p.parent/cfg['shotgun'][k]).resolve())
    if cfg['route'] not in ['amplicon_vsearch97','shotgun_kraken_bracken','matrix']:raise ValueError('unknown route')
    if cfg['matrix_unit'] not in ['counts','relative']:raise ValueError('matrix_unit must be counts or relative (0-1)')
    if not 1<=int(cfg['threads'])<=os.cpu_count():raise ValueError('threads exceeds available CPU or <1')
    if cfg['memory_gb']<=0 or cfg['min_free_disk_gb']<0:raise ValueError('invalid resource request')
    if any(l not in lev.RANKS for l in cfg['levels']):raise ValueError('unknown taxonomy level')
    if len(cfg['levels'])!=len(set(cfg['levels'])):raise ValueError('duplicate levels')
    if cfg['native_level'] and cfg['native_level'] not in cfg['levels']:raise ValueError('native_level must be selected in levels')
    if cfg['diversity']['input'] not in ['filtered_counts','filtered_tss','counts','relative']:raise ValueError('invalid diversity input')
    if cfg['diversity']['design'] not in ['descriptive','independent','within_subject']:raise ValueError('unknown statistical design')
    if int(cfg['diversity']['permutations'])<1:raise ValueError('permutations must be positive')
    for k in ['min_prevalence','min_mean_relative']:
        if not 0<=cfg['abundance'][k]<=1:raise ValueError(k+' must lie in [0,1]')
    if cfg['shotgun']['variant'] not in ['original','v2']:raise ValueError('shotgun variant must be original or v2')
    return cfg

def read_samples(cfg):
    if not cfg['samples']:return {}
    p=Path(cfg['samples']);rows=list(csv.DictReader(p.open(),delimiter='\t'));groups={}
    required={'dataset_id','sample_id','run_id','layout','r1','r2'}
    if not rows or not required<=set(rows[0]):raise ValueError('sample manifest requires '+','.join(sorted(required)))
    for r in rows:
        for k in ['r1','r2']:
            if r[k]:r[k]=str((p.parent/r[k]).resolve())
        groups.setdefault(r['dataset_id'],[]).append(r)
    return groups

def check_dataset(cfg,samples,ds,deep=True):
    if not ds or ds in ['.','..'] or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in ds):raise ValueError('unsafe/empty dataset ID')
    result={'dataset_id':ds,'route':cfg['route'],'warnings':[]}
    if cfg['route']=='matrix':
        if not cfg['matrix']:raise ValueError('matrix path required')
        x=frame(cfg['matrix'])
        if cfg['matrix_unit']=='relative' and cfg['diversity']['input'] in ['counts','filtered_counts']:raise ValueError('relative-only input has no counts; choose diversity.input=relative or filtered_tss')
        if cfg['matrix_unit']=='counts' and not np.equal(x.values,np.floor(x.values)).all():result['warnings'].append('noninteger counts; Chao1 will be NA')
        if not cfg['taxonomy'] and not cfg['native_level']:raise ValueError('matrix input requires taxonomy or native_level; no invented taxonomy')
        if cfg['taxonomy']:
            t=taxonomy(cfg['taxonomy']);result['missing_taxonomy']=len(x.index.difference(t.index))
        ids=list(x.columns)
    else:
        if not samples:raise ValueError('no explicit manifest samples for dataset')
        ids=[s['sample_id'] for s in samples]
        for sid in ids:
            if not sid or sid in ['.','..'] or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in sid):raise ValueError(f'unsafe/empty sample ID: {sid}')
        if len(ids)!=len(set(ids)):raise ValueError('duplicate sample IDs: provide explicit lane/run merge upstream; no silent merging')
        runs=[s['run_id'] for s in samples]
        if len(runs)!=len(set(runs)):raise ValueError('multiple rows per run; lane merging requires explicit upstream decision')
        for key in ['specimen_id']:
            vals=[s.get(key,'') for s in samples if s.get(key,'')]
            if len(vals)!=len(set(vals)):raise ValueError('multiple runs for same specimen: explicitly resolve technical replicates before independent analysis')
        subjects=[s.get('subject_id','') for s in samples if s.get('subject_id','')]
        if len(subjects)!=len(set(subjects)) and cfg['diversity']['subject_column']!='subject_id':
            raise ValueError('manifest contains repeated subject_id; set subject_column=subject_id and an appropriate design; no independent inference')
        files=[]
        for s in samples:
            if s['layout']!=cfg['layout']:raise ValueError(f"{s['sample_id']}: layout differs from explicit dataset layout")
            if s['layout'] not in ['SE','PE']:raise ValueError('layout must be SE/PE')
            if not s['run_id']:raise ValueError('run_id required')
            if s['layout']=='PE' and not s['r2']:raise ValueError(f"{s['sample_id']}: PE missing R2")
            if s['layout']=='SE' and s['r2']:raise ValueError(f"{s['sample_id']}: SE contains R2; record selection explicitly")
            for k in ['r1']+(['r2'] if s['layout']=='PE' else []):
                if not Path(s[k]).is_file():raise ValueError(f"{s['sample_id']}: input missing {k}: {s[k]}")
                files.append(s[k])
            if deep:qc_pair(s['r1'],s['r2'] or None,1000)
        if len(files)!=len(set(files)):raise ValueError('same FASTQ assigned to multiple manifest records')
    if cfg['metadata']:
        meta=pd.read_csv(cfg['metadata'],sep='\t',index_col=0,dtype=str,keep_default_na=False)
        if meta.index.duplicated().any():raise ValueError('duplicate metadata sample IDs')
        if not set(ids)&set(meta.index):raise ValueError('metadata has zero ID overlap')
        result['metadata_missing_samples']=sorted(set(ids)-set(meta.index));result['metadata_extra_samples']=sorted(set(meta.index)-set(ids))
        for k in ['group_column','subject_column','block_column']:
            c=cfg['diversity'][k]
            if c and c not in meta:raise ValueError(f'metadata missing configured {k}: {c}')
        for c in cfg['diversity']['covariates']:
            if c not in meta:raise ValueError(f'metadata missing covariate {c}')
        # Known subject information must be configured; it must never be ignored.
        cols=[c for c in ['subject_id','participant_ID','participant_id','SubjectID'] if c in meta]
        for c in cols:
            v=meta[c].replace('',np.nan).dropna()
            if v.duplicated().any() and cfg['diversity']['subject_column']!=c:raise ValueError(f'known repeated subjects in {c}; configure subject_column and design')
    elif cfg['diversity']['group_column'] or cfg['diversity']['covariates']:raise ValueError('group/covariates configured but metadata missing')
    tools=['Rscript']
    if cfg['route']=='amplicon_vsearch97':
        tools+=['qiime','biom','vsearch']
        if not cfg['classifier'] or not Path(cfg['classifier']).is_file():raise ValueError('classifier QZA missing; see docs/databases_zh.md')
        from workflow.common import validate_file
        validate_file(cfg['classifier'])
    if cfg['route']=='shotgun_kraken_bracken':
        if cfg['layout']!='PE' and cfg['shotgun']['variant']!='original':raise ValueError('shotgun v2 requires PE; original branch supports SE')
        tools+=['fastp','kraken2','bracken','samtools'];p=cfg['shotgun']
        for n in ['hash.k2d','opts.k2d','taxo.k2d',f"database{p['read_length']}mers.kmer_distrib"]:
            if not p['kraken_db'] or not (Path(p['kraken_db'])/n).is_file():raise ValueError('missing Kraken/Bracken database file '+n)
        if not p['skip_host_removal']:
            tools+=['bowtie2']
            if not p['host_index'] or len(list(Path(p['host_index']).parent.glob(Path(p['host_index']).name+'.*.bt2*')))<6:raise ValueError('missing complete host Bowtie2 index')
        if p['variant']=='original' and not p['correction']:result['warnings'].append('original source enables fastp correction; requested configuration differs')
    for tool in tools:
        if not shutil.which(tool):raise ValueError(f'missing executable {tool}; activate route environment (see envs)')
    if deep:
        import subprocess
        p=subprocess.run(['Rscript','-e','stopifnot(requireNamespace("vegan",quietly=TRUE),requireNamespace("permute",quietly=TRUE))'],capture_output=True,text=True)
        if p.returncode:raise ValueError('R dependencies unavailable: '+p.stderr[-1000:])
    mem={line.split(':')[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines()}
    available=mem.get('MemAvailable',mem['MemFree'])/1024**2
    if available<cfg['memory_gb']:raise ValueError(f'requested memory {cfg["memory_gb"]} GiB > current available {available:.2f} GiB')
    result.update(n_samples=len(ids),available_memory_gb=available,tools={k:shutil.which(k) for k in tools})
    return result

def materialize_metadata(cfg,samples,ids,out):
    meta=pd.read_csv(cfg['metadata'],sep='\t',index_col=0,dtype=str,keep_default_na=False) if cfg['metadata'] else pd.DataFrame(index=ids)
    original_metadata_ids=set(meta.index) if cfg['metadata'] else set()
    # Include explicit run/specimen/subject mapping, never derive biological groups.
    for s in samples:
        for k in ['run_id','specimen_id','subject_id']:
            if s.get(k):
                old=meta.loc[s['sample_id'],k] if k in meta and s['sample_id'] in meta.index else ''
                if pd.notna(old) and old and old!=s[k]:raise ValueError(f'manifest/metadata {k} mismatch: {s["sample_id"]}')
                meta.loc[s['sample_id'],k]=s[k]
    gc=cfg['diversity']['group_column']
    if gc and gc in meta:
        meta['source_group_value']=meta[gc]
        meta[gc]=meta[gc].where(~meta[gc].isin(cfg['diversity']['missing_group_values']),'')
    audit=pd.DataFrame({'sample_id':list(dict.fromkeys([*ids,*meta.index]))});audit['in_matrix']=audit.sample_id.isin(ids);audit['in_metadata']=audit.sample_id.isin(original_metadata_ids);audit['has_manifest_mapping']=audit.sample_id.isin([s['sample_id'] for s in samples])
    tsv(audit,out/'sample_audit.tsv',False);tsv(meta.reindex(ids).fillna('').rename_axis('sample_id'),out/'metadata.tsv')


def mark_not_applicable(out,level,folder):
    import time
    previous=out/'diversity'/level
    if previous.exists():
        archive=out/'history'/('not_applicable_'+str(time.time_ns()))
        archive.mkdir(parents=True,exist_ok=True)
        shutil.move(str(previous),str(archive/level))
    for prefix in ['metadata.','diversity.','statistics.','plot.']:
        cp=out/'state'/(prefix+level+'.json')
        if cp.exists():
            record=json.loads(cp.read_text());record['status']='invalidated_not_applicable';record['reason']='source no longer contains this named taxonomy rank';dump(cp,record)
    dump(folder/'not_applicable.json',{'level':level,'reason':'no named rank in source taxonomy'})

def execute(cfg,samples,ds,out,mode,start,end):
    from workflow.common import sha
    runner=Runner(out,cfg,ds,mode);runner.start=start;runner.end=end
    if mode!='dry-run':
        out.mkdir(parents=True,exist_ok=True)
        if shutil.disk_usage(out).free/1024**3<cfg['min_free_disk_gb']:raise ValueError('insufficient free output disk')
        dump(out/'effective_config.json',cfg) if mode!='validate' else None
    if cfg['route']=='amplicon_vsearch97':sources=amplicon(runner,samples,cfg)
    elif cfg['route']=='shotgun_kraken_bracken':sources=shotgun(runner,samples,cfg)
    else:
        tb=out/'input/features.tsv';tx=out/'input/taxonomy.tsv'
        outputs=[tb]+([tx] if cfg['taxonomy'] else [])
        def copy_inputs():
            tsv(frame(cfg['matrix']),tb)
            if cfg['taxonomy']:shutil.copyfile(cfg['taxonomy'],tx)
        runner.step('features',[cfg['matrix']]+([cfg['taxonomy']] if cfg['taxonomy'] else []),outputs,{'unit':cfg['matrix_unit'],'native_level':cfg['native_level']},copy_inputs,lambda:frame(tb))
        sources=[(tb,tx if cfg['taxonomy'] else None,cfg['native_level'] or None,cfg['matrix_unit'],'features')]
    if STAGES.index(end)<STAGES.index('levels'):return runner
    for tb,tx,native,unit,dep in sources:
        selected=[native] if native else cfg['levels']
        for level in selected:
            if level not in cfg['levels']:continue
            d=out/'abundance'/level;tag='levels.'+level
            applicable=True
            if mode!='dry-run' and tx and not native:
                applicable=any(lev.lineage_at(s,level)[3]=='assigned' for s in taxonomy(tx).values)
            if not applicable:
                marker=d/'not_applicable.json';runner.step(tag,[tb,tx],[marker],cfg['abundance'],lambda level=level,d=d:mark_not_applicable(out,level,d),deps=[dep]);continue
            files=['relative.tsv','filtered_original_denominator.tsv','filtered_tss.tsv','feature_to_taxon.tsv','taxonomy.tsv','filter_audit.tsv','mass_balance.tsv','contract.json']+(['counts.tsv','filtered_counts.tsv'] if unit=='counts' else ['relative_input.tsv','filtered_relative_input.tsv'])
            runner.step(tag,[tb]+([tx] if tx else []),[d/f for f in files],{**cfg['abundance'],'unit':unit,'level':level,'native':native},lambda tb=tb,tx=tx,level=level,d=d,native=native,unit=unit:lev.aggregate(tb,tx,level,d,cfg['abundance'],unit,native),lambda d=d:lev.verify(d),[dep])
            dc=cfg['diversity'];choice=dc['input'];matrix=d/(choice+'.tsv')
            if unit=='relative' and choice in ['counts','filtered_counts']:raise ValueError('relative input cannot produce counts')
            div=out/'diversity'/level;meta=div/'metadata.tsv';audit=div/'sample_audit.tsv';tagm='metadata.'+level
            def metadata_action(matrix=matrix,div=div):materialize_metadata(cfg,samples,list(frame(matrix).columns),div)
            runner.step(tagm,[matrix]+([cfg['metadata']] if cfg['metadata'] else [])+([cfg['samples']] if cfg['samples'] else []),[meta,audit],{'metadata':cfg['metadata']},metadata_action,deps=['levels.'+level])
            outputs=[div/f for f in ['alpha.tsv','gamma.tsv','gamma_tests.tsv','definitions.json']]+[div/(name+'.'+suffix+'.tsv') for name in ['BrayCurtis','Jaccard','Aitchison'] for suffix in ['distance','coordinates','eigenvalues']]
            tagd='diversity.'+level
            selected_unit='counts' if choice in ['counts','filtered_counts'] else 'relative'
            runner.step(tagd,[matrix,meta],outputs,dc,lambda matrix=matrix,div=div,selected_unit=selected_unit,meta=meta:compute(matrix,div,selected_unit,dc,meta),lambda div=div:verify_diversity(div),[tagm])
            cmd=['Rscript',str(ROOT/'scripts/diversity/statistics.R'),str(div),str(meta),dc['group_column'],dc['subject_column'],dc['block_column'],dc['design'],','.join(dc['covariates']),str(dc['permutations']),str(dc['seed'])]
            statoutputs=[div/f for f in ['alpha_tests.tsv','permanova.tsv','dispersion.tsv','statistics_sample_audit.tsv','statistics_session.txt']]
            runner.step('statistics.'+level,[meta,*outputs],statoutputs,dc,lambda cmd=cmd:runner.cmd(cmd),deps=[tagd],plan=[cmd])
            stems=['alpha_'+m for m in ['Observed','Shannon','Simpson','Chao1']]+['PCoA_'+n for n in ['BrayCurtis','Jaccard','Aitchison']]+['gamma_'+m for m in ['Observed','Shannon','Simpson']]
            plotoutputs=[div/'captions.tsv']+[div/'figures'/fmt/(stem+'.'+fmt) for fmt in ['pdf','png'] for stem in stems]
            runner.step('plot.'+level,[meta,*outputs,*statoutputs],plotoutputs,{'dataset':ds,'level':level},lambda div=div,meta=meta,level=level:plot(div,meta,dc['group_column'],ds,level),deps=['statistics.'+level])
    if mode not in ['dry-run','validate']:
        dump(out/'run_status.json',{'dataset_id':ds,'status':'success','scope':[start,end],'config_sha256':digest(cfg)})
        (out/'last_failure.json').unlink(missing_ok=True)
    return runner

def verify_diversity(folder):
    a=pd.read_csv(folder/'alpha.tsv',sep='\t',index_col=0)
    if not len(a) or not {'Observed','Shannon','Simpson','Chao1'}<=set(a):raise ValueError('invalid alpha table')
    for name in ['BrayCurtis','Jaccard','Aitchison']:
        d=frame(folder/(name+'.distance.tsv'))
        if list(d.index)!=list(d.columns) or list(d.index)!=list(a.index):raise ValueError('distance IDs differ')
        if not np.allclose(d,d.T) or not np.allclose(np.diag(d),0):raise ValueError('distance not symmetric/hollow')

def main(argv=None):
    parser=argparse.ArgumentParser(description='AREK FASTQ → abundance → diversity; strict portable project workflow')
    sub=parser.add_subparsers(dest='command',required=True)
    s=sub.add_parser('scan');s.add_argument('--root',required=True);s.add_argument('--out',required=True)
    for name in ['check','dry-run','run','resume','validate']:
        p=sub.add_parser(name);p.add_argument('--config',required=True);p.add_argument('--dataset',default='all');p.add_argument('--out',default='outputs');p.add_argument('--classifier');p.add_argument('--kraken-db');p.add_argument('--host-index');p.add_argument('--from-stage',choices=STAGES,default='preprocess');p.add_argument('--to-stage',choices=STAGES,default='plot')
    p=sub.add_parser('redraw');p.add_argument('--folder',required=True);p.add_argument('--metadata');p.add_argument('--group-column',default='');p.add_argument('--dataset',default='');p.add_argument('--level',default='')
    args=parser.parse_args(argv)
    try:
        if args.command=='scan':return scan(args.root,args.out)
        if args.command=='redraw':plot(args.folder,args.metadata,args.group_column,args.dataset,args.level);return 0
        cfg=load_config(args.config)
        for arg,key in [('classifier','classifier')]:
            if getattr(args,arg,None):cfg[key]=str(Path(getattr(args,arg)).resolve())
        if args.kraken_db:cfg['shotgun']['kraken_db']=str(Path(args.kraken_db).resolve())
        if args.host_index:cfg['shotgun']['host_index']=str(Path(args.host_index).resolve())
        datasets=read_samples(cfg)
        if cfg['route']=='matrix':datasets={args.dataset if args.dataset!='all' else 'matrix_input':[]}
        elif args.dataset!='all':datasets={ds:datasets.get(ds,[]) for ds in args.dataset.split(',')}
        if not datasets:raise ValueError('no datasets in sample manifest')
        if STAGES.index(args.from_stage)>STAGES.index(args.to_stage):raise ValueError('from-stage is after to-stage')
        failures=[]
        for ds,samples in datasets.items():
            try:
                # Dry run only parses declarations and prints actual commands; no subprocess/no output mutation.
                if args.command!='dry-run':print(json.dumps(check_dataset(cfg,samples,ds),ensure_ascii=False))
                out=Path(args.out).resolve()/ds
                if args.command=='check':continue
                if args.command=='dry-run':execute(cfg,samples,ds,out,args.command,args.from_stage,args.to_stage);continue
                if cfg['route']!='matrix':
                    for s in samples:
                        for k in ['r1','r2']:
                            if s[k] and (out==Path(s[k]) or out in Path(s[k]).parents):raise ValueError('output root overlaps raw inputs')
                # Existing formal project directories are never overwritten by this release entry.
                if (out/'04_export').exists() or (out/'03_vsearch_otu97').exists():raise ValueError('formal legacy result directory detected; choose a new isolated --out')
                if args.command=='validate':
                    lockpath=out/'pipeline.lock'
                    if lockpath.exists():
                        with lockpath.open('r') as lock:
                            try:fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
                            except BlockingIOError:raise ValueError('dataset/tool still running; validation requires a quiescent checkpoint')
                            execute(cfg,samples,ds,out,args.command,args.from_stage,args.to_stage)
                    else:execute(cfg,samples,ds,out,args.command,args.from_stage,args.to_stage)
                else:
                    out.mkdir(parents=True,exist_ok=True)
                    with (out/'pipeline.lock').open('a+') as lock:
                        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                        except BlockingIOError:raise ValueError('dataset already running: pipeline.lock held')
                        execute(cfg,samples,ds,out,args.command,args.from_stage,args.to_stage)
            except Exception as e:
                failures.append({'dataset_id':ds,'error':str(e)});print(f'ERROR {ds}: {e}',file=sys.stderr)
                if args.command not in ['dry-run','check','validate']:dump(Path(args.out)/ds/'last_failure.json',{'error':str(e),'traceback':traceback.format_exc()})
        if failures:print(json.dumps({'partial_failures':failures},ensure_ascii=False),file=sys.stderr);return 2
        return 0
    except Exception as e:print('ERROR:',e,file=sys.stderr);return 2
