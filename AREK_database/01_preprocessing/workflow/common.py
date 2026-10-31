"""Portable strict readers and content-addressed local execution."""
from pathlib import Path
import csv, hashlib, json, os, subprocess, sys, time, zipfile, fcntl, shutil, signal
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,default=str).encode()).hexdigest()
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.tmp');q.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n');q.replace(p)
def tsv(df,p,index=True):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,sep='\t',index=index,na_rep='NA')
def frame(p):
    with Path(p).open() as f:
        skip=0
        for line in f:
            if line.startswith('# Constructed') or line.startswith('# Generated'):skip+=1
            else:break
    # Reject duplicate headers before pandas can silently rename them.
    with Path(p).open() as f:
        for _ in range(skip):next(f)
        header=next(csv.reader(f,delimiter='\t'))
    if len(header)!=len(set(header)):raise ValueError(f'{p}: duplicate column IDs')
    x=pd.read_csv(p,sep='\t',skiprows=skip,index_col=0,keep_default_na=False,dtype=str)
    x.index=x.index.astype(str)
    if not len(x) or not len(x.columns):raise ValueError(f'{p}: empty matrix')
    if x.index.duplicated().any() or (x.index=='').any():raise ValueError(f'{p}: duplicate/empty row IDs')
    try:x=x.apply(pd.to_numeric,errors='raise')
    except Exception as e:raise ValueError(f'{p}: missing/structural/non-numeric cell (never filled with zero): {e}')
    if not np.isfinite(x.values).all() or (x.values<0).any():raise ValueError(f'{p}: non-finite or negative values')
    return x

def taxonomy(p):
    x=pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
    if x.shape[1]<2:raise ValueError(f'{p}: taxonomy needs feature ID and Taxon columns')
    if x.iloc[:,0].duplicated().any():raise ValueError(f'{p}: duplicate taxonomy IDs')
    taxcol=next((c for c in x if c.lower() in ['taxon','taxonomy','consensus lineage']),None)
    if not taxcol:raise ValueError(f'{p}: missing Taxon column')
    return pd.Series(x[taxcol].values,index=x.iloc[:,0].astype(str))

def version_snapshot():
    import importlib.metadata as m
    names=['numpy','pandas','scipy','matplotlib','qiime2','q2-vsearch','q2-feature-classifier','scikit-learn','biom-format']
    result={'python':sys.version,'executables':{n:shutil.which(n) for n in ['qiime','Rscript','fastp','bowtie2','kraken2','bracken','samtools']}}
    for n in names:
        try:result[n]=m.version(n)
        except m.PackageNotFoundError:result[n]='missing'
    # Conda metadata records tools whose wrappers lack consistent --version.
    prefix=Path(sys.prefix)/'conda-meta'
    for p in prefix.glob('*.json'):
        try:
            d=json.loads(p.read_text())
            if d['name'] in ['vsearch','fastp','bowtie2','kraken2','bracken','samtools','r-base','r-vegan','r-permute']:result[d['name']]=d['version']+'='+d.get('build','')
        except Exception:pass
    return result

def validate_file(p):
    p=Path(p)
    if not p.is_file() or p.stat().st_size==0:raise ValueError(f'missing/empty output: {p}')
    if p.suffix in ['.qza','.qzv']:
        with zipfile.ZipFile(p) as z:
            if z.testzip() or not any('/data/' in n for n in z.namelist()):raise ValueError(f'invalid QIIME archive: {p}')
    if p.suffix=='.json':json.loads(p.read_text())

class Runner:
    def __init__(self,out,config,dataset,mode):
        self.out=Path(out);self.config=config;self.dataset=dataset;self.mode=mode
        self.versions=version_snapshot();self.current=None;self.rebuilt=set()
        self.code={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.suffix in ['.py','.R','.sh'] and not any(x in p.parts for x in ['validation','examples'])}
        self.record=None
    def fingerprint(self,tag,inputs,params):
        if tag.startswith('plot.'):paths=['scripts/diversity/plot.py']
        elif tag.startswith('statistics.'):paths=['scripts/diversity/statistics.R']
        elif tag.startswith('diversity.'):paths=['scripts/diversity/compute.py','scripts/diversity/historical_metrics.py','scripts/diversity/historical_gamma.py']
        elif tag.startswith('levels.'):paths=['scripts/abundance/levels.py']
        elif tag.startswith('metadata.'):paths=[]
        else:paths=['workflow/routes.py','workflow/scanner.py','scripts/preprocess/host_filter.sh','scripts/abundance/historical_bracken_merge.py']
        paths+=['workflow/common.py','workflow/cli.py','run_pipeline.py']
        code={p:self.code[p] for p in paths if p in self.code}
        return digest({'inputs':{str(p):sha(p) for p in inputs},'parameters':params,'resources':{'threads':self.config['threads'],'memory_gb':self.config['memory_gb']},'code':code,'tools':self.versions})
    def event(self,x):
        p=self.out/'execution.jsonl';p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('a') as f:f.write(json.dumps({'dataset':self.dataset,'time':time.time(),**x},ensure_ascii=False,default=str)+'\n')
    def cmd(self,args):
        args=list(map(str,args));self.record['commands'].append(args);self.event({'stage':self.current,'command':args,'event':'command_start'})
        env=os.environ.copy();env.update(PYTHONNOUSERSITE='1',OMP_NUM_THREADS=str(self.config['threads']),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        tmp=self.out/'tmp';tmp.mkdir(exist_ok=True);env.update(TMPDIR=str(tmp),TMP=str(tmp),TEMP=str(tmp))
        log=self.out/'logs'/f'{self.current}.log';log.parent.mkdir(exist_ok=True)
        with log.open('a') as f:
            f.write('\nCOMMAND '+json.dumps(args)+'\n');f.flush()
            # Keep the dataset flock alive in the actual tool if the coordinator is killed.
            # Linux flock follows the inherited open-file description, including across exec.
            lockfds=[]
            for fd in Path('/proc/self/fd').iterdir():
                try:
                    if Path(os.readlink(fd)).resolve()==(self.out/'pipeline.lock').resolve():lockfds.append(int(fd.name))
                except (OSError,ValueError):pass
            p=subprocess.Popen(args,stdout=f,stderr=subprocess.STDOUT,env=env,pass_fds=tuple(lockfds),start_new_session=True)
            self.record['child_pid']=p.pid;dump(self.out/'state'/f'{self.current}.running.json',self.record)
            try:rc=p.wait()
            except BaseException:
                # An interactive interrupt must not leave workers writing into an unlocked output.
                try:os.killpg(p.pid,signal.SIGTERM)
                except ProcessLookupError:pass
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:os.killpg(p.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    p.wait()
                raise
        self.event({'stage':self.current,'command':args,'exit_code':rc,'event':'command_end'})
        if rc:raise RuntimeError(f'{self.current}: exit {rc}; see {log}')
    def step(self,tag,inputs,outputs,params,action,validator=None,deps=(),plan=None):
        group='plot' if tag.startswith('plot.') else 'diversity' if tag.startswith(('metadata.','diversity.','statistics.')) else 'levels' if tag.startswith('levels.') else 'taxonomy' if tag in ['classify','taxonomy'] else 'features' if tag in ['features','shotgun_merge_P','shotgun_merge_F','shotgun_merge_G','shotgun_merge_S'] else 'preprocess'
        stages=['preprocess','features','taxonomy','levels','diversity','plot']
        if stages.index(group)>stages.index(getattr(self,'end','plot')):return
        upstream_only=stages.index(group)<stages.index(getattr(self,'start','preprocess'))
        effective_mode='validate' if upstream_only and self.mode!='dry-run' else self.mode
        inputs=list(map(Path,inputs));outputs=list(map(Path,outputs));state=self.out/'state'/f'{tag}.json'
        if self.mode=='dry-run':
            print(json.dumps({'dataset':self.dataset,'stage':tag,'inputs':list(map(str,inputs)),'outputs':list(map(str,outputs)),'commands':plan or ['internal Python: '+tag],'parameters':params},ensure_ascii=False,default=str));return
        fp=self.fingerprint(tag,inputs,params);valid=False;reason='checkpoint missing'
        if state.exists():
            try:
                old=json.loads(state.read_text());assert old['fingerprint']==fp,'input/code/config/tool fingerprint changed'
                assert old['status']=='success','not successful'
                for p in outputs:validate_file(p);assert old['outputs'].get(str(p))==sha(p),f'output changed: {p}'
                if validator:validator()
                assert not set(deps)&self.rebuilt,'upstream was restored'
                valid=True
            except Exception as e:reason=str(e)
        if effective_mode=='validate':
            if not valid:raise ValueError(f'{self.dataset}/{tag}: invalid or stale: {reason}')
            print('VALID',self.dataset,tag);return
        if self.mode=='resume' and valid:
            self.event({'stage':tag,'event':'reused','fingerprint':fp});print('REUSE',self.dataset,tag);return
        self.rebuilt.add(tag);self.current=tag
        self.record={'dataset':self.dataset,'stage':tag,'fingerprint':fp,'input_sha256':{str(p):sha(p) for p in inputs},'expected_outputs':list(map(str,outputs)),'parameters':params,'configuration_sha256':digest(params),'effective_configuration_sha256':digest(self.config),'code_sha256':self.code,'tools':self.versions,'resources':{'threads':self.config['threads'],'memory_gb':self.config['memory_gb']},'started':time.time(),'pid':os.getpid(),'commands':[],'status':'running'}
        state.parent.mkdir(parents=True,exist_ok=True);dump(state.with_name(tag+'.running.json'),self.record)
        print('RUN',self.dataset,tag,flush=True)
        try:
            # Remove obsolete declared products when a stage switches applicable/not-applicable output contracts.
            if state.exists():
                previous=json.loads(state.read_text())
                for name in (previous.get('outputs') or previous.get('expected_outputs',[])):
                    stale=Path(name)
                    if stale not in outputs and self.out.resolve() in stale.resolve().parents and stale.is_file():stale.unlink()
            for p in outputs:
                p.parent.mkdir(parents=True,exist_ok=True)
                if p.exists():p.unlink() # only declared outputs under isolated --out, lock held
            action()
            for p in outputs:validate_file(p)
            if validator:validator()
            self.record.update(status='success',exit_code=0,validation='content+semantic validation passed',outputs={str(p):sha(p) for p in outputs})
        except BaseException as e:
            self.record.update(status='failed',exit_code=1,error=str(e));raise
        finally:
            self.record['ended']=time.time();dump(state,self.record);self.event(self.record)
            state.with_name(tag+'.running.json').unlink(missing_ok=True)
