"""Conservative discovery; only explicit mate syntax, never fuzzy pairing/grouping."""
from pathlib import Path
from collections import defaultdict
import csv,gzip,bz2,lzma,re,json
from workflow.common import tsv
import pandas as pd
EXT=re.compile(r'\.(?:fastq|fq)(?:\.(?:gz|bz2|xz))?$',re.I)
MATE=re.compile(r'^(.*?)(?:[_.](?:R|read)?([12]))(?:_001)?$',re.I)

def opener(p):
    return gzip.open if str(p).endswith('.gz') else bz2.open if str(p).endswith('.bz2') else lzma.open if str(p).endswith('.xz') else open

def reads(p,limit=None):
    n=0
    with opener(p)(p,'rt') as f:
        while True:
            h=f.readline()
            if not h:break
            s=f.readline().rstrip();plus=f.readline();q=f.readline().rstrip()
            if not h.startswith('@') or not plus.startswith('+') or not s or len(s)!=len(q):raise ValueError(f'{p}: invalid FASTQ record {n+1}')
            if any(c.isdigit() for c in s):raise ValueError(f'{p}: color-space reads detected; do not send SOLiD to base-space amplicon/shotgun route')
            if any(ord(c)<33 or ord(c)>126 for c in q):raise ValueError(f'{p}: invalid quality at record {n+1}')
            yield h.strip(),s,q;n+=1
            if limit and n>=limit:return
    if not n:raise ValueError(f'{p}: empty FASTQ')

def qc_pair(a,b=None,limit=None):
    import itertools
    n=0;lens=[]
    if b:
        for x,y in itertools.zip_longest(reads(a,limit),reads(b,limit)):
            if x is None or y is None:raise ValueError(f'{a}, {b}: unequal mate record counts')
            norm=lambda s:re.sub(r'/[12]$','',s.split()[0])
            if norm(x[0])!=norm(y[0]):raise ValueError(f'{a}, {b}: mate ID mismatch at {n+1}: {x[0]} vs {y[0]}')
            lens.append((len(x[1]),len(y[1])));n+=1
    else:
        for x in reads(a,limit):lens.append((len(x[1]),0));n+=1
    if not n:raise ValueError(f'{a}: empty FASTQ')
    return {'records_checked':n,'median_R1':float(pd.Series([x[0] for x in lens]).median()),'median_R2':float(pd.Series([x[1] for x in lens]).median())}

def scan(root,out):
    root=Path(root).resolve();out=Path(out);groups=defaultdict(lambda:defaultdict(list));conflicts=[];files=[];artifacts=[]
    if not root.is_dir():raise ValueError(f'root does not exist: {root}')
    for p in sorted(root.rglob('*')):
        if not p.is_file() or out==p or out in p.parents:continue
        rel=p.relative_to(root)
        if EXT.search(p.name):
            # Nearest study ancestor; unconstrained local IDs also supported when explicitly rooted.
            studies=[s for s in rel.parts[:-1] if re.fullmatch(r'(?:PRJ[A-Z]*|SRP|ERP|DRP|CRA)\d+.*',s)]
            ds=studies[-1] if studies else (rel.parts[0] if len(rel.parts)>1 and rel.parts[0]!='fastq' else root.name)
            stem=EXT.sub('',p.name);m=MATE.fullmatch(stem)
            if not m and stem.lower().endswith('_f1'):
                sid=stem[:-3];mate='1'
            else:sid=m[1] if m else stem;mate=m[2] if m else 'SE'
            lane=re.search(r'_L\d{3}$',sid);run=sid[:lane.start()] if lane else sid
            groups[(ds,sid)][mate].append(str(p));files.append({'dataset_id':ds,'sample_id':sid,'run_id':run,'lane':lane[0] if lane else '', 'mate':mate,'path':str(p),'bytes':p.stat().st_size})
        elif p.suffix.lower() in ['.qza','.biom','.tsv','.json'] and any(t in p.name.lower() for t in ['table','tax','alpha','beta','gamma','manifest','contract']):
            artifacts.append({'path':str(p),'bytes':p.stat().st_size,'kind':p.suffix,'status':'present_unverified'})
    samples=[];runlanes=defaultdict(list)
    for (ds,sid),d in groups.items():
        errs=[]
        if any(len(v)>1 for v in d.values()):errs.append('duplicate mate/sample ID: '+json.dumps(dict(d)))
        if 'SE' in d and len(d)>1:errs.append('SE and paired files share ID')
        if ('1' in d)!=('2' in d):errs.append('orphan mate; explicit SE decision required')
        if re.search(r'_L\d{3}$',sid):runlanes[(ds,re.sub(r'_L\d{3}$','',sid))].append(sid)
        row={'dataset_id':ds,'sample_id':sid,'run_id':re.sub(r'_L\d{3}$','',sid),'specimen_id':'','subject_id':'','layout':'PE' if '1' in d and '2' in d else 'SE','r1':(d.get('1') or d.get('SE') or [''])[0],'r2':(d.get('2') or [''])[0],'status':'conflict' if errs else 'discovered_mapping_unconfirmed'}
        samples.append(row)
        for e in errs:conflicts.append({'dataset_id':ds,'sample_id':sid,'error':e})
    for (ds,run),lanes in runlanes.items():
        if len(lanes)>1:conflicts.append({'dataset_id':ds,'sample_id':run,'error':'multiple lanes; provide explicit sample/run/specimen map and premerge decision: '+','.join(lanes)})
    out.mkdir(parents=True,exist_ok=True)
    for name,rows,cols in [('files',files,['dataset_id','sample_id','run_id','lane','mate','path','bytes']),('samples',samples,['dataset_id','sample_id','run_id','specimen_id','subject_id','layout','r1','r2','status']),('conflicts',conflicts,['dataset_id','sample_id','error']),('artifacts',artifacts,['path','bytes','kind','status'])]:tsv(pd.DataFrame(rows,columns=cols),out/(name+'.tsv'),False)
    states=[]
    from workflow.common import sha,validate_file
    for cp in root.rglob('state/*.json'):
        if out in cp.parents or cp.name.endswith('.running.json'):continue
        try:
            rec=json.loads(cp.read_text());status=rec.get('status','unknown')
            if status=='success':
                try:
                    assert rec.get('outputs'),'no declared outputs'
                    for name,expected in rec['outputs'].items():validate_file(name);assert sha(name)==expected,'output hash mismatch'
                    status='completed_content_verified_config_unchecked'
                except Exception as e:status='invalid_output: '+str(e)
            states.append({'dataset_id':rec.get('dataset',cp.parent.parent.name),'stage':rec.get('stage',cp.stem),'status':status})
        except Exception as e:states.append({'dataset_id':cp.parent.parent.name,'stage':cp.stem,'status':'invalid_checkpoint: '+str(e)})
    for ds in sorted({s['dataset_id'] for s in samples}):
        for stage in ['import','features','taxonomy','levels','diversity','plot']:
            states.append({'dataset_id':ds,'stage':stage,'status':'conflict' if any(c['dataset_id']==ds for c in conflicts) else 'not_verified'})
    tsv(pd.DataFrame(states,columns=['dataset_id','stage','status']),out/'step_status.tsv',False)
    print(f'{len(samples)} samples, {len(conflicts)} conflicts; {out}')
    return 2 if conflicts else 0
