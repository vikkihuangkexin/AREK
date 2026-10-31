#!/usr/bin/env python3
"""Audit a source tree or freshly extracted ZIP without external project source."""
from pathlib import Path
import ast,csv,hashlib,json,re,sys,zipfile
root=Path(sys.argv[1]).resolve();manifest=root/'FILE_MANIFEST.tsv';seen=set();count=0
if manifest.exists():
 for row in csv.DictReader(manifest.open(),delimiter='\t'):
  p=root/row['path'];assert p.is_file(),p;assert p.stat().st_size==int(row['bytes']),p
  assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],p;seen.add(row['path']);count+=1
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='FILE_MANIFEST.tsv'}
 assert seen==actual,{'missing':sorted(seen-actual),'extra':sorted(actual-seen)}
code=[]
for p in root.rglob('*'):
 assert not p.is_symlink(),p
 if p.suffix not in ['.py','.R','.sh'] or 'validation' in p.relative_to(root).parts:continue
 code.append(str(p.relative_to(root)));s=p.read_text()
 assert '/data3/khuang6' not in s and '/home/khuang6' not in s,('hardcoded execution path',p)
 if p.suffix=='.py':
  tree=ast.parse(s)
  for node in ast.walk(tree):
   names=[]
   if isinstance(node,ast.Import):names=[n.name for n in node.names]
   elif isinstance(node,ast.ImportFrom) and node.module:names=[node.module]
   for name in names:
    if name.split('.')[0] not in ['workflow','scripts']:continue
    q=root.joinpath(*name.split('.'));assert q.is_dir() or q.with_suffix('.py').is_file(),(p,name)
# Direct ROOT-relative subprocess/script references in the canonical workflow.
for rel in ['scripts/preprocess/host_filter.sh','scripts/diversity/statistics.R','scripts/diversity/compute.py','scripts/diversity/plot.py','scripts/abundance/levels.py','config/historical.json']:
 assert (root/rel).is_file(),rel
for p in (root/'config').glob('*.json'):json.loads(p.read_text())
print(json.dumps({'status':'passed','manifest_payload_files':count,'local_code_files':len(code),'absolute_execution_paths':'none','local_imports':'all resolved','required_script_paths':'all present'},ensure_ascii=False))
