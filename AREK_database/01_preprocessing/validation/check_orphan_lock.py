#!/usr/bin/env python3
"""Kill only a private validation coordinator; verify its live tool retains the dataset lock."""
from pathlib import Path
import sys,subprocess,time,json,fcntl,os
P=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True
w=Path(sys.argv[1]).resolve();w.mkdir(parents=True,exist_ok=True)
helper=w/'private_lock_probe.py'
helper.write_text('''from pathlib import Path
import sys,fcntl,os
sys.dont_write_bytecode=True
sys.path.insert(0,sys.argv[1])
from workflow.common import Runner
out=Path(sys.argv[2]);out.mkdir(exist_ok=True)
with (out/'pipeline.lock').open('a+') as f:
 fcntl.flock(f,fcntl.LOCK_EX)
 r=Runner(out,{'threads':1,'memory_gb':1},'PRIVATE_LOCK_TEST','run');r.current='private_probe';r.record={'commands':[],'pid':os.getpid(),'status':'running'}
 r.cmd([sys.executable,'-c','import time;time.sleep(5)'])
''')
out=w/'output';parent=subprocess.Popen([sys.executable,str(helper),str(P),str(out)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
checkpoint=out/'state/private_probe.running.json';deadline=time.time()+30
while not checkpoint.exists() and time.time()<deadline:
 if parent.poll() is not None:raise AssertionError(parent.stderr.read().decode())
 time.sleep(.05)
assert checkpoint.exists();record=json.loads(checkpoint.read_text());parent.kill();parent.wait()
with (out/'pipeline.lock').open('r') as f:
 try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);raise AssertionError('orphan tool lost dataset lock')
 except BlockingIOError:pass
# Wait only for this test's five-second child to finish naturally.
deadline=time.time()+15;released=False
while time.time()<deadline:
 with (out/'pipeline.lock').open('r') as f:
  try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);released=True;break
  except BlockingIOError:pass
 time.sleep(.1)
assert released,'private tool did not release lock on exit'
result={'status':'passed','only_private_coordinator_SIGKILL':True,'recorded_tool_pid':record['child_pid'],'lock_held_after_coordinator_death':True,'duplicate_execution_blocked':True,'lock_released_after_tool_exit':True}
(w/'orphan_lock_test.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
