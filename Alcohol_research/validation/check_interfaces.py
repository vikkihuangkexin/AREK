"""Exercise overwrite protection, resume, dry-run and blocked raw preflights."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE/"src"))
from ac.io import dump, sha

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", required=True,
        help="External original primary verification input, such as the source archive's AC/demo.")
    ap.add_argument("--output-root", required=True)
    a = ap.parse_args()
    out = Path(a.output_root).resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    checks = []
    command = [sys.executable, str(PACKAGE/"scripts/11_alcohol_vs_health_M1.py"),
        "--input-root", str(Path(a.input_root).resolve()), "--output-root", str(out/"main")]
    def run(name, cmd, predicate):
        result = subprocess.run(cmd, env=env, text=True, capture_output=True)
        passed = bool(predicate(result))
        checks.append(dict(check=name, status="PASS" if passed else "FAIL", command=cmd,
            returncode=result.returncode, output=(result.stdout+result.stderr)[-2000:]))
        if not passed:
            raise AssertionError(name+": "+result.stdout+result.stderr)
        return result
    run("fresh_external_primary_replay", command, lambda r:r.returncode == 0)
    run("overwrite_refused", command, lambda r:r.returncode != 0 and "Output exists" in r.stderr)
    run("verified_resume", command+["--resume"], lambda r:r.returncode == 0 and "REUSED_HASH_CHECKED" in r.stdout)
    run("explicit_overwrite", command+["--overwrite"], lambda r:r.returncode == 0)
    receipt = next((out/"main").glob("*/execution.json"))
    before = sha(receipt)
    run("existing_node_dry_run_no_writes", command+["--dry-run", "--overwrite"],
        lambda r:r.returncode == 0 and "DRY_RUN_NO_WRITES" in r.stdout)
    assert sha(receipt) == before
    driver = [sys.executable, str(PACKAGE/"run.py"), "--modules", "11",
        "--input-root", str(Path(a.input_root).resolve()), "--output-root", str(out/"driver_dry"), "--dry-run"]
    result = run("dependency_dry_run", driver, lambda r:r.returncode == 0 and r.stdout.count("COMMAND") == 4)
    assert not (out/"driver_dry").exists()
    for module, source in [("03","source_amplicon.example.json"),("04","source_amplicon.example.json"),("05","source_shotgun.example.json")]:
        script = next((PACKAGE/"scripts").glob(module+"_*.py"))
        run("raw_preflight_"+module, [sys.executable, str(script), "--input-root", str(Path(a.input_root).resolve()),
            "--output-root", str(out/("blocked_"+module)), "--source-config", str(PACKAGE/"configs"/source)],
            lambda r:r.returncode != 0 and "BLOCKED" in r.stderr)
    dump(checks, out/"interface_checks.json")
    print("INTERFACE_CHECKS_PASS", len(checks))

if __name__ == "__main__":
    main()
