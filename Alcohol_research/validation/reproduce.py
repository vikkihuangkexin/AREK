"""Execute syntax checks and real-data module replays in an isolated output root."""
from pathlib import Path
import argparse
import ast
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE/"src"))
from ac.cli import NAMES
from ac.io import dump, read, save

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", required=True)
    ap.add_argument("--clinical-root", required=True)
    ap.add_argument("--output-root", required=True)
    ap.add_argument("--v2-root", help="Optional original genus/function v2 root for module 16.")
    ap.add_argument("--v3-root", help="Optional original v3 source input root for full-composition alpha.")
    ap.add_argument("--database-root", help="Optional original FAPROTAX-mapper_v2.0 directory.")
    ap.add_argument("--modules", default="all", help="Comma-separated IDs or all.")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    output = Path(a.output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    selected = set(NAMES) if a.modules == "all" else {x.zfill(2) for x in a.modules.split(",")}
    status = []
    previous = output/"validation_status.json"
    if previous.exists():
        prior = json.loads(previous.read_text())
        replayed = {x["module"] for x in prior if x.get("check") == "scientific_replay" and x["status"] == "PASS"}
        status += [x for x in prior if
            x.get("check") == "scientific_replay" and x["module"] not in selected
            and x["module"] not in {"03", "04", "05"}
            and not (x["status"] == "NOT_RUN" and x["module"] in replayed)]
    for path in PACKAGE.rglob("*.py"):
        ast.parse(path.read_text(encoding="utf-8"))
    for path in (PACKAGE/"src/ac").glob("*.sh"):
        subprocess.run(["bash", "-n", str(path)], check=True)
    for path in (PACKAGE/"src/ac").glob("*.R"):
        subprocess.run(["Rscript", "-e", "parse(file="+json.dumps(str(path))+")"],
            check=True, stdout=subprocess.DEVNULL)
    for module, name in NAMES.items():
        command = [sys.executable, str(PACKAGE/"scripts"/(module+"_"+name+".py")), "--help"]
        start = time.perf_counter()
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        status.append(dict(module=module, check="syntax_import_help", status="PASS" if result.returncode == 0 else "FAIL",
            elapsed_seconds=time.perf_counter()-start, command=command))
        if result.returncode:
            raise RuntimeError(result.stderr)
    print("All 22 syntax/import/help checks passed.", flush=True)
    jobs = [
        ("01", "audit", []),
        ("02", "metadata", []),
        ("06", "rank_export", []),
        ("07", "profile", ["--source-root", a.v3_root] if a.v3_root else []),
        ("08", "membership", []),
        ("09", "alpha", ["--source-root", a.v3_root] if a.v3_root else []),
        ("11", "main", []),
    ]
    jobs += [("12", "disease_"+str(i), ["--comparison", c]) for i, c in enumerate(
        ["02_AUD_or_ALD_vs_Health", "03_AUD_vs_Health", "04_ALD_vs_Health", "05_ALD_vs_non_ALD"], 1)]
    jobs += [
        ("13", "subset_supplemented", []),
        ("13", "subset_original", ["--scope", "M1_clinical_complete"]),
        ("14", "clinical_full", ["--variable", "all", "--panel", "all"]),
        ("14", "clinical_single", ["--variable", "sex", "--body", "gut", "--scope", "frozen",
            "--outcome", "taxa", "--panel", "gut__full_source_reference", "--coverage", "base"]),
        ("15", "source_main", []),
        ("15", "source_function", ["--outcome", "function"]),
        ("17", "function_base", []),
        ("17", "function_global", ["--coverage", "global_coverage"]),
        ("17", "function_trait", ["--coverage", "trait_coverage"]),
        ("18", "contributions", []),
        ("10", "beta_shared", []),
        ("19", "PAM_ALD", []),
        ("19", "PAM_Alcohol", ["--comparison", "01_Alcohol_vs_Health"]),
        ("20", "overlaps", []),
        ("21", "figures", []),
        ("22", "presentation", []),
    ]
    if a.v2_root:
        extra = ["--source-root", a.v2_root]
        if a.database_root:
            extra += ["--database-root", a.database_root, "--rerun-mapper"]
        jobs.append(("16", "annotation", extra))
    elif "16" in selected:
        status.append(dict(module="16", check="scientific_replay", status="NOT_RUN",
            reason="Original FAPROTAX source root was not supplied."))
    for module in ("03", "04", "05"):
        status.append(dict(module=module, check="scientific_replay", status="BLOCKED",
            reason="Original reads, databases and required profiling executables are not available in the supplied completed archives."))
    for module, label, extra in jobs:
        if module not in selected:
            continue
        command = [sys.executable, str(PACKAGE/"scripts"/(module+"_"+NAMES[module]+".py")),
            "--input-root", str(Path(a.input_root).resolve()), "--clinical-root", str(Path(a.clinical_root).resolve()),
            "--output-root", str(output/label)] + extra
        if a.overwrite:
            command.append("--overwrite")
        start = time.perf_counter()
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        elapsed = time.perf_counter()-start
        log = output/(label+".log")
        log.write_text(result.stdout+"\n"+result.stderr, encoding="utf-8")
        record = dict(module=module, check="scientific_replay", label=label,
            status="PASS" if result.returncode == 0 else "FAIL", validation_mode="RECOMPUTED",
            elapsed_seconds=elapsed, command=command, log=str(log))
        status.append(record)
        dump(status, output/"validation_status.json")
        print(module, label, record["status"], round(elapsed, 2), "seconds", flush=True)
        if result.returncode:
            print(result.stderr[-2500:], flush=True)
    environment = dict(python=sys.version, platform=platform.platform(),
        packages={p: importlib.metadata.version(p) for p in
            ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "matplotlib", "seaborn",
             "beautifulsoup4", "pytest"]}, thread_limits={k: env[k] for k in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]})
    dump(environment, output/"environment.json")
    dump(status, output/"validation_status.json")
    save(__import__("pandas").DataFrame([{k: v for k, v in x.items() if k != "command"} for x in status]), output/"validation_status.tsv")
    return int(any(x["status"] == "FAIL" for x in status))

if __name__ == "__main__":
    raise SystemExit(main())
