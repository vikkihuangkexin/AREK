"""Dependency-ordered module driver; profile inputs do not require a raw pipeline."""
from pathlib import Path
import argparse
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"src"))
from ac.cli import NAMES, parser

DEPS = {"07": ["01"], "08": ["01"], "09": ["07", "08"], "10": ["07", "08"],
    "11": ["07", "08"], "12": ["07", "08"], "13": ["07", "08"],
    "14": ["01"], "15": ["01"], "17": ["01"], "18": ["01"],
    "19": ["07", "08"], "20": ["01"], "21": ["01"], "22": ["01"]}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modules", default="11", help="Comma-separated module numbers; upstream profile dependencies are added.")
    ap.add_argument("--dry-run", action="store_true", help="Print exact dependency-ordered commands without running.")
    ap.add_argument("--continue-on-error", action="store_true", help="Finish independent modules and record failed/blocked nodes.")
    options, remainder = ap.parse_known_args()
    requested = [x.zfill(2) for x in options.modules.split(",")]
    order = []
    def visit(node):
        if node not in NAMES:
            ap.error("Unknown module: "+node)
        for upstream in DEPS.get(node, []):
            visit(upstream)
        if node not in order:
            order.append(node)
    for node in requested:
        visit(node)
    supported = set().union(*(set(parser(node)._option_string_actions) for node in order))
    for token in remainder:
        if token.startswith("--") and token not in supported:
            ap.error("Unknown module option: "+token)
    records = []
    failed = set()
    for node in order:
        known_options = parser(node)._option_string_actions
        passed = []
        i = 0
        while i < len(remainder):
            key = remainder[i]
            if key in known_options:
                passed.append(key)
                action = known_options[key]
                if action.nargs != 0 and i+1 < len(remainder):
                    passed.append(remainder[i+1])
                    i += 1
            elif not key.startswith("--"):
                ap.error("Unrecognized standalone argument: "+key)
            elif i+1 < len(remainder) and not remainder[i+1].startswith("--"):
                i += 1
            i += 1
        command = [sys.executable, str(ROOT/"scripts"/(node+"_"+NAMES[node]+".py"))]+passed
        if any(x in failed for x in DEPS.get(node, [])):
            failed.add(node)
            records.append(dict(module=node, status="BLOCKED", reason="upstream_failure"))
            continue
        print("COMMAND", json.dumps(command), flush=True)
        if options.dry_run:
            continue
        result = subprocess.run(command)
        records.append(dict(module=node, status="PASS" if result.returncode == 0 else "FAIL",
            returncode=result.returncode))
        if result.returncode:
            failed.add(node)
            if not options.continue_on_error:
                break
    if not options.dry_run:
        p = argparse.ArgumentParser(add_help=False)
        p.add_argument("--output-root", required=True)
        output, _ = p.parse_known_args(remainder)
        dest = Path(output.output_root)
        dest.mkdir(parents=True, exist_ok=True)
        (dest/"driver_status.json").write_text(json.dumps(records, indent=2)+"\n")
    return int(bool(failed))

if __name__ == "__main__":
    raise SystemExit(main())
