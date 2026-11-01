"""Verify the GitHub upload increment; optionally replay external frozen primary inputs."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import zipfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", required=True)
    ap.add_argument("--extract-to", required=True, help="A fresh, nonexistent directory.")
    ap.add_argument("--receipt", required=True, help="External JSON receipt, outside the archive.")
    ap.add_argument("--input-root", help="Optional external original primary verification input.")
    ap.add_argument("--analysis-output", help="Fresh output directory, required with --input-root.")
    a = ap.parse_args()
    if bool(a.input_root) != bool(a.analysis_output):
        ap.error("Supply --input-root and --analysis-output together.")
    archive, target = Path(a.archive).resolve(), Path(a.extract_to).resolve()
    analysis_out = Path(a.analysis_output).resolve() if a.analysis_output else None
    if target.exists() or (analysis_out and analysis_out.exists()):
        raise FileExistsError("Extraction and analysis destinations must be fresh.")
    invalid = re.compile(r'[<>:"\\|?*\x00-\x1f]')
    reserved = {"CON", "PRN", "AUX", "NUL"} | {
        prefix+str(i) for prefix in ("COM", "LPT") for i in range(1, 10)}
    root_files = {"README.md", ".gitignore"}
    prefix = "Alcohol_research/"
    manifest_name = prefix + "FILE_MANIFEST.tsv"
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        names = [x.filename for x in entries]
        if len(names) != len(set(names)) or len(names) != len({n.casefold() for n in names}):
            raise ValueError("Duplicate or case-colliding ZIP names.")
        for entry in entries:
            path = PurePosixPath(entry.filename)
            if path.is_absolute() or ".." in path.parts or not path.parts:
                raise ValueError("Unsafe ZIP path: " + entry.filename)
            if entry.filename not in root_files and path.parts[0] != "Alcohol_research":
                raise ValueError("Unexpected increment path: " + entry.filename)
            if entry.is_dir() or entry.compress_type != zipfile.ZIP_DEFLATED:
                raise ValueError("Expected regular DEFLATE payloads only.")
            if (entry.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Symlinks are not allowed.")
            for part in path.parts:
                if invalid.search(part) or part.endswith((" ", ".")) or part.split(".")[0].upper() in reserved:
                    raise ValueError("Not a Windows-safe member: " + entry.filename)
            if len(entry.filename) > 180:
                raise ValueError("Archive path exceeds the release limit.")
        if z.testzip() is not None:
            raise ValueError("ZIP CRC failure.")
        rows = list(csv.DictReader(io.StringIO(z.read(manifest_name).decode()), delimiter="\t"))
        expected = {row["path"]: row for row in rows}
        actual = {n[len(prefix):] for n in names if n.startswith(prefix) and n != manifest_name}
        if len(expected) != len(rows) or set(expected) != actual:
            raise ValueError("Research manifest membership mismatch.")
        if set(names) != {prefix+p for p in expected} | {manifest_name} | root_files:
            raise ValueError("Unexpected or missing increment payloads.")
        z.extractall(target)
    module = target / "Alcohol_research"
    for relative, row in expected.items():
        path = (module / relative).resolve()
        if module not in path.parents or not path.is_file():
            raise ValueError("Unsafe or missing manifest path: " + relative)
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
            raise ValueError("Research payload mismatch: " + relative)
    preparation = json.loads((module / "history/repository_preparation.json").read_text())
    for row in preparation["root_update_files"]:
        path = target / row["path"]
        if row["path"] not in root_files or path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
            raise ValueError("Root update mismatch: " + row["path"])
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
        PYTHONDONTWRITEBYTECODE="1")
    scripts = sorted((module / "scripts").glob("*.py"))
    for script in [module / "run.py", *scripts]:
        result = subprocess.run([sys.executable, str(script), "--help"],
            capture_output=True, text=True, env=env, cwd=target)
        if result.returncode:
            raise RuntimeError(script.name + ": " + result.stdout + result.stderr)
    replay = {"status": "NOT_RUN", "reason": "No external frozen verification input supplied."}
    if analysis_out:
        command = [sys.executable, str(module / "run.py"), "--modules", "11",
            "--input-root", str(Path(a.input_root).resolve()), "--output-root", str(analysis_out)]
        result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=target)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        node = next(analysis_out.glob("11__*/execution.json"))
        detail = json.loads(node.read_text())["detail"]
        for field, value in {"fixed_family": 360, "estimable": 182, "candidates": 39,
                "positive": 15, "negative": 24, "full_reference_features": 3201}.items():
            if detail[field] != value:
                raise AssertionError(field)
        driver = json.loads((analysis_out / "driver_status.json").read_text())
        if [x["module"] for x in driver] != ["01", "07", "08", "11"] or any(x["status"] != "PASS" for x in driver):
            raise AssertionError("Extracted dependency driver did not pass.")
        detail["numerical_reference"] = "External frozen primary reference table"
        replay = {"status": "PASS", "driver_modules": ["01", "07", "08", "11"], "summary": detail}
    receipt = {"status": "PASS", "archive": archive.name, "bytes": archive.stat().st_size,
        "sha256": digest(archive), "compression": "DEFLATE", "research_root": prefix,
        "research_payload_files": len(expected), "root_update_files": sorted(root_files),
        "zip_members": len(names), "all_hashes_and_sizes_verified": True,
        "crc": "PASS", "safe_paths": "PASS", "duplicates": 0, "unlisted_payloads": 0,
        "max_member_path_length": max(map(len, names)), "independent_extraction": True,
        "entrypoint_help_checks": len(scripts)+1, "external_primary_replay": replay,
        "python": sys.version}
    receipt_path = Path(a.receipt).resolve()
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print("UPLOAD_ARCHIVE_PASS", len(names), "files", receipt["sha256"])


if __name__ == "__main__":
    main()
