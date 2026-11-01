"""Tabular I/O, immutable inputs, and content-addressed execution receipts."""
from pathlib import Path
import hashlib
import json
import os
import time
import numpy as np
import pandas as pd

def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)

def save(frame, path, index=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=index, na_rep="NA", float_format="%.17g")

def dump(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + "\n", encoding="utf-8")

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            h.update(block)
    return h.hexdigest()

def hash_ids(ids):
    return hashlib.sha256(("\n".join(map(str, ids)) + "\n").encode()).hexdigest()

def framehash(frame):
    return hashlib.sha256(pd.util.hash_pandas_object(frame, index=True).values.tobytes()
        + "|".join(map(str, frame.columns)).encode()).hexdigest()

def truth(series):
    return series.astype(str).str.lower().eq("true")

def require_columns(frame, columns):
    absent = set(columns) - set(frame.columns)
    if absent:
        raise ValueError("Missing required columns: " + ", ".join(sorted(absent)))
    if not frame.index.is_unique:
        raise ValueError("Duplicate specimen identifiers are not permitted.")

def file_manifest(root):
    root = Path(root)
    return [dict(path=p.relative_to(root).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(root.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]

def execute(args, module, inputs, callback):
    """Protect inputs and reuse only receipts with matching code and input hashes."""
    root = Path(args.output_root).resolve()
    for field in ("input_root", "clinical_root", "source_root", "database_root"):
        value = getattr(args, field, None)
        if value:
            source = Path(value).resolve()
            if root == source or source in root.parents:
                raise ValueError("Output must be outside each authoritative input root.")
    key = {k: v for k, v in vars(args).items()
        if k not in {"output_root", "resume", "overwrite", "dry_run"}}
    paths = sorted(set(Path(p).resolve() for p in inputs))
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("BLOCKED: missing input files: " + "; ".join(missing))
    code_root = Path(__file__).parents[2]
    code = {p.relative_to(code_root).as_posix(): sha(p)
        for folder in ("src", "scripts", "configs")
        for p in sorted((code_root/folder).rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts}
    signature = dict(arguments=key, inputs={str(p): sha(p) for p in paths}, code=code)
    signature_hash = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    target = root/module
    receipt = target/"execution.json"
    if getattr(args, "dry_run", False):
        print(json.dumps(dict(module=module, inputs=[str(p) for p in paths], output=str(target),
            existing_output=target.exists(), action="DRY_RUN_NO_WRITES")))
        return None
    if target.exists():
        if args.resume and receipt.exists():
            previous = json.loads(receipt.read_text())
            intact = all((target/x["path"]).is_file() and
                (target/x["path"]).stat().st_size == x["bytes"] and
                sha(target/x["path"]) == x["sha256"] for x in previous.get("outputs", []))
            if previous["signature_hash"] == signature_hash and intact:
                print("REUSED_HASH_CHECKED", module)
                return previous
        if not args.overwrite:
            raise FileExistsError("Output exists; choose --resume, --overwrite, or a new output root: " + str(target))
        import shutil
        shutil.rmtree(target)
    target.mkdir(parents=True)
    start = time.perf_counter()
    try:
        detail = callback(target) or {}
    except Exception as error:
        dump(dict(status="FAIL", error=str(error)), target/"failure.json")
        raise
    record = dict(status="PASS", validation_mode="RECOMPUTED", module=module,
        signature_hash=signature_hash, signature=signature, elapsed_seconds=time.perf_counter()-start,
        detail=detail, outputs=file_manifest(target))
    dump(record, receipt)
    print("PASS", module, round(record["elapsed_seconds"], 3), "seconds")
    return record
