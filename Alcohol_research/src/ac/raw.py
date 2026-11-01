"""Fail-closed launcher for the recovered source-specific read pipelines."""
from .io import *
import shutil
import subprocess

def run(args, out, stage):
    if not args.source_config:
        raise ValueError("BLOCKED: provide --source-config with original reads, database paths and output roots.")
    conf = json.loads(Path(args.source_config).read_text())
    source = conf["adapter"]
    if (stage == "04" and source != "amplicon_v5") or (stage == "05" and source != "shotgun_v2"):
        raise ValueError("The source adapter does not match this profiling stage.")
    required = ["RAW_ROOT", "RESULT_ROOT", "CLASSIFIER"] if source == "amplicon_v5" else [
        "RAW_ROOT", "OUT_ROOT", "MISS_ROOT", "HOST_BT2_DB", "KRAKEN_DB"]
    env = os.environ.copy()
    for field in required:
        if field not in conf:
            raise ValueError("BLOCKED: missing source configuration: "+field)
        env[field] = str(Path(conf[field]).resolve())
    env.update({k: str(v) for k, v in conf.get("settings", {}).items()})
    if str(conf.get("settings", {}).get("FORCE", "0")) != "0":
        raise ValueError("FORCE must be zero; use a new source-specific output directory.")
    if source == "amplicon_v5":
        commands = ["python", "qiime", "biom"]
        if not Path(env["CLASSIFIER"]).is_file():
            raise FileNotFoundError("BLOCKED: original SILVA classifier is absent.")
    elif source == "shotgun_v2":
        commands = ["python3", "fastp", "bowtie2", "samtools", "pigz", "kraken2", "bracken"]
        if not Path(env["KRAKEN_DB"]).is_dir():
            raise FileNotFoundError("BLOCKED: original Kraken2/Bracken database is absent.")
        if not list(Path(env["HOST_BT2_DB"]).parent.glob(Path(env["HOST_BT2_DB"]).name+"*.bt2*")):
            raise FileNotFoundError("BLOCKED: original hg38 Bowtie2 index is absent.")
    else:
        raise ValueError("Unsupported recovered source adapter.")
    missing = [x for x in commands if not shutil.which(x)]
    if missing:
        raise RuntimeError("BLOCKED: missing executables: "+", ".join(missing))
    raw = Path(env["RAW_ROOT"])
    if not raw.is_dir():
        raise FileNotFoundError("BLOCKED: raw FASTQ directory is absent.")
    destination = Path(env["RESULT_ROOT"] if source == "amplicon_v5" else env["OUT_ROOT"])
    if destination == raw or raw in destination.parents or destination in raw.parents:
        raise ValueError("Source output must not overlap the raw input directory.")
    for field in ("RESULT_ROOT", "OUT_ROOT", "MISS_ROOT"):
        if field in env and field in conf and not Path(env[field]).is_relative_to(out.resolve()):
            raise ValueError("Configure every source output root inside this module's isolated output directory.")
    env["AC_STOP_AFTER_QC"] = "1" if stage == "03" else "0"
    script = Path(__file__).with_name(source+".sh")
    command = ["bash", str(script)]
    if source == "amplicon_v5":
        command += [conf["project"]]
    else:
        env["PROJECT"] = conf["project"]
    dump(dict(command=command, config=conf, stage=stage), out/"source_execution.json")
    subprocess.run(command, env=env, check=True)
    return dict(adapter=source, stage=stage, raw_FASTQ_validation="RECOMPUTED")
