#!/usr/bin/env python3
"""Verify relocation and matrix numerics without QIIME2, R or reference databases."""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from workflow.cli import load_config, read_samples
from workflow.common import frame
from scripts.abundance.levels import aggregate, verify
from scripts.diversity.compute import compute


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', required=True)
    args = parser.parse_args()
    checks = []

    def record(name, **details):
        checks.append({'check': name, 'status': 'passed', **details})

    # The supplied file map preserves original hashes even after relocation.
    rows = list(csv.DictReader((ROOT / 'docs/provenance/file_relocation.tsv').open(), delimiter='\t'))
    identical = [r for r in rows if r['status'] == 'identical']
    for r in identical:
        p = ROOT / r['repository_module_path']
        assert hashlib.sha256(p.read_bytes()).hexdigest() == r['source_sha256'], p
    record('source_bytes_preserved', files=len(identical))

    py = list(ROOT.rglob('*.py'))
    shell = list(ROOT.rglob('*.sh'))
    for p in py:
        ast.parse(p.read_text())
    for p in shell:
        subprocess.run(['bash', '-n', str(p)], check=True, capture_output=True)
    record('python_and_shell_syntax', python_files=len(py), shell_files=len(shell))

    configs = sorted((ROOT / 'config').glob('*.json'))
    for p in configs:
        cfg = load_config(p)
        for key in ['samples', 'matrix', 'taxonomy', 'metadata']:
            if cfg[key]:
                assert Path(cfg[key]).is_file(), (p, key)
    pe = read_samples(load_config(ROOT / 'config/real_pe_example.json'))
    assert len(pe['PRJDB10075']) == 4
    for row in pe['PRJDB10075']:
        assert Path(row['r1']).is_file() and Path(row['r2']).is_file()
    record('config_and_manifest_relative_paths', configs=len(configs), real_pe_samples=4)

    with tempfile.TemporaryDirectory(prefix='arek_repository_check_') as temporary:
        work = Path(temporary)

        def cli(*argv, expected=0):
            p = subprocess.run([sys.executable, str(ROOT / 'run_pipeline.py'), *map(str, argv)],
                               cwd=work, capture_output=True, text=True)
            assert p.returncode == expected, p.stdout + p.stderr
            return p.stdout

        cli('--help')
        for command in ['scan', 'check', 'dry-run', 'run', 'resume', 'validate', 'redraw']:
            cli(command, '--help')
        record('cli_help_from_external_working_directory', commands=8)

        for name, ds in [('matrix_example.json', 'REAL_MATRIX'),
                         ('real_pe_example.json', 'PRJDB10075'),
                         ('real_se_example.json', 'PRJDB10075'),
                         ('shotgun_v2_example.json', 'PRJCA021196'),
                         ('shotgun_original_example.json', 'PRJCA021196')]:
            output = work / ('planned_' + name)
            text = cli('dry-run', '--config', ROOT / 'config' / name,
                       '--dataset', ds, '--out', output)
            assert text.strip() and not output.exists()
        record('dry_run_all_routes_no_writes', configurations=5)

        cli('scan', '--root', ROOT / 'examples/real_fastq', '--out', work / 'scan')
        samples = pd.read_csv(work / 'scan/samples.tsv', sep='\t')
        assert len(samples) == 4 and samples.layout.eq('PE').all()
        record('real_fastq_scan', paired_samples=4)

        cfg = load_config(ROOT / 'config/matrix_example.json')
        for level in ['Phylum', 'Genus']:
            abundance = work / level / 'abundance'
            abundance.mkdir(parents=True)
            aggregate(cfg['matrix'], cfg['taxonomy'], level, abundance, cfg['abundance'], 'counts', None)
            verify(abundance)
            matrix = abundance / 'filtered_counts.tsv'
            x = frame(matrix)
            diversity = work / level / 'diversity'
            diversity.mkdir()
            compute(matrix, diversity, 'counts', cfg['diversity'])
            a = pd.read_csv(diversity / 'alpha.tsv', sep='\t', index_col=0)
            for sid in x.columns:
                v = x[sid].to_numpy()
                probabilities = v[v > 0] / v.sum()
                f1, f2 = int(sum(v == 1)), int(sum(v == 2))
                reference = [sum(v > 0), -sum(float(p) * math.log(float(p)) for p in probabilities),
                             1 - sum(float(p) ** 2 for p in probabilities),
                             sum(v > 0) + f1 * (f1 - 1) / (2 * (f2 + 1))]
                assert np.allclose(a.loc[sid, ['Observed', 'Shannon', 'Simpson', 'Chao1']],
                                   reference, atol=1e-12, rtol=0)
            u, v = x.iloc[:, 0].to_numpy(), x.iloc[:, 1].to_numpy()
            reference = sum(abs(float(i) - float(j)) for i, j in zip(u, v)) / sum(
                float(i) + float(j) for i, j in zip(u, v))
            d = pd.read_csv(diversity / 'BrayCurtis.distance.tsv', sep='\t', index_col=0)
            assert abs(d.iloc[0, 1] - reference) < 1e-12
            record('matrix_aggregation_and_independent_numerics', level=level,
                   samples=x.shape[1], features=x.shape[0], absolute_tolerance=1e-12)

        relative = frame(cfg['matrix'])
        relative = relative / relative.sum()
        relative_path = work / 'relative.tsv'
        relative.to_csv(relative_path, sep='\t')
        out = work / 'relative'
        out.mkdir()
        aggregate(relative_path, cfg['taxonomy'], 'Genus', out, cfg['abundance'], 'relative', None)
        verify(out)
        assert not list(out.glob('*counts.tsv'))
        div = work / 'relative_diversity'
        div.mkdir()
        compute(out / 'filtered_tss.tsv', div, 'relative', cfg['diversity'])
        assert pd.read_csv(div / 'alpha.tsv', sep='\t').Chao1.isna().all()
        record('relative_only_contract', no_counts=True, Chao1='NA')

    report = {
        'date_America_Chicago': '2026-10-31', 'status': 'passed', 'checks': checks,
        'external_tools_available': {t: bool(shutil.which(t)) for t in ['Rscript', 'qiime', 'kraken2']},
        'not_rerun': ['FASTQ external-tool chains', 'R statistical tests', 'full study cohorts',
                      'clean environment installation'],
    }
    target = Path(args.report).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'checks': len(checks), 'report': str(target)}))


if __name__ == '__main__':
    main()
