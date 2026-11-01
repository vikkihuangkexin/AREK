# Validation evidence and reproduction

## Evidence origins

| Location | Origin and meaning |
| --- | --- |
| [coverage.tsv](coverage.tsv), [../CODE_INDEX.tsv](../CODE_INDEX.tsv) | Biomni's reported validation boundary for each module |
| `biomni_summary/` | Compact summaries copied from the supplied archive; historical evidence, not a new full rerun |
| [repository_checks.json](repository_checks.json) | Checks performed for the `Alcohol_research/` repository layout |

Detailed per-node execution receipts and their large signatures remain in the original Biomni archive. Its original payload manifest was checked before reorganization. Current code, configuration, documentation, and compact summaries have new hashes in the module-root `FILE_MANIFEST.tsv`.

The repository checks cover Python/shell syntax, all module help entries, dependency ordering, English text, links, source preservation, manifest integrity, and a primary M1 replay against the external original verification input. They do not rerun the 406-node clinical chain, R Beta/PAM, full functions, raw-read profiling, or official mapper execution. Installed versions and unavailable test dependencies are recorded separately from Biomni's pinned environment.

## External-data checks

The original real-data verification input was excluded from the code release. Extract `AC/demo/` from the original Biomni archive or supply the complete frozen M1 root:

```bash
python validation/check_interfaces.py \
  --input-root /path/to/original/AC/demo --output-root /path/to/new/interface_checks

AC_TEST_INPUT_ROOT=/path/to/original/AC/demo python -m pytest -q tests
```

The statistical tests independently compare fixed-family BH and real-project CR1 with statsmodels. They require the original verification inputs and test dependencies; no external-input setting produces explicit skips.

Full scientific replay requires the original external roots:

```bash
python validation/reproduce.py --input-root "$M1" --clinical-root "$CLIN" \
  --v2-root "$V2" --v3-root "$V3" --output-root "$OUT/validation"
python validation/compare_outputs.py --input-root "$M1" --clinical-root "$CLIN" \
  --validation-root "$OUT/validation"
```

## Verify an upload archive

The upload archive contains root `README.md`, root `.gitignore`, and `Alcohol_research/`. Existing `AREK_database/` and `LICENSE` files are preserved in the target repository rather than duplicated in this increment.

```bash
python validation/verify_archive.py --archive /path/to/Alcohol_research_GitHub_upload.zip \
  --extract-to /path/to/new/extracted_release --receipt /path/to/archive_check.json
```

This checks safe paths, CRC, the research manifest, and extracted module help. To additionally replay the original primary M1 input, supply both `--input-root /path/to/original/AC/demo` and `--analysis-output /path/to/new/analysis_output`. The archive receipt stays outside the ZIP to avoid a circular self-hash.

## Unvalidated boundaries

Original raw reads, matching raw-tool environments and references, and optional official-mapper pins are unavailable here. Biomni's raw preflight rejections are not successful raw-read runs. Statistical support failures remain `not_estimable`; see [KNOWN_GAPS.md](../docs/KNOWN_GAPS.md).
