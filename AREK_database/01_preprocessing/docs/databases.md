# Reference databases and reproducibility boundaries

## SILVA classifier

The project classifier is `silva-138-99-nb-classifier.qza`. Its size and SHA256 are recorded in [database_inventory.tsv](database_inventory.tsv). Recovered action provenance records training on 2024-05-30 with q2-feature-classifier 2024.5.0; reference-read UUID `3d95f079-b580-4bce-8c28-8b8a437be0d6` and reference-taxonomy UUID `0e9ef7c4-65a0-408a-859a-b9ce2548e1f8`. Source-package validation loaded and classified it with sklearn 1.4.2 / QIIME2 2024.10.1.

The original archive's `classifier_training_provenance.json` records 7-mers, 8,192 features, alpha .001, and fit_prior=false. The resource belongs to the SILVA 138 99% series. A local filename alone does not establish the exact source URL, primer interval, or region specificity.

The [QIIME2 2024.10 resource page](https://docs.qiime2.org/2024.10/data-resources/) and [QIIME2 Library](https://library.qiime2.org/) are resource entry points. A different downloaded classifier cannot be assumed identical to the historical file. For exact reproduction, restore the verified project backup:

```bash
mkdir -p databases
cp /path/to/verified_backup/silva-138-99-nb-classifier.qza databases/
sha256sum databases/silva-138-99-nb-classifier.qza
# Compare the complete hash with the inventory entry.
qiime tools validate databases/silva-138-99-nb-classifier.qza --level max
```

If the matching reference QZA files are recovered, the recorded method can be rebuilt. Retraining does not guarantee the same UUID or byte hash:

```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /path/to/silva138_99_reference_reads.qza \
  --i-reference-taxonomy /path/to/silva138_99_reference_taxonomy.qza \
  --p-classify--alpha 0.001 --p-classify--chunk-size 20000 --p-no-classify--fit-prior \
  --p-feat-ext--analyzer char_wb --p-feat-ext--n-features 8192 \
  --p-feat-ext--ngram-range '[7, 7]' --p-feat-ext--norm l2 \
  --o-classifier databases/rebuilt-silva138-99.qza
```

Compare all remaining training settings with the archived provenance and locked environment. Missing references should fail explicitly. Large-classifier retraining was not rerun during repository preparation.

## Kraken2 and Bracken

Recovered files include `hash.k2d`, `opts.k2d`, `taxo.k2d`, and `database150mers.kmer_distrib`. File timestamps are 2022-12-05, which does not establish a database release date. `inspect.txt` records nucleotide, k35/l31, 38,236 taxonomy nodes, and table size 1,400,809,609. The database build name, composition, and original URL were not recovered. Exact reproduction requires the same backup and matching hashes.

A different database is a new configuration and needs its own source/build/checksum record. The [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) and [Bracken repository](https://github.com/jenniferlu717/Bracken) describe building resources. This example builds a new database:

```bash
# Requires the route environment and access to the relevant NCBI resources.
export KRAKEN_DB=/path/to/new_kraken_database
kraken2-build --standard --threads 8 --db "$KRAKEN_DB"
bracken-build -d "$KRAKEN_DB" -t 8 -k 35 -l 150
# Retain taxonomy, build logs, source versions, dates, and checksums.
sha256sum "$KRAKEN_DB"/hash.k2d "$KRAKEN_DB"/opts.k2d \
  "$KRAKEN_DB"/taxo.k2d "$KRAKEN_DB"/database150mers.kmer_distrib
```

Bracken distributions must match the configured read length; preflight checks the expected file. P/F/G/S estimates are independent. The workflow retains `bracken_native_relative.tsv`; standardized TSS divides by the rank's total estimated reads. Interpret native fraction/TSS differences using effective reads within each rank. No cross-rank conservation is imposed. Intermediate taxonomy ranks such as G1/F1 are not parsed as G/F; exact-rank parsing is a recorded correction.

## Bowtie2 host reference

The recovered prefix is hg38, with six `.bt2` files and `hg38.fa`; hashes are in the inventory. The exact FASTA snapshot and original URL were not recovered. Restore and verify the same files for historical reproduction. The [UCSC hg38 download directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) contains multiple assembly/analysis/patch resources that are not automatically interchangeable.

```bash
# A newly selected FASTA produces a newly documented index.
bowtie2-build /path/to/verified/hg38.fa /path/to/new_index/hg38
sha256sum /path/to/new_index/hg38.*.bt2
```

Establish the host from study metadata. Human hg38 is unsuitable as a general mouse/rat host reference. Set `host_index` explicitly, or document the historical `skip_host_removal` choice.

## Partial SOLiD/BFAST branch

PRJEB18041 uses ABI SOLiD 5500 color-space shotgun reads. The branch provides diagnostics, hg38 color-space indexing, and chunked host removal. Historical successful evidence is limited to small BFAST/SAM validation logs in `test_solid_bfast`. A complete downstream base-space abundance route has not been demonstrated.

`envs/bfast-0.5.6-source.tar.gz` contains recovered third-party source, excluding server ELF/object files and generated build state. Compilation needs a system C compiler, make, zlib, and libm. The separate legacy-tool environment is recorded in `envs/solid_tools.linux-64.explicit.txt`. BFAST itself must be built from source:

```bash
mkdir -p /path/to/build
cd /path/to/build
tar -xzf /path/to/package/envs/bfast-0.5.6-source.tar.gz
cd bfast-0.5.6
./configure --prefix=/path/to/bfast-0.5.6
make -j2
make install
```

The legacy source dates to the 2010s. Recompilation and full hg38 color-space indexing were not rerun during repository preparation. `legacy/solid/build_index.sh` retains the historical layout, offsets, and scoring generation. Set `BFAST_HOME`, `SOLID_ENV`, `HG38_FASTA`, `BFAST_DB_ROOT`, `RAW_DATA_ROOT`, and `RESULT_ROOT`, and enable `AREK_LEGACY_SOURCE_MODE=isolated_research` only with a new output directory. Missing indices/tools fail explicitly.
