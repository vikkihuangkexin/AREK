# Parameter sources and method differences

`config/historical.json` stores parseable defaults. Historical method evidence and portable resource settings are distinguished below; portable configurations are not exact parameter-by-parameter replays of historical runs.

| Parameter | Value | Source or status |
| --- | --- | --- |
| OTU identity | 0.97 | v5 shell and PRJDB10075 clustered QZA provenance |
| Median read length | >=100 bp | v5 layout checks sampled the first 3,000 reads; main preparation checks the full file and requires an explicit mate |
| PE minovlen / maxdiffs / minlen | 10 / 10 / 1 | Recovered merged-QZA action; truncqual/maxns/maxee/minmergelen/maxmergelen were null; allowmergestagger=false |
| Derep min_seq_length / min_unique_size | 1 / 1 | Recovered derep QZA; derep_prefix=false |
| UCHIME dn / mindiffs / mindiv / minh / xn | 1.4 / 3 / 0.8 / 0.28 / 8 | Recovered nonchimeras QZA |
| Classification | sklearn, confidence=0.7, read_orientation=auto | Taxonomy QZA; 16 historical jobs versus two threads in verification settings |
| SILVA training | 138 99%; q2-feature-classifier 2024.5; 7-mer; 8,192 features; alpha=.001; fit_prior=false | Classifier provenance; complete recovered settings in `classifier_training_provenance.json` in the original archive; primer/region identity not established |
| fastp v2 PE | Minimum length 50; correction=false; duplication evaluation disabled when supported; compression=1 | PRJNA613834 v2 source; archived `project_run_parameters.tsv` records actual settings including threads=8, fastp=4, read_len=150; source threads=16 is not an executed value |
| fastp original PE | Minimum length 50; detect_adapter_for_pe; correction=true | PRJCA021196 source; original SE omits correction and the PE adapter flag |
| fastp quality settings | Locked tool defaults | Original commands did not explicitly record all quality cutoffs; no complete historical per-option snapshot was recovered |
| Bowtie2 PE | very-sensitive; samtools -f12 -F256; name sort; fastq -n | Original/v2 source; both mates unmapped and nonsecondary; portable nonhost outputs may be uncompressed with unchanged read content |
| Bowtie2 SE | very-sensitive; samtools -f4 -F256 | Original source; v2 does not support SE |
| Kraken2 | paired/use-names; confidence default 0 | Original command omitted confidence; zero is the locked tool default, not a recovered custom value |
| Bracken | Read length 150; threshold 10; independent P/F/G/S | Shell and real logs; ranks are independently estimated rather than summed from Species; no forced cross-rank total equality |
| Filtering | Prevalence .10; mean relative abundance 1e-4 | `02_run_taxonomy_analysis_from_standardized.py`; full-lineage grouping fixes old same-name collisions and can change taxon counts |
| TSS | Original all-feature denominator plus a separate retained-taxon denominator | Both retained and audited |
| Alpha | Observed; natural-log Shannon; 1-D Simpson; bias-corrected Chao1 | Recovered analysis_core functions; Chao1 restricted to integer counts |
| Beta | Direct Bray on selected values; binary Jaccard; Euclidean CLR of ln(x+.5) | Recovered analysis_core; matrix units determine pseudocount meaning |
| Gamma | Sum group columns, then S/H/1-D; independent label permutations with plus-one p values | Recovered `run_manifest_microbiome_analysis.py` functions; unsupported repeated/blocked/adjusted tests return NA |
| Permutations | 999; seed=1 | Historical Python defaults; 99/19 in verification are new small-test settings |
| Rarefaction | None | No rarefaction in the main Python workflow |
| Phylogeny / NMDS | None | No implementation in the selected executed standardized route |
| Existing merged-matrix branch | Seed 20261005; 999 permutations; unfiltered Genus relative Alpha and filtered-TSS Beta; Aitchison as defined in the R source | Separate from the main counts route |
| Threads / memory | Defaults two / 8 GiB; matrix example 1 GiB | Portable preflight settings, not historical resources or measured RSS |

Input diagnostics, full-lineage aggregation, descriptive missing-group output, content checkpoints, distance/coordinate/eigenvalue export, true PERMDISP, and repeated-subject guards are recorded portable-workflow corrections. The old invalid betadisper function is retained only for audit.

No defaults are invented for unrecovered primers, regions, read-to-subject mappings, or database build identities. A historical longer-mate decision can be represented by an explicit SE manifest only with supporting read-layout evidence. Discovery does not infer the biological mate from length.

The main shotgun workflow retains declared clean/nonhost intermediates for checkpoint validation. This differs from the old v2 `CLEAN_INTERMEDIATE=1` operational policy without changing read content or classification algorithms; account for the disk requirement. Cleanup/FORCE controls are operational settings.

The historical standardized caller also required `min-samples=6` and `min-group-size=3`. The portable workflow retains descriptions for smaller or unlabeled inputs rather than canceling the module. Its corrected statistical entry requires at least two effective observations per group as an estimability threshold. To retain a study's historical inclusion rules, supply verified inclusion metadata and inspect group sizes. Four-sample software tests do not reproduce the original study analysis. These differences, asymptotic R p values, and the PERMDISP correction are explicitly recorded behavior changes from the historical source, predating the English-language update.
