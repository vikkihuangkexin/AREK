# Environment inventory

| Environment | Purpose |
| --- | --- |
| `analysis` | Matrix processing, main-workflow diversity, and plotting |
| `amplicon` | QIIME2/vsearch FASTQ-to-diversity route |
| `shotgun` | fastp/Bowtie2/Kraken2/Bracken-to-diversity route |
| `integration` | Standalone R route for existing merged Genus matrices, including ggplot2 |
| `solid_tools` | Recovered SOLiD Python/samtools/make/zlib dependencies; BFAST built separately |

The explicit locks contain required conda metadata seed packages and their transitive dependency closure for Linux x86_64. Each `*.versions.json` records actual versions, builds, and package SHA256/MD5. Recovered missing-dependency arrays are empty. These files describe installed environments; clean network restoration was not rerun during repository preparation. Download failures should be explicit, without silently substituting newer packages.

The system must provide bash, GNU coreutils, awk/sed/grep, gzip, tar, and time. BFAST compilation additionally needs a C compiler and zlib development headers. [system_tools.tsv](system_tools.tsv) records recovered system tool versions. Third-party BFAST retains its bundled license; recompilation was not rerun.

Keep BFAST's legacy dependencies separate from QIIME2. Activate the corresponding environment at its new installation location; code does not depend on the old server prefix. Environments, references, and paths contribute to checkpoint fingerprints, so changes require validation or rerunning affected stages.
