# 参考库与可复现边界

## 实际分类器

`silva-138-99-nb-classifier.qza` 来自当前项目；SHA256/bytes 见 database_inventory.tsv。其自身 action 显示 2024-05-30 用 q2-feature-classifier 2024.5.0 训练，reference_reads UUID 为 3d95f079-b580-4bce-8c28-8b8a437be0d6，reference_taxonomy UUID 为 0e9ef7c4-65a0-408a-859a-b9ce2548e1f8；不是依据文件名猜训练版本。本轮在 sklearn1.4.2/QIIME2 2024.10.1 成功加载并分类。

原始 ZIP 中的 classifier_training_provenance.json 保存历史恢复参数（未上传），包括 7-mer、8192 features、alpha .001、fit_prior false。源数据是 SILVA138 99% 系列；精确原始下载 URL/引物区间不能仅由本地文件名证明，未恢复内容不补默认。

原 [QIIME2 2024.10 资源页](https://docs.qiime2.org/2024.10/data-resources/) 目前将分类器导向新的 [QIIME2 Library](https://library.qiime2.org/)；不能把当前新版本文件当作历史同一分类器。历史精确复现的首选是从本项目参考库备份复制该文件，在目标机器执行：

```bash
mkdir -p databases
cp /path/to/verified_backup/silva-138-99-nb-classifier.qza databases/
sha256sum databases/silva-138-99-nb-classifier.qza
# 与 docs/database_inventory.tsv 的该行逐字比对
qiime tools validate databases/silva-138-99-nb-classifier.qza --level max
```

若恢复了对应 reference QZA，可以按已恢复方法重新训练；这些命令是**重建方法**，不承诺重新训练得到同一 UUID/字节 hash，不使用待证实的网络地址：

```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /path/to/silva138_99_reference_reads.qza \
  --i-reference-taxonomy /path/to/silva138_99_reference_taxonomy.qza \
  --p-classify--alpha 0.001 --p-classify--chunk-size 20000 --p-no-classify--fit-prior \
  --p-feat-ext--analyzer char_wb --p-feat-ext--n-features 8192 \
  --p-feat-ext--ngram-range '[7, 7]' --p-feat-ext--norm l2 \
  --o-classifier databases/rebuilt-silva138-99.qza
```

训练命令与训练参考文件缺失时准确报错；不下载一个不同数据库假装恢复。训练其余参数需查原始 ZIP 的 provenance，固定环境内缺省与原记录比对。本轮未重新训练大型分类器。

## Shotgun Kraken2/Bracken

当前实际库有 `hash.k2d`、`opts.k2d`、`taxo.k2d`、`database150mers.kmer_distrib`。mtime 为 2022-12-05，不等于已证实的库发布日期。inspect.txt 显示 nucleotide、k35/l31、38236 taxonomy nodes、table size1400809609；**库名称、组成及原下载 URL 未恢复**，不能命名为某个 standard 最新版。已计算实际文件 SHA256，可从同一备份复制后核验。

要在新机器使用不同库，应单独记录来源/构建版本/校验并作为新配置；官方的[Kraken2 手册](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)和[Bracken 源码说明](https://github.com/jenniferlu717/Bracken)提供构建方式。下列是新库构建示例，**不是历史精确复现**：

```bash
# 路线环境已安装 kraken2-build/bracken-build，网络可访问对应 NCBI 来源时
export KRAKEN_DB=/path/to/new_kraken_database
kraken2-build --standard --threads 8 --db "$KRAKEN_DB"
bracken-build -d "$KRAKEN_DB" -t 8 -k 35 -l 150
# 保留下载库与 taxonomy、构建日志、下载日期、来源版本及文件 hash
sha256sum "$KRAKEN_DB"/hash.k2d "$KRAKEN_DB"/opts.k2d \
  "$KRAKEN_DB"/taxo.k2d "$KRAKEN_DB"/database150mers.kmer_distrib
```

不能把不同 read-length 的 Bracken distribution 替换成150；预检明确查找对应文件。P/F/G/S 是独立重估，不拿 species 表推导其他所有 rank。保存 native `bracken_native_relative.tsv`，规范 TSS 用本 rank 估计 reads 总和，native fraction 与 TSS 的差异应依据各 rank 的有效 reads 解释，不强行跨 rank 守恒。source taxonomy 的 G1/F1 等中间 rank 不再误认作 G/F，已记录为修复。

## 宿主 Bowtie2

当前 prefix 为 hg38，有6个 bt2及 hg38.fa，SHA256见清单。精确 FASTA 快照/原下载 URL 未恢复；历史复现需复制同一文件并比对 hash。[UCSC hg38 下载目录](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)可用于获取新的参考，目录中 assembly/analysis set/patch版本不同，不能无证据互换。

```bash
# 如采用新的已明确版本的 FASTA，记录该选择；不是历史同一索引。
bowtie2-build /path/to/verified/hg38.fa /path/to/new_index/hg38
sha256sum /path/to/new_index/hg38.*.bt2
```

宿主类型必须由研究 metadata 确认。人 hg38 不是鼠/大鼠数据的通用去宿主库；规范入口不按目录自动换宿主，也不能用示例 human 配置处理其他物种。显式配置 host_index 或说明历史未去宿主的 skip_host_removal 选择。

## SOLiD/BFAST 部分支路

这是 PRJEB18041 ABI SOLiD5500 的 color-space shotgun。提供诊断、hg38 CS建索引及分块去宿主源码；成功证据仅为 test_solid_bfast 中的小型 BFAST/SAM validation 日志。没有证明可用的后续 base-space 丰度链，因此其分类/多样性从 FASTQ 分支为**未完成历史支路**，不能拿标准 ASV/OTU 入口替代。

`envs/bfast-0.5.6-source.tar.gz` 为当前项目可恢复的第三方源代码，去掉服务器 ELF/objects/生成的构建状态。库依赖为系统 C编译器、make、zlib、libm；独立旧工具环境精确包闭包见 `solid_tools.linux-64.explicit.txt`。BFAST 自身不是 conda 包，需单独从所附源码编译：

```bash
mkdir -p /path/to/build
cd /path/to/build
tar -xzf /path/to/package/envs/bfast-0.5.6-source.tar.gz
cd bfast-0.5.6
./configure --prefix=/path/to/bfast-0.5.6
make -j2
make install
```

原源码为2010年代版本；本轮未重新编译及重建巨大 hg38 CS索引。`build_index.sh` 保留原布局、offset和scoring生成方法。运行前设置 BFAST_HOME、SOLID_ENV、HG38_FASTA、BFAST_DB_ROOT、RAW_DATA_ROOT、RESULT_ROOT，并只在新的输出目录开启 `AREK_LEGACY_SOURCE_MODE=isolated_research`。索引/网络/工具缺失会明确报错，不能用服务器外部旧脚本路径代替源码。
