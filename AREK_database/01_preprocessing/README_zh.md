# AREK：FASTQ → 分类层级丰度 → 多样性

这是 AREK 数据库的第一阶段代码，来源于 2026-10-30 的 FASTQ 到丰度与多样性源码包，仓库整理日期为 2026-10-31（America/Chicago）。简明入口见 [README.md](README.md)，步骤与模块对应见 [docs/steps_zh.md](docs/steps_zh.md)。所有数据、参考库和输出目录均从配置或命令行传入。历史验证记录与当前仓库检查分别登记，不能把历史执行视为本次重跑。

## 范围与来源

规范入口是 `python run_pipeline.py`。它串联已恢复源码，提供本地串行调度，每个工具按线程配置运行。没有新增集群框架，不会提交或停止服务器其他任务。默认拒绝向旧 `04_export`/`03_vsearch_otu97` 正式结果目录写入。

| 路线/入口 | 真实方法 | 支持范围 |
|---|---|---|
| `route=amplicon_vsearch97`，`layout=PE` | QIIME2 导入 → vsearch merge-pairs → dereplicate → 97% de novo OTU → UCHIME de novo → 去嵌合表/序列 → sklearn SILVA 分类 | 完整 FASTQ 至多样性；本轮真实生物 PE 子集验证 |
| 同一路线，`layout=SE` | 明确的生物学 read → 导入 → 相同 OTU97/分类 | 完整 SE 链；本轮以真实 PE 数据中**显式选取 R1**验证，不冒充原研究历史 SE 结果 |
| `route=shotgun_kraken_bracken`，`shotgun.variant=v2` | fastp → Bowtie2 去宿主 → Kraken2 → 各 rank 独立 Bracken | 原 PRJNA613834 PE 路线；代码/环境/数据库检查，未全链重算 |
| `shotgun.variant=original` | 原 PRJCA021196 fastp（PE correction）及 PE/SE 去宿主路线 | PE/SE 源码已恢复；本轮未全链执行 |
| `route=matrix` | 现成 feature 表+taxonomy，或明确 `native_level` 的本级矩阵 → 丰度/多样性 | counts 和 0–1 相对丰度；真实矩阵已验证；**不是从该研究 FASTQ 重现** |
| `legacy/pooled_diversity/integrated_diversity.R` | 原有合并 Genus 矩阵的独立多样性路线 | 只读已核实的合并矩阵；不自动混合研究/物种/部位；未对完整合并队列重跑 |
| `legacy/amplicon/selected_se_batch.sh`、`legacy/amplicon/PRJCA022064.sh` | 原批量 SE derep/特殊命名及 best-mate 原生支路 | 保留原源码和路径适配，供隔离研究；不具备规范入口的内容续跑保证，未执行 |
| `scripts/preprocess/solid/` | ABI SOLiD 5500 color-space 诊断、BFAST 0.5.6 hg38 去宿主 | **部分路线**，没有可证明成功的后续 base-space 分类链；不能声称 FASTQ 至丰度完整复现 |

`docs/source_script_map.tsv` 列源码 SHA256、当前入口、依赖和修订。历史日志、QZA provenance、全量研究清单和结果快照仅保留在原始 ZIP，不上传代码仓库。

版本选择不依据 v3/latest 字样：PRJDB10075 的 v5 日志有成功结束，但目前 table/taxonomy QZA 的实际 provenance 为 2026-06-26/27，晚于 5 月日志。PRJNA909793 的 v5 日志有 SIGSEGV，且原批次混有 WGS/AMPLICON；不能沿用混合输入。实际策略必须由样本表人工确认。只把项目已有的 OTU97 路线打通，未把 shotgun 套成 ASV，也未引入新的去噪、去引物算法。

源码中缺少某研究的引物、扩增区域或平台证据时，明确为**未恢复**，不从文件名推测。当前统一 OTU97 源码没有独立 cutadapt/质量截断步骤，PE 的质量相关参数为已恢复的 merge 默认（见下文），SE 不额外加造质量筛选。若输入仍含需去除的引物，必须另行确认该研究的真实前处理；本包不声称已去引物。

不包含差异菌、回归模块、功能预测、FAPROTAX、PICRUSt2、MICOM、患者分型或预测模型。多样性内部的协变量校正只用于指定的组间距离检验。

## 安装、版本与数据库

已验证平台：Linux x86_64。环境按路线拆开，`envs/*.linux-64.explicit.txt` 是从实际 conda-meta 提取所需种子包及**传递依赖闭包**，包含下载 URL、构建和校验，不是整个服务器环境。`*.versions.json` 保存版本/构建及包 SHA256/MD5。不会要求旧 BFAST 与 QIIME2 安装在同一环境。

```bash
cd AREK_database/01_preprocessing
# 根据使用范围选一个环境，名称可修改。
conda create -n arek_amplicon --file envs/amplicon.linux-64.explicit.txt
conda activate arek_amplicon
# 只做矩阵/多样性时：
# conda create -n arek_analysis --file envs/analysis.linux-64.explicit.txt
# shotgun：
# conda create -n arek_shotgun --file envs/shotgun.linux-64.explicit.txt
python run_pipeline.py --help
```

关键已安装版本：Python 3.10.20，QIIME2 2024.10.1，q2-vsearch/q2-feature-table/q2-feature-classifier 2024.10.0，scikit-learn 1.4.2，numpy 1.26.4，pandas 2.2.2，scipy 1.13.0，matplotlib 3.8.4，R vegan 2.7-1，permute 0.9-8。R 和系统工具精确版本以 `envs/*.versions.json` 为准。旧 create_env 脚本引用的 2026.4 链接**不是**目前运行环境的证据，故没有直接采用。

环境锁定文件来源经过核对，但本轮没有联网新建这些环境。远端包下架/网络失败应在 conda 明确报错，不会自动替换版本。其他平台需重新求解并验证，不直接使用 Linux-64 explicit 文件。

参考库不随 ZIP 分发。详见 `docs/databases_zh.md` 和 `docs/database_inventory.tsv`，其中记录本机文件 SHA256、来源证据与未恢复信息。

```bash
# 使用自己放置的已核对 SILVA 分类器（无需改代码）
export AREK_CLASSIFIER=/path/to/silva-138-99-nb-classifier.qza
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
```

示例资源为**新配置**：2 线程、8 GiB 内存预检请求、输出盘至少剩余 1 GiB。它不限制进程 RSS，不能当成实测峰值。源包真实小样本耗时见原始 ZIP 的 `validation/final_tests/real_chain_measurements.json`（未上传）；峰值 RSS 未测量。全量扩增子聚类可能需几十至数百 GiB，与 unique 序列数量有关；shotgun 还需加载 Kraken 库与宿主索引。这些是估计，应先小规模测试并按实际资源修改配置。未因验证调整其他正在运行的作业。

## 输入、扫描与映射

支持 `.fastq/.fq` 及 `.gz/.bz2/.xz` 扫描；原 QIIME2/fastp 工具是否接受相应压缩格式仍以预检/工具为准，推荐未压缩或 gzip。明确 mate 后缀如 `_R1_001/_R2_001`、`_R1/_R2`、`_1/_2`、`.1/.2`、`_read1/_read2`、`_f1/_r2`。无 mate 后缀视为 SE 候选。多 lane (`_L001` 等) 不会自动合并。孤立 mate、重复 ID、同 ID 同时有 SE/PE 都输出冲突。扫描不能判断 run 是否对应同一标本/个体，必须补充明确映射；发现结果仅标为待确认。

```bash
python run_pipeline.py scan --root /path/to/data_root --out /path/to/scan_report
# 本包真实小数据扫描：
python run_pipeline.py scan --root examples/real_fastq --out /tmp/arek_scan
```

扫描输出：`files.tsv`、`samples.tsv`、`conflicts.tsv`、`artifacts.tsv`、`step_status.tsv`。旧文件仅标“存在未验证”；找到规范 checkpoint 时核验其内容 hash，并注明尚未核对当前配置。不能凭文件存在判断完成，最终要执行 `validate`。扫描有冲突返回 2，但仍输出所有独立数据集发现结果。

实际运行样本表为 TSV，**相对路径基于该样本表位置**。dataset_id 是按真实研究、宿主物种、部位、处理流程分开的分析单元；例如同研究的 saliva/feces 应分别登记，不把 scan 候选清单直接当成已核实的生物学映射：

```tsv
dataset_id	sample_id	run_id	specimen_id	subject_id	layout	r1	r2
PRJDB10075	DRR234531	DRR234531			PE	real_fastq/PRJDB10075/DRR234531_1.fastq.gz	real_fastq/PRJDB10075/DRR234531_2.fastq.gz
```

表中的 `\t` 表示制表符，真实文件见 `examples/real_samples.tsv`。必需列为 dataset_id/sample_id/run_id/layout/r1/r2；SE 的 r2 必须空，PE 必须有两个 mate。specimen_id、subject_id 在有证据时填写，空值不是“已确认独立”。不同原始 run 对应同一标本或多个 lane 时，规范入口要求先明确技术重复处理，不擅自合并或重命名消冲突。独立统计分析前必须确认 subject，不能以 run ID 代替个体 ID。

本地不限制 PRJNA 编号，支持 PRJDB、PRJCA、SRP、CRA 及安全的自定义 ID；允许字母、数字、点、横线、下划线。不是网页 API 的 PRJ 编号白名单。

Metadata 用 TSV，第一列唯一 sample_id，其他列由配置指定；保留未匹配/多余样本审计，不要求所有样本有组别。示例 `examples/metadata_template.tsv` 是字段示意，不是真实分组。可用 group、subject_id、study、host_species、body_site、age、sex 等列。manifest 与 metadata 的 run/specimen/subject 信息冲突时定位报错。配置所指列不存在、ID 无交集或重复会报错。

统计字段位于 `diversity`：`group_column`、`subject_column`、`block_column`、`covariates`（列名数组）、`design`。`design` 可为 descriptive/independent/within_subject。未知组值通过 `missing_group_values` 显式列出，保留 source_group_value 后转为缺失，**绝不转成 Control**。未指定组别仍输出全部描述结果。

配置文件相对路径基于**配置文件位置**；`--out`/命令行数据库路径基于当前工作目录；代码模块相对包 root。新配置以 `config/historical.json` 合并缺省键，未知键报错。`historical.json` 的方法参数来自历史实现/provenance，资源与空路径是可移植新配置，详细区分见 `docs/parameters_zh.md`。固定指标数学契约另列 `config/metric_contract.tsv`，与可选择的矩阵/分组/置换参数分开。

## 完整命令

先安装/激活环境并设置 AREK_CLASSIFIER，然后执行：

```bash
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
python run_pipeline.py dry-run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
python run_pipeline.py run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
python run_pipeline.py validate --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
```

`dry-run` 只打印实际外部命令、内部函数阶段及参数，不启动工具、不写目录、不提交任务。不存在的未来 taxonomy 无法提前断定各 level 是否适用，会在真实执行时判断。

同一配置样本表可含多数据集，分别执行；`--dataset all` 为批量，或传逗号列表。一个失败不阻断其他数据集，批末有任何失败返回 2，全部成功返回 0。`--dataset PRJDB10075` 为单数据集。不同布局/路线或研究设计请使用各自配置，不能用同一设置强行混跑。

```bash
# 原样续跑：只复用 hash、内容验证与当前输入/代码/配置/工具指纹一致的阶段。
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe
# 已有有效前处理 checkpoint，仅执行指定阶段范围：
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_pe \
  --from-stage levels --to-stage diversity
# 真实 R1 显式 SE 入口，独立输出：
python run_pipeline.py run --config config/real_se_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_example_se
```

范围按 preprocess/features/taxonomy/levels/diversity/plot 排序。跳过的上游仍需通过现有 checkpoint 核验；若手头只有表，不要伪造 checkpoint，使用矩阵入口。`run` 对选定阶段重新执行，`resume` 按内容复用。输出必须为独立目录，未完成的同数据集用 flock 防重复运行；退出/中断自动释放内核锁，running JSON 保存 PID 和子进程 PID，成功才生成成功 checkpoint。工具子进程继承数据集 flock，调度进程异常死亡时仍阻止重复启动；交互中断会终止其独立工具进程组。`validate` 也会拒绝仍有运行锁的目录。不要根据残留 PID/空 running 文件推断作业活着。中断子进程应确认整个进程树退出后再续跑；本包本地等待子进程，不自动重新提交。

现成 feature 矩阵无需分类器：

```bash
python run_pipeline.py check --config config/matrix_example.json --dataset REAL_MATRIX
python run_pipeline.py run --config config/matrix_example.json --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py validate --config config/matrix_example.json --dataset REAL_MATRIX --out /tmp/arek_matrix
```

改配置的 matrix/taxonomy 指向用户表，即可从该阶段接入。有 taxonomy 按实际提供的 rank 聚合；只有 Genus 表时设 `native_level: "Genus"`、`levels: ["Genus"]` 并保留源分类键，不编造上下级。只有相对丰度时设 `matrix_unit: "relative"`、`diversity.input: "filtered_tss"` 或 `"relative"`。输入必须是明确的 0–1、列和为 1 的完整分母；百分比或截断表会报错，需根据真实分母在上游**显式**处理，不能静默 TSS 掩盖缺项。

Shotgun 配置见 `config/shotgun_v2_example.json` 与 `shotgun_original_example.json`；用户填入真实样本表，用命令行指定参考库：

```bash
python run_pipeline.py check --config config/shotgun_v2_example.json \
  --dataset PRJCA021196 --kraken-db /path/to/kraken_db --host-index /path/to/hg38/hg38
python run_pipeline.py dry-run --config config/shotgun_v2_example.json \
  --dataset PRJCA021196 --kraken-db /path/to/kraken_db --host-index /path/to/hg38/hg38 \
  --out /path/to/new_shotgun_output
# check 通过后将 dry-run 换成 run，再用 validate；示例表中的 user_fastq 是明确待提供的输入。
```

不依赖 FASTQ 的重画入口（读真实 TSV，重新输出 PDF/PNG/caption；不改变指标）：

```bash
python run_pipeline.py redraw --folder /tmp/arek_matrix/REAL_MATRIX/diversity/Genus \
  --dataset REAL_MATRIX --level Genus
# 要按已有 metadata 标色：增加 --metadata /path/to/metadata.tsv --group-column group
```

直接 redraw 会更新图/caption，可能使原 plot checkpoint 失效，这是预期；需按原配置 `resume --from-stage plot` 恢复规范图及其指纹。

## 输出与指标

所有输出位于 `--out/DATASET_ID/`：

- `upstream/`：原始导入、拼接/去重复/聚类/嵌合 QZA、最终 feature BIOM/TSV、代表序列 FASTA 与 taxonomy；shotgun 按样本保存 clean/nonhost FASTQ、分类报告及 P/F/G/S 的独立 Bracken 表。
- `abundance/LEVEL/`：feature×sample 的 counts（Bracken 为估计 reads）、0–1 relative、filtered_counts 或 filtered_relative_input、filtered_tss、filtered_original_denominator、完整谱系分类键、feature→taxon 映射、filter_audit、mass_balance 和 contract。无 counts 的输入绝不产生 counts 文件。未支持层级只有 `not_applicable.json`。
- `diversity/LEVEL/`：逐样本 alpha、三种距离矩阵、全轴坐标与全部 eigenvalue、Gamma 值/检验、Alpha 检验、独立的 PERMANOVA 和 dispersion、metadata/样本交集审计、指标定义及 R session。
- `diversity/LEVEL/figures/{pdf,png}` 与 `captions.tsv`：每个 Alpha、三种 PCoA、三个 Gamma 图；PNG 300 dpi。未知组保留描述性标识，图中样本数来自相同表。
- `state/*.json`、`execution.jsonl`、`logs/`：输入/输出 SHA256、配置和代码指纹、工具版本、真实 argv、资源、起止时间、PID、退出码及内容验证。`effective_config.json` 是实际解析配置。

各级以**到该 rank 的完整实际谱系**作为分类键，同名但不同谱系不合并。空 rank 保留空标签，另列 classification_status，保留实际上级谱系；不把未分类项丢掉或伪造 genus/species 标签。原始 feature 表留作追溯，不把不同研究的 OTU/ASV ID 当成同一 feature。原始 taxonomy 缺失项在 contract/映射列记录；源矩阵结构缺失、非数字、负值或非有限值报错，不填成零。shotgun 中明确缺席的 taxid 在完整样本报告间外连接为真实零；缺报告不补零。

聚合前后总量必须守恒。默认过滤为 prevalence≥0.10 且 mean relative≥1e-4（原标准化分析入口参数）。`relative.tsv` 分母为原表全部 features，含未分类；`filtered_original_denominator.tsv` 沿用此分母，列和可小于 1；`filtered_tss.tsv` 分母为过滤后保留的总量，列和为 1。二者同时保留并核验，不能用后者掩盖前者的损失。零文库被明确排除出可计算多样性，排除列表写入 contract/definitions。

规范入口默认多样性输入 `filtered_counts`，沿用原标准化流程；也可显式选择 counts/relative/filtered_tss。所有 level 分别计算，不把不同 level 指标混在一起。

| 指标 | 定义与输入 |
|---|---|
| Observed/richness | 所选过滤状态下严格 >0 的 taxon 数；是 rank richness，不称 ASV richness |
| Shannon | `-Σp ln(p)`，自然对数，p 为该样本所选列除以其总量 |
| Simpson | `1-Σp²`，不是 D 或倒数 |
| Chao1 | `S + f1(f1-1)/(2(f2+1))`；只对整数 counts 计算，relative/非整数估计量为 NA |
| Bray–Curtis | `Σ|x-y| / Σ(x+y)`；直接使用配置选定矩阵，不偷偷转成 TSS |
| Jaccard | 正值 presence/absence 的 `1-交集/并集` |
| Aitchison | 原实现 `ln(x+0.5)` 后按样本中心化，再算欧氏距离；0.5 加在所选矩阵单位上，counts 与 TSS 的意义不同 |
| PCoA | 完整特征分解；导出负特征值，图解释率分母为正特征值之和；没有静默距离校正 |
| Gamma | 同组样本列相加后计算 Observed/Shannon/1-D；counts 时深度加权，TSS 时每样本等权，All_samples 是描述汇总，不是新增组 |

本规范流程**没有稀释抽样**，没有深度/抽样次数；seed=1 是原置换函数的默认值（示例置换次数 99 为新验证配置，历史默认 999）。没有树依赖指标、UniFrac 或 NMDS 实现，不声称已计算它们，亦不虚造系统发育树。另一个旧 microeco 脚本虽存在，但未找到相应已运行产物，未为了增添指标把它当成成功路线。

原 Gamma 双组检验为保持组大小的独立标签置换，双侧绝对差、加一 p 值和 BH 校正。重复/blocked/含协变量设计下，通用 pooled Gamma 检验不具备现成正确估计量，输出 NA 原因，保留汇总。默认不为此另创新方法。

## 统计限制和修复

旧 `analysis_core.betadisper_test` 使用前两条 PCoA 轴和 Kruskal，**不是真正 PERMDISP**。精确旧函数仅留作审计，规范入口使用 vegan `betadisper(type="median")` + `permutest`，与 `adonis2(by="margin")` 分开导出。规范入口保存全部坐标/解释率，旧代码未导出的距离 TSV 已补全。这些统计结果不能冒充旧结果逐数值复现。

independent 设计的 Alpha 为双组 Mann–Whitney（R 渐近版本）或多组 Kruskal，BH 校正；其 p 值可能与旧 SciPy 自动 exact 选择不同。within_subject 且恰好两组、每人每组一个观测时用配对 Wilcoxon。其他重复/协变量 Alpha 设计明确 NA，保留描述。PERMANOVA 支持配置协变量并做完整案例/设计矩阵秩检查；分组和 block 缺失、组内少于 2、设计秩不足或无残差自由度时记录原因。

已知 subject 重复不能忽略：manifest/常见 metadata subject 字段发现重复且未配置时预检报错；配置后仍无正确设计时只做描述及 NA 检验。within_subject 要求 subject 与置换 block 一致且各人的组别可交换；组别完全属于个体（例如病例/对照且每人多访视）不能做组内 label 置换，会报告 group_invariant_within_subject。应由研究方案明确选标本或另设适合的研究级模型，不能把这些重复 run 当独立样本。源码不能识别未提供的个体映射，必须根据原 metadata 补充。

单臂/未知饮酒组别不取消 Alpha、距离、坐标、Gamma 描述。每个检验的 NA 原因及实际完整案例样本数均导出；Unknown 从不视为 Control。

## 原有合并矩阵支路

常规入口不做自动跨研究合并。原有 `integrated_diversity.R` 仅接收已有 integration root（包含 FINAL 或 02_merged_cohorts 的约定 Genus 矩阵/metadata），保留原过滤前 Alpha、过滤后 Beta、Dataset 调整和 seed=20261005/999 permutations。它不是规范入口的 counts 参数模式，不能把其数值当成同一种输入处理。

```bash
AREK_POOLED_SCOPE=verified_same_host_site_pipeline \
Rscript legacy/pooled_diversity/integrated_diversity.R /path/to/new_diversity_output /path/to/existing_integration_root
```

仅在核实同宿主、同部位、可比处理方法、真正 Genus 分类身份及明确纳入范围后使用该单独入口。新增的缺失值/重复个体保护会拒绝不合适输入。原代码根据“count-like”数值四舍五入算 Chao1 不可靠；整理版默认不重建 Chao1，只有用户明确设 `AREK_INTEGRATION_INTEGER_COUNTS=1` 且输入本身为整数 counts 时计算，不从 relative 伪造。完整合并队列与来源矩阵构建未在本轮重跑；建议新数据使用按研究的规范入口。

保留的 legacy shell 原生支路须设置真实 RAW_DATA_ROOT、RESULT_ROOT、DATABASE_ROOT 等，并显式 `AREK_LEGACY_SOURCE_MODE=isolated_research`；该模式用于已知历史批处理/诊断的独立目录，非规范 resume 保证。不要把它们的文件存在判断当成本包通过验证。SOLiD 还需独立 SOLID_ENV/BFAST_HOME/BFAST_DB_ROOT，源码与 BFAST 原始源归档已收齐，详见参考库说明。

## 验证与常见错误

原始 ZIP 的 `validation/VALIDATION_zh.md`、`validation/final_tests/` 记录源包历史验证，未上传代码仓库。真实示例来自 PRJDB10075 原始导入 QZA 内的原 FASTQ，各 run 前 2,000 对 reads；公开 accession 与示例文件 hash 见 `examples/provenance.json`。没有使用合成 reads 充当生物验证；metadata 模板和统计故障夹具中的分组仅为软件测试。当前仓库检查见 `validation/README.md`。

```bash
# 可在新解压目录复做完整测试（需要分类器与安装环境）
python validation/validate_release.py --work-root /tmp/arek_recheck_full \
  --classifier "$AREK_CLASSIFIER"
# 无参考库时，真实现成矩阵及故障/统计设计验证：
python validation/validate_release.py --work-root /tmp/arek_recheck_matrix --matrix-only
```

验证包含 help、所有配置解析、真实扫描配对、多 lane/孤立 mate 冲突、dry-run 无写入、预检依赖/缺库、首次运行、原样 resume、删除 Genus 距离只恢复其多样性/统计/图分支、坏样本与独立任务继续、真实 SE、矩阵/相对输入、重复个体保护与独立数值核对。完整 ZIP 会另目录解压后校验 CRC、payload hash、内部路径并运行代表入口，记录包外解压报告。

常见错误及处理：

- `classifier QZA missing`：用 `--classifier` 指向已核验文件；不是旧服务器绝对路径或联网占位下载。
- `missing executable/R dependencies`：激活对应 explicit 环境；不要为了避错跨用旧环境。
- `orphan/duplicate/multiple lanes`：根据真实 run/lane/specimen 表确认映射，写显式 manifest；不要按相似文件名猜。
- `mate ID mismatch/invalid FASTQ record`：定位具体文件和记录；修复副本/重新获取损坏输入后 resume，不修改原始数据。
- `color-space reads detected`：SOLiD 不能送入 base-space 路线；查独立诊断支路，不能把数字序列当 DNA。
- `relative input does not sum to 1`：核实是百分比、截断类群还是不同分母，明确单位和缺项后再接入，不能自动重新归一化掩盖问题。
- `invalid or stale`：查看 state/日志中缺失输出或指纹变化，使用 resume 恢复受影响分支；`--from-stage` 不会绕过上游核验。
- `already running`：另一个本地任务持有该数据集锁；查 running JSON 中 PID 和子进程，不重复提交。
- `not_applicable`：源分类没有该 rank 的命名结果，是适用性结果，不是 Species 完成。
- 资源不足：按 MemAvailable 预检，降低验证线程/调整资源请求、选择独立输出，勿停止其他作业；大数据库也需要足够磁盘。

`FILE_MANIFEST.tsv` 列当前模块文件的相对路径、bytes 和 SHA256，不包含自身。原 ZIP 的完整清单不上传；保留源码的去向和源包 SHA256 见 `docs/provenance/reorganization.json`。`docs`/`validation` 中旧服务器路径仅为不可执行的来源证据，规范代码/示例配置不依赖旧脚本或旧账户目录。没有原始全量 FASTQ、完整历史环境或大型参考库随包发送；不会把未测试支路写成“全通过”。
