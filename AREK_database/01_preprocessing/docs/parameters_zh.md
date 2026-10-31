# 参数来源与方法差异

`config/historical.json` 以可解析 JSON 保存规范缺省值。方法阈值与资源设置的证据不同，不把便携配置称为旧运行逐参数重放。

| 参数 | 值 | 来源/状态 |
|---|---|---|
| OTU identity | 0.97 | v5 shell + PRJDB10075 clustered QZA provenance |
| median read length | ≥100 bp | v5 样本布局判断；旧探测前3000 reads，本规范 prepare 检查全文件；明确选择 mate，不再自动挑长端 |
| PE minovlen / maxdiffs / minlen | 10 / 10 / 1 | 实际 merged QZA action；其 truncqual/maxns/maxee/minmergelen/maxmergelen 为 null，allowmergestagger=false |
| derep min_seq_length / min_unique_size | 1 / 1 | 实际 derep QZA，derep_prefix=false |
| UCHIME dn/mindiffs/mindiv/minh/xn | 1.4 / 3 / 0.8 / 0.28 / 8 | 实际 nonchimeras QZA |
| 分类 | sklearn, confidence=0.7, read_orientation=auto | taxonomy QZA；16 jobs 为旧记录，2 为本轮验证资源 |
| SILVA 训练 | 138 99%；2024.5 q2-feature-classifier；7-mer，8192 features，alpha=.001，fit_prior=false | 分类器自身 provenance；完整训练参数见 classifier_training_provenance.json。扩增区域/引物未从已恢复记录证实，不能称 V3–V4 专用分类器 |
| fastp v2 PE | min length50；correction=false；dup eval关闭（工具支持时）；compression1 | PRJNA613834 v2 源码；project_run_parameters.tsv 记录实际线程8、fastp4、read_len150等，未把源码16当实际8 |
| fastp original PE | min length50；detect_adapter_for_pe；correction=true | PRJCA021196 原始源码；original SE 无 correction 和 PE adapter flag |
| fastp 质量参数 | 使用已锁定工具缺省 | 原命令未显式记录 quality cutoff 等；并无可核对的历史逐项配置快照。版本锁保存实际工具版本，但不能伪称历史所有默认均恢复 |
| Bowtie2 PE | very-sensitive；samtools -f12 -F256；name sort；fastq -n | original/v2 源码；仅保留双端均未比对且非secondary；整理后 nonhost 可为未压缩 FASTQ，不改变 read 内容 |
| Bowtie2 SE | very-sensitive；samtools -f4 -F256 | original 源码；v2 不支持 SE |
| Kraken2 | paired/use-names；confidence缺省0 | 原命令未给 confidence；0 为锁定工具标准缺省，未恢复独立历史自定义值 |
| Bracken | r150，t10，P/F/G/S分别估计 | 原 shell 与真实日志；不是先 Species 再向上求和，rank 间不强制总量一致 |
| filtering | prevalence .10；mean relative 1e-4 | 02_run_taxonomy_analysis_from_standardized.py 参数。旧同名不同谱系合并已修复，可能改变 taxa 数，不称旧矩阵精确复现 |
| TSS | 原始全 feature 分母；过滤后另输出 retained 分母 | 不静默改变分母，二者同时保存 |
| Alpha | Observed；自然log Shannon；1-D Simpson；bias-corrected Chao1 | analysis_core 精确函数提取；Chao1 只允许整数counts，非整数估计或相对输入为 NA |
| Beta | counts直接 Bray；二元Jaccard；ln(x+.5) CLR欧氏 | analysis_core；可显式切换输入矩阵，不能混淆 counts 与 TSS 的伪计数单位 |
| Gamma | 同组列相加后 S/H/1-D；独立标签置换加一p | run_manifest_microbiome_analysis.py 精确函数提取。重复/blocked/adjusted设计不作不当独立置换 |
| 置换 | 999，seed1 | 历史 Python 函数默认。验证99/19是明确的新小测试配置 |
| 稀释 | 无 | 规范 Python 多样性实现没有 rarefaction，不加造深度和抽样次数 |
| phylogeny/NMDS | 无 | 所选择的已运行标准化路线没有相关实现；没有虚造树 |
| 原合并矩阵支路 | seed20261005，999置换；Alpha未过滤Genus relative、Beta过滤TSS；Aitchison实现见原R | 与规范 counts 入口独立；既有成功结果，不将其指标混称为同一配置 |
| threads/memory | 缺省2/8 GiB，matrix示例1 GiB | 便携验证配置；内存为预检估计，不是历史值/实测峰值 |

规范入口的错误定位、谱系聚合、缺组别描述性输出、内容 checkpoint、距离/坐标/eigenvalue导出、真实 PERMDISP、重复个体保护属于整理版修复。旧指标函数保留在 `historical_*` 供数值对照，错误的旧 betadisper 不被规范入口调用。

本包没有为未恢复的引物、区域、read-to-subject 映射和数据库 build 名称填默认。原脚本自动取长端的决策可以通过明确 SE 样本表重现，但必须有相应 read_layout/fastq_used 证据；扫描不会根据长度自动冒充 biological mate。

规范 shotgun 入口保留 clean/nonhost 等已声明中间产物以供内容续跑核验，与旧 v2 的 CLEAN_INTERMEDIATE=1 操作策略不同；不改变 FASTQ 内容或分类算法，但全量磁盘需求应相应估计。原 cleanup/FORCE 控制不被冒充成生物学参数。

历史caller还有 min-samples=6、min-group-size=3（02标准化入口）。规范入口为满足无组/小样本描述保留输出，不沿用“少于6样本整个模块取消”的策略；修正统计入口目前要求每组至少2个有效观测，是新入口的可估计性下限，不能冒充旧caller的3样本纳入规则。要严格保留某研究旧统计纳入范围，应提供该研究已核实的纳入metadata并核对组数，不把本包4样本软件统计测试解释为旧研究分析重现。该差异与R渐近p值/PERMDISP修复一样，属于明确登记的新入口行为。
