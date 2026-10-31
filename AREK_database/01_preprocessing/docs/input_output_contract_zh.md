# 输入输出与后续模块接口

## 输入

| 文件 / 参数 | 最低要求 |
| --- | --- |
| 配置 JSON | 以 `config/historical.json` 为默认值合并；未知键报错；路径以配置文件所在目录为基准 |
| FASTQ 样本表 | TSV，必需 `dataset_id, sample_id, run_id, layout, r1, r2`；相对读段路径以样本表位置为基准；SE 的 r2 为空，PE 必须成对 |
| 生物身份 | `specimen_id`、`subject_id` 仅在有证据时填写；多 lane / 同标本多 run 需明确上游处理，不自动合并 |
| Metadata | TSV，首列唯一 sample_id，其他列由配置指定；未知组别保留缺失，不能归入 Control |
| Counts 矩阵 | feature × sample，行列 ID 唯一，非负有限数值；结构缺失不补零；非整数估计 counts 的 Chao1 为 NA |
| Relative 矩阵 | 明确 0–1 单位，各列和为 1 的完整分母；百分比 / 截断矩阵应在上游显式处理 |
| Taxonomy | 唯一 feature ID 与 `Taxon` 列；仅本级矩阵可改用 `native_level`，不编造缺失的其他 rank |
| 外部参考库 | 扩增子需 SILVA 分类器；shotgun 需 Kraken / Bracken 库，按配置需 Bowtie2 宿主索引；详见 databases_zh.md |

dataset_id 是明确的研究 / 宿主 / 部位 / 处理路线分析单元。同一研究不同部位应登记为不同单元。扫描只提出候选文件映射，不能替代研究 metadata 确认。

## 标准输出

所有正式结果写入 `--out/DATASET_ID/`。

| 相对路径 | 内容 / 单位 | 后续用途 |
| --- | --- | --- |
| `upstream/` 或矩阵路线的 `input/` | QZA / 原始 feature / taxonomy，或接入矩阵副本 | 来源追溯；不把跨研究 OTU ID 当同一 feature |
| `abundance/LEVEL/counts.tsv` | 该 rank 的 counts；Bracken 是估计 reads | 需要计数的整合 / 分析；relative-only 不生成此文件 |
| `abundance/LEVEL/relative.tsv` | 原表全部 features 为分母的 0–1 丰度 | 完整丰度谱 |
| `abundance/LEVEL/filtered_counts.tsv` | 按配置过滤后的 counts | 默认多样性输入；relative-only 不生成此文件 |
| `abundance/LEVEL/filtered_original_denominator.tsv` | 过滤后保留类群，分母仍为原表总量 | 核对过滤造成的丰度损失；列和可小于 1 |
| `abundance/LEVEL/filtered_tss.tsv` | 过滤后重新归一化的 0–1 丰度 | 明确需要此分母的下游分析 |
| `abundance/LEVEL/taxonomy.tsv`、`feature_to_taxon.tsv` | 完整谱系分类键、源 feature 映射 | 以真实分类身份对齐 |
| `abundance/LEVEL/contract.json`、`mass_balance.tsv`、`filter_audit.tsv` | 单位、总量、过滤与分类审计 | 后续读表前核验 |
| `diversity/LEVEL/` | Alpha / Gamma，三种 distance / coordinates / eigenvalues，统计与样本审计 | 数据库结果展示与分析 |
| `diversity/LEVEL/figures/{pdf,png}/`、`captions.tsv` | 图与来源一致的图注 | 展示 |
| `effective_config.json`、`state/*.json`、`execution.jsonl`、`logs/` | 参数、输入输出 hash、源码 / 工具指纹、命令与退出状态 | 执行与复现记录 |

某 rank 没有命名分类时输出 `not_applicable.json`，不应当作完成了该 rank 的丰度与多样性。零文库按现有契约排除出多样性，并记录原因。

后续模块应同时读取数值表与 contract / taxonomy / metadata，明确使用原分母还是 filtered TSS。不能从 relative 反推 counts，不能静默把缺报告当作零丰度，不能自动推断饮酒 / 疾病标签。
