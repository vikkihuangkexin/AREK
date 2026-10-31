# 预处理步骤与源码对应

从 `AREK_database/01_preprocessing/` 执行命令。生产运行均通过 `run_pipeline.py`；不要单独运行依赖内部函数的 Python 模块。

## 步骤表

| 顺序 | 输入 → 操作 → 输出 | 调用入口 / 源码 | 阶段范围选择器 |
| --- | --- | --- | --- |
| 00 | 原始目录 → 明确命名配对扫描 → files / samples / conflicts TSV | `scan`；`workflow/scanner.py` | 独立命令 |
| 01 | 配置、样本表、参考库 → 预检与计划 → 检查记录 / 命令计划 | `check`、`dry-run`；`workflow/cli.py` | 独立命令 |
| 02 | FASTQ → 研究路线对应的读段处理与分类准备 → QZA 或各样本分类报告 | `workflow/routes.py`；shotgun 调用 `scripts/preprocess/host_filter.sh` | `preprocess` |
| 03 | 扩增子过滤后的 table / 代表序列 → feature TSV / BIOM / FASTA；或各样本 Bracken → rank 矩阵；或复制外部矩阵 | `workflow/routes.py`、`workflow/cli.py`；`scripts/abundance/historical_bracken_merge.py` 的谱系解析函数 | `features` |
| 04 | 扩增子代表序列 + 分类器 → taxonomy QZA / TSV | `workflow/routes.py` | `taxonomy` |
| 05 | feature 表 + taxonomy，或 native-level 表 → 层级聚合、过滤、分母审计 | `scripts/abundance/levels.py` | `levels` |
| 06 | 选定层级矩阵 + metadata → Alpha、距离、PCoA、Gamma、组间统计 | `workflow/cli.py`；`scripts/diversity/compute.py`、`historical_metrics.py`、`historical_gamma.py`、`statistics.R` | `diversity`，含 metadata 与 statistics 子阶段 |
| 07 | 多样性 TSV → PDF / PNG / captions | `scripts/diversity/plot.py` | `plot`；另有独立 `redraw` |
| 08 | 已有结果 + 原输入 / 配置 / 工具 → 核验成功与否 | `validate`；`workflow/common.py`、`workflow/cli.py` | 独立命令，可限定阶段范围 |

`workflow/common.py` 负责各阶段的命令执行、日志、内容指纹与 checkpoint；不属于单独的生物计算步骤。`run` 重做选定阶段，`resume` 按指纹复用，`validate` 只核验。

## 三条正式路线

**扩增子：**明确 PE / SE → QIIME2 导入 → PE merge（SE 跳过）→ dereplicate → 97% OTU cluster → UCHIME → 过滤 table 与代表序列 → feature 导出 → SILVA sklearn → taxonomy → 分类丰度 → 多样性 / 统计 → 图。现有 CLI 把 import 至嵌合过滤归入 `preprocess`，不要把表中的六个选择器误当成六条外部命令。

**Shotgun：**fastp → 去宿主（依配置）→ Kraken2 → 各 rank 的 Bracken → 同 rank 合并 → 分类丰度 → 多样性 / 统计 → 图。现有 CLI 把样本级外部处理归入 `preprocess`，rank 合并归入 `features`，没有另行训练或调用 SILVA。P / F / G / S 各自估计，不从 species 求和冒充 genus Bracken 结果。`variant=v2` 要求 PE；`original` 支持 PE / SE。

**矩阵：**严格读入 → 保存原 feature 与 taxonomy → 层级聚合或 `native_level` → 多样性 / 统计 → 图。没有原始 FASTQ 时不宣称重现了 FASTQ 前处理；只有 Genus 时不生成缺失的 Phylum / Species。

## 数据库模块之间的接口

后续模块读取 `--out/DATASET_ID/` 下的已验证产物，至少记录 dataset_id、宿主 / 部位、分类层级、输入单位、过滤规则、样本身份和方法来源。正式接口见 [input_output_contract_zh.md](input_output_contract_zh.md)。

每个数据集单独处理；不同研究的 feature ID 不能直接合并。若后续统一 Genus 矩阵，应在新的整合模块中明确同名分类身份、宿主、部位、处理方法与原分母，而不在本阶段隐式完成。

本说明的机器可读索引为 [step_registry.tsv](step_registry.tsv)。历史单独入口见 [../legacy/README.md](../legacy/README.md)。
