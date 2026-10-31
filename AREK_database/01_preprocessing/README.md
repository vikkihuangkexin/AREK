# 01 · 预处理：FASTQ → 分类丰度 → 多样性

AREK 数据库的第一个代码模块。处理明确登记的数据集，生成可追溯的分类层级丰度表、多样性、统计表与图。整理保留原计算方法；正式计算源码与上传包逐字节一致。

## 目录

| 目录 / 文件 | 内容 |
| --- | --- |
| `run_pipeline.py` | 统一命令入口 |
| `workflow/` | 配置、批次调度、路由、内容 checkpoint 与输入检查 |
| `scripts/preprocess/` | 正式流程的去宿主辅助脚本 |
| `scripts/abundance/` | Bracken 合并、分类层级聚合与丰度过滤 |
| `scripts/diversity/` | Alpha / Beta / Gamma、组间统计与绘图 |
| `config/` | 方法默认值、各路线示例与指标契约 |
| `envs/` | 各路线独立的 Linux conda 锁定文件与版本记录 |
| `examples/` | 样本表、metadata 模板、真实小样本 FASTQ 和矩阵 |
| `validation/` | 可重新运行的核验程序；结果写入指定工作目录 |
| `legacy/` | 历史扩增子、SOLiD 和合并矩阵多样性支路 |
| `docs/` | 步骤、参数、参考库、数据契约和源码来源索引 |

## 处理步骤

| 顺序 | 操作 | CLI / 阶段名 | 主要结果 |
| --- | --- | --- | --- |
| 00 | 扫描并确认样本映射 | `scan` | 候选样本表、配对冲突；需人工确认 |
| 01 | 预检并查看执行计划 | `check`、`dry-run` | 输入、依赖与数据库检查 / 真实命令计划 |
| 02 | 读段处理 | `preprocess` | 扩增子导入 / PE 拼接 / 去重复 / 聚类 / 去嵌合；或 shotgun 清洗 / 去宿主 / 分类报告 |
| 03 | feature 表生成 | `features` | OTU feature 表，或各 rank 独立 Bracken 表，或接入现成矩阵 |
| 04 | taxonomy 注释 | `taxonomy` | 扩增子分类注释；shotgun 使用其已有分类身份 |
| 05 | 分类层级丰度 | `levels` | counts / relative、过滤前后表、映射和质量守恒审计 |
| 06 | 多样性与统计 | `diversity` | Alpha、距离 / PCoA、Gamma、PERMANOVA、PERMDISP 与样本审计 |
| 07 | 出图 | `plot` | PDF / PNG 与图注 |
| 08 | 结果核验 | `validate` | 内容、参数、源码和工具指纹核验 |

这里的 `preprocess/features/taxonomy/levels/diversity/plot` 是已有 CLI 的六个阶段选择器，具体外部工具子步骤由路由决定，见 [完整步骤对应表](docs/steps_zh.md)。

## 选择路线

| 输入 | `route` | 方法 / 适用性 |
| --- | --- | --- |
| 扩增子 PE / SE FASTQ | `amplicon_vsearch97` | QIIME2 + vsearch 97% OTU + SILVA sklearn；不是 ASV 去噪流程 |
| Shotgun FASTQ | `shotgun_kraken_bracken` | fastp → 可配置去宿主 → Kraken2 → 各 rank 独立 Bracken |
| 现成矩阵 | `matrix` | feature + taxonomy，或明确 `native_level` 的本级矩阵 |

SOLiD color-space 仅保留在独立历史支路，尚无完整分类链成功证据。Shotgun 全链和完整合并队列未在此次整理中重跑。原包的扩增子小样本验证与当前检查范围见 [validation/README.md](validation/README.md)。

## 快速开始

以下命令从仓库根目录开始。环境锁定文件支持 Linux x86_64；参考库需另行准备。

```bash
cd AREK_database/01_preprocessing

# 已有矩阵 / 多样性路线
conda create -n arek_analysis --file envs/analysis.linux-64.explicit.txt
conda activate arek_analysis
python run_pipeline.py --help
python run_pipeline.py check --config config/matrix_example.json --dataset REAL_MATRIX
python run_pipeline.py dry-run --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py run --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
python run_pipeline.py validate --config config/matrix_example.json \
  --dataset REAL_MATRIX --out /tmp/arek_matrix
```

扩增子使用 `envs/amplicon.linux-64.explicit.txt` 环境及外部分类器：

```bash
export AREK_CLASSIFIER=/path/to/silva-138-99-nb-classifier.qza
python run_pipeline.py check --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER"
python run_pipeline.py run --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_pe
```

原样续跑使用 `resume`；指定阶段使用 `--from-stage` / `--to-stage`，跳过的上游仍需有有效 checkpoint。只有外部矩阵时应选择 `matrix` 入口。

```bash
python run_pipeline.py resume --config config/real_pe_example.json \
  --dataset PRJDB10075 --classifier "$AREK_CLASSIFIER" --out /tmp/arek_pe \
  --from-stage levels --to-stage plot
```

样本表内相对 FASTQ 路径以样本表为基准；配置内相对路径以配置文件为基准。样本身份与输入要求见 [input_output_contract_zh.md](docs/input_output_contract_zh.md)。

## 文档与边界

- [README_zh.md](README_zh.md)：完整安装、命令、指标定义和限制。
- [docs/steps_zh.md](docs/steps_zh.md)：每一步的入口、源码、依赖与产物。
- [docs/parameters_zh.md](docs/parameters_zh.md) / [docs/databases_zh.md](docs/databases_zh.md)：方法参数与外部参考库。
- [validation/README.md](validation/README.md)：本次检查和历史验证的区别。
- [legacy/README.md](legacy/README.md)：历史支路的用途与限制。
- [docs/provenance/README.md](docs/provenance/README.md)：保留源码的去向与原包 SHA256。

本模块不执行跨研究合并，不推断饮酒 / 疾病组别，不包含差异菌、功能预测、分型、预测模型或数据库服务部署。完整研究数据与大型参考库放在仓库外，`examples/` 仅用于验证。
