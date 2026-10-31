# 独立历史支路

这些代码保留来源和特定研究处理能力，未接入 `run_pipeline.py` 的完整 checkpoint 契约。正式常规分析使用模块根目录的统一入口。

| 目录 | 原脚本 / 用途 | 当前限制 |
| --- | --- | --- |
| `amplicon/selected_se_batch.sh` | 原批量 SE derep 路线 | 隔离研究入口；未在此次整理执行 |
| `amplicon/PRJCA022064.sh` | 特殊命名与 best-mate 支路 | 需确认该研究的输入设计；未执行 |
| `solid/` | SOLiD color-space 准备、BFAST 索引与分块去宿主 | 部分路线，缺少后续完整分类成功证据 |
| `pooled_diversity/integrated_diversity.R` | 对既有 merged Genus 矩阵计算多样性 | 不负责构建 / 合并矩阵；完整队列未重跑 |

历史 shell 入口保留原有的显式环境开关 `AREK_LEGACY_SOURCE_MODE=isolated_research`，并要求真实 RAW_DATA_ROOT、RESULT_ROOT、DATABASE_ROOT 等目录。SOLiD 还需单独的 BFAST 环境和参考索引。参数及来源见模块根目录的 [README_zh.md](../README_zh.md)、[docs/source_script_map.tsv](../docs/source_script_map.tsv)。

合并矩阵入口的调用示例（从模块根目录执行）：

```bash
AREK_POOLED_SCOPE=verified_same_host_site_pipeline \
Rscript legacy/pooled_diversity/integrated_diversity.R \
  /path/to/new_diversity_output /path/to/existing_integration_root
```

调用前需确认同宿主、同部位、可比处理方法、真实 Genus 身份和纳入范围。它保留历史的过滤前 Alpha / 过滤后 Beta 约定，与正式入口默认输入不同；不能把两者视为同一分析配置。
