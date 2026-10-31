# 来源与仓库整理记录

权威输入为用户上传的 `AREK_FASTQ_to_Abundance_Diversity_pipeline_20261030.zip`。整理日期：2026-10-31（America/Chicago）。

| 文件 | 用途 |
| --- | --- |
| [reorganization.json](reorganization.json) | 原 ZIP SHA256、来源文件数、整理范围与代码变动说明 |
| [file_relocation.tsv](file_relocation.tsv) | 本仓库保留的源文件去向、原 hash 与原样 / 适配标记 |
| 模块根目录 `FILE_MANIFEST.tsv` | 当前整理后模块的文件与 hash；不包含清单自身 |

正式计算源码、配置、环境文件和示例数据与上传包逐字节一致。移动六个独立历史入口到 `legacy/`，调整三个验证程序的证据输出目录，更新使用说明与源码路径映射，并从说明清单中移除历史服务器路径。

全部历史执行日志、历史验证快照、队列 / 结果清单、项目扫描、QZA 中间产物、测试输出表和图均未纳入代码仓库；它们保留在用户原始 ZIP 中。仓库只保留源码、配置、依赖清单、说明和公开小样本测试夹具，以及当前整理检查的汇总。

历史函数名称中的 `historical_*` 有些仍属于正式算法依赖，例如 Bracken 谱系解析和多样性数学实现；这些文件保留在 `scripts/`，与可独立执行的 `legacy/` 支路区分。
