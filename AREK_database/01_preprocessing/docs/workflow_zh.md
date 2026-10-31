# 规范调用关系与状态契约

```text
run_pipeline.py
  workflow/cli.py             参数/配置/批次/metadata/阶段范围
    workflow/scanner.py       命名扫描、FASTQ格式/配对检查
    workflow/common.py        flock外围、执行日志、内容/指纹checkpoint
    workflow/routes.py
      amplicon()              QIIME2/vsearch原始方法命令
      shotgun()               fastp/host_filter.sh/Kraken2/Bracken
        scripts/abundance/historical_bracken_merge.py (恢复的谱系解析函数，rank修复)
    scripts/abundance/levels.py
    scripts/diversity/compute.py
      historical_metrics.py   原数学实现
      historical_gamma.py     原Gamma定义与独立置换
    scripts/diversity/statistics.R  修正的统计入口
    scripts/diversity/plot.py       TSV独立重画
```

原始脚本 SHA256 与当前入口见 source_script_map.tsv；历史执行日志和来源路径仅保留在原始 ZIP。规范入口所需的 Python/R/shell 模块都在包内，不调用服务器旧源码目录。外部依赖仅为指定环境里的命令行工具和用户传入数据/参考库，工具用 explicit 锁定文件恢复。

扩增子 DAG：prepare → import → (PE merge) → dereplicate → cluster → chimera → filter_table/filter_seqs → features；filter_seqs → classify → taxonomy；features+taxonomy → 各 level → 各 level metadata → diversity → statistics → plot。

shotgun DAG：每样本 fastp → host（可显式按历史配置跳过）→ kraken → 各 rank Bracken；同 rank 报告合并 → 对应 native level。P/F/G/S 独立估计，不把 S 向上求和冒充 F/G 的 Bracken结果；各 rank 内再执行相同多样性模块。

矩阵 DAG：严格读表/保留原feature → taxonomy或明确native_level → level/diversity。只给Genus表不会生成Phylum、Species。

所有阶段成功需：输入/方法参数/所用代码/工具版本匹配；输出非空、结构/数值/ID核验通过；每个声明输出 hash 一致。QZA 检查ZIP CRC/data结构，feature表与代表序列ID、taxonomy ID匹配。空目录、退出0、单独done文件都不够。长任务记录父/子 PID，规范入口同步等待退出；同数据集外层 flock 不允许重复启动。

resume 使用本阶段相关源码指纹，图代码改变不影响FASTQ阶段。通用调度/读表代码改变可能影响所有阶段，会保守作废。某产物损坏后，该阶段及显式依赖下游重做，其他 level 保持复用。本轮删除 Genus 距离的日志验证只重做 Genus diversity/statistics/plot。

`--from-stage` 前的上游不重新执行，但会核对 checkpoint；不能拿陌生旧输出直接跳过指纹。外部既有表请从 matrix 路线正式接入。`validate` 自身不写产物、不提交任务；检查整个选定DAG（含上游）内容。`check`预检错误会保留可定位信息；批量逐数据集捕获错误，末尾返回2，而不是中途把后续数据集全部取消。

日志中 effective_configuration_sha256 是完整实际配置指纹，本阶段 fingerprint 使用相关参数/输入内容/代码/工具/资源，不因完全无关绘图参数误使上游过期。工具版本来自实际 Python包、conda-meta、可执行路径及R session；不是只写一个预设版本字符串。

既有 merged Genus R 和 legacy shell 是单独登记的历史支路，不参与上述 checkpoint。它们的限制与未验证范围在 README 明确，不把原生脚本的exists判断误称为规范恢复保证。
