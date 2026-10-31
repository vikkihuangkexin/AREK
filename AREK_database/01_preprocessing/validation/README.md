# 验证入口与证据范围

验证程序放在本目录，生成结果写入各命令的工作目录。历史日志、结果快照和完整验证输出仅保留在原始 ZIP，不上传代码仓库。

## 可重新执行的检查

从 `AREK_database/01_preprocessing/` 执行：

```bash
# 当前模块文件清单、源码语法、内部路径与导入静态核验
python validation/audit_package.py .

# 无需 R / QIIME2 的目录搬迁、各路线 dry-run 与真实矩阵数值检查
python validation/check_repository.py --report /tmp/arek_repository_checks.json

# 真实现成矩阵 / 故障 / 统计设计验证：需 analysis 环境及 R vegan / permute
python validation/validate_release.py --work-root /tmp/arek_validation_matrix --matrix-only

# 包含真实 PE / 显式 R1-SE 小样本：另需 amplicon 环境与外部分类器
python validation/validate_release.py --work-root /tmp/arek_validation_full \
  --classifier "$AREK_CLASSIFIER"

# 补充契约与统计设计检查，仍需相应运行环境
python validation/check_contract_transitions.py /tmp/arek_contract
python validation/check_extra_statistics.py /tmp/arek_statistics

# Linux 私有进程锁测试；仅终止它自己创建的测试进程
python validation/check_orphan_lock.py /tmp/arek_lock
```

执行后，validate_release 的证据默认位于 `--work-root/evidence/`；contract / extra statistics 也写入各自工作目录的 `evidence/`；锁测试结果写入其工作目录的 `orphan_lock_test.json`。不覆盖源包历史验证记录。

## 历史验证

源包记录过真实 PRJDB10075 小样本的扩增子 PE 完整链、显式选择 R1 的 SE 链、真实矩阵与相对丰度入口、续跑 / 故障恢复、数值核验和统计设计检查。详细来源与未验证范围需查原始 ZIP 的 `validation/VALIDATION_zh.md`，不将其视为当前仓库重新执行的证据。仓库中的 [公开样例来源](../examples/README.md) 仅保留 accession、提取方法和文件 hash，不包含队列表型。

## 本次仓库整理

本次检查与原包历史验证分开登记于 [repository_checks_20261031.json](repository_checks_20261031.json)。重点是目录移动后的入口、内部路径、配置、扫描、矩阵计算和文件完整性；不以此宣称重新跑过 FASTQ 全链。当前运行环境没有 QIIME2、Rscript 或大型参考库，完整 FASTQ 及 R 统计链没有在此环境重跑。
