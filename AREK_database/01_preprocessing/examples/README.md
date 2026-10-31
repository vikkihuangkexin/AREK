# 公开小样本与软件模板

`real_fastq/PRJDB10075/` 包含公开 accession PRJDB10075 下 DRR234531–DRR234534 的小样本，每个 run 保留前 2,000 对 reads。这些读段由原包从已导入的原 FASTQ 归档恢复；没有合成或更改序列。来源标识、提取方式和文件 SHA256 见 [provenance.json](provenance.json)。

`matrix/counts.tsv` 和 `matrix/taxonomy.tsv` 是这套小样本的 feature / taxonomy 验证夹具，不包含全研究队列或真实饮酒 / 疾病 / 临床分组。`real_samples.tsv` 明确 PE 映射，`real_se_samples.tsv` 显式选择相同数据的 R1，用于检验 SE 入口。

`metadata_template.tsv` 的 RUN_A / RUN_B、组别、年龄和性别均为字段格式示意，不是真实研究参与者记录。Shotgun 样本表中的路径是待用户填写的示例，不包含真实 shotgun 数据。

这些文件仅用于安装和入口验证，不能据此给出酒精相关研究的生物学结论。完整研究数据和大型参考库不随代码上传。
