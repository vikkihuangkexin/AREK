# 环境清单

- analysis：矩阵处理/规范多样性/作图所需闭包。
- amplicon：规范 QIIME2/vsearch FASTQ至多样性所需闭包。
- shotgun：规范 fastp/Bowtie2/Kraken2/Bracken至多样性所需闭包。
- integration：独立既有 Genus 合并矩阵 R 路线所需闭包，含ggplot2。
- solid_tools：原服务器 SOLiD 配套Python/samtools/make/zlib工具闭包；BFAST另从随包源码构建。

explicit文件均由实际conda-meta种子+传递依赖生成，Linux x86_64；versions.json 列实际版本、build、包SHA256/MD5。缺失依赖记录数组为空。它们锁定的是实际可用安装，不代表本轮在干净机器重新下载安装过。网络/远端包失效应明确报错，不能自动切换新版冒充同环境。

仍需 Linux 系统提供 bash、GNU coreutils、awk/sed/grep、gzip、tar、time；编译 BFAST 需系统C编译器和zlib开发头。`system_tools.tsv` 记录当前服务器工具版本，不把系统工具打进代码ZIP。BFAST源码包含自身第三方授权文本，重新编译未在本轮验证。

不要将 BFAST 的旧工具要求强装进 QIIME2 环境。使用新安装位置激活对应环境即可，代码不依赖原服务器环境prefix。环境/参考库和路径都会参与实际checkpoint指纹，切换环境后应重新验证/运行。
