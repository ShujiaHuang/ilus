## Update

- [2023-11-01] 添加基于 SOAPnuke 的原始数据过滤模块
- [2023-10-31] 完成 Sentieon 模块的添加，并且流程中不需要有额外的 CombineGVCF 这个步骤。
- [2023-10-31] 为 genotypeGVCFs 添加 --dbsnp 参数，用于注释已知位点
- [2023-10-30] 比对统计信息的计算和数据污染计算从 BQSR 前移到 markduplicate 这一步之后。并将 `.recal.table` 后缀改为 `.BQSR.recal.table`
- [2023-10-25] 删除配置文件中 vqsr_options 中的 `--max-gaussians`，改为 SNP(`vqsr_snp_max_gaussians`) 和 Indel (`vqsr_indel_max_gaussians`)单独设置.
- [2023-10-24] 添加 sentieon 模块，并基于 sentieon 构造 WGS/WES 分析流程.
- [2023-10-22] 修改配置文件的格式，将 `resources` 从 `gatk` 拎出来，单独成一块，方便 sentieon 和 GATK 共同调用.  
