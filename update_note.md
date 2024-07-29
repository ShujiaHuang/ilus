## Update


- [2023-12-05] 版本 v2.0.0 发布
- [2023-12-05] WES 改名为 capseq，反正 WES 也是 capseq 的一种，这样改好听
- [2023-11-23] 改写了 `genotype-joint-calling` 模块，输入参数 `gvcf.list` 第一列可以是任意区间信息，包括 `chromosome ID`，`chr:start-end` 或者 `.bed` 文件
- [2023-11-17] 为流程添加一个参数，自动将 WGS/WES 无必要的子 vcf 删除
- [2023-11-17] 设计出了 WES 模块，并将它和 WGS 融合，成为其子集
- [2023-11-09] 为 WGS 增加 --interval 参数作为 variant-calling-interval 
- [2023-11-08] 用 -I 代表输入参数，而不是原来的 -L
- [2023-11-06] 将 gatk 模块封装为 GATK class
- [2023-11-01] 添加基于 SOAPnuke 的原始数据过滤模块
- [2023-10-31] 完成 Sentieon 模块的添加，并且流程中不需要有额外的 CombineGVCF 这个步骤。
- [2023-10-31] 为 genotypeGVCFs 添加 --dbsnp 参数，用于注释已知位点
- [2023-10-30] 比对统计信息的计算和数据污染计算从 BQSR 前移到 markduplicate 这一步之后。并将 `.recal.table` 后缀改为 `.BQSR.recal.table`
- [2023-10-25] 删除配置文件中 vqsr_options 中的 `--max-gaussians`，改为 SNP(`vqsr_snp_max_gaussians`) 和 Indel (`vqsr_indel_max_gaussians`)单独设置.
- [2023-10-24] 添加 sentieon 模块，并基于 sentieon 构造 WGS/WES 分析流程.
- [2023-10-22] 修改配置文件的格式，将 `resources` 从 `gatk` 拎出来，单独成一块，方便 sentieon 和 GATK 共同调用.  
