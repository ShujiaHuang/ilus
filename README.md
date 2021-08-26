Ilus
====

**ilus** 是一个轻量的、可拓展的、易用的 **半自动化** 全基因组（Whole
genome sequencing, WGS）和 全外显子（Whole exom
sequencing，WES）分析流程生成器.

[English](./README_EN.md) | 简体中文

简介
----

**ilus** 是一个流程生成器，用途是生成 WGS/WES 分析流程，但 **ilus**
不执行流程的具体步骤。用户需要自己手动操作，而且执行过程不再 依赖
**ilus**，这也是为何称之为半自动化的原因。

目前 **ilus** 含有三个功能模块，分别是：

-   第一、`WGS` 全基因组数据分析流程模块。该模块基于 [GATK
    最佳实践](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)，采用
    `bwa-mem + GATK`
    这个目前最主流的方式构建分析流程。模块中集成了比对、排序、同一个样本多
    lane 数据合并、标记重复序列（Markduplicates）、HaplotypeCaller gvcf
    生成、多样本联合变异检测（Joint-calling）和变异质控（Variant quality
    score recalibrator, VQSR） 这 5 个完整的过程。这个模块同样适用于
    `WES` 数据的分析，只需要将配置文件 `variant_calling_interval` 设置为
    WES 的外显子捕获区间即可（下文详述）。

-   第二、`genotype-joint-calling` 联合变异检测模块。这个模块是从 **ilus
    WGS** 中分离出来的，目的是为了满足对已有 **gvcf**
    的样本直接进行变异检测的需要。或者碰到需要分多批次完成 WGS/WES
    数据分析时，可以分批次生成 **gvcf**，最后再整理一个总的 **gvcf**
    文件列表，然后用该功能完成后续步骤就可以了，这可以增加跑流程的灵活度。

-   第三、`VQSR` 变异质控模块，同样是从 **ilus WGS**
    中分离出来的，目的是方便我们对变异结果进行 `VQSR` 变异质控。

    > 需要注意的是 ilus 不包含对原始 fastq 数据的质控，ilus
    > 流程默认你输入的测序数据都是清洗好的 clean data。

**ilus** 不直接运行任务有如下两个方面的考虑：

-   第一，避免在 **ilus**
    的核心功能之外做关于任务调度方面的优化。不同的计算集群（本地和云上），作业被调度的方式是多种多样的，如果将这些情况都一一考虑进去，**ilus**
    会变得臃肿复杂，并且还不一定能够做得好，这反而会导致一部分人无法有效使用
    **ilus**，甚至也容易在管理任务投递的过程中丢失 **ilus**
    的重点。作为一个轻量级的工具，我在设计 **ilus**
    时从一开始就没将自动投递和运行任务的功能考虑在内。我希望它作为一个框架程序，严格依据你的输入数据和配置文件的信息，生成符合你分析需求的流程脚本。

    > 如果要实现作业的自动投递和监控，**ilus**
    > 希望将来可以通过外部插件的形式来实现，但目前需要手动去投递这些脚本。

-   第二，增加灵活性和可操控性。**ilus**
    所生成的流程都有一个特点，那就是它们都完全独立于
    **ilus**，不会再依赖 **ilus** 的任何功能，这个时候即使你将 **ilus**
    彻底卸载删除掉都没关系。另外，**ilus
    所生成的执行脚本(shell)中的每一行都是可以独立运行的，彼此之间互不影响**。

这样做的好处是，你可以按照计算机集群的特点将脚本拆分成若干个子脚本（如果你的样本不多或者集群资源不充足，也可以不拆分），然后再分别投递
任务，这样可以大大加快任务的完成速度，**而且，这是目前并行投递ilus任务的唯一方式**。

至于要拆成多少个子脚本取决于你自己的具体情况。举个例子，比如你一共有 10 个样本，第一步的 `xxx.step1.bwa.sh` 比对脚本中一共有 10 个比对命令，每一行都是一个样本的 `bwa`。

由于这 10 个命令彼此独立互不依赖，因此你可以将该步骤步拆分为 10（或者小于10）个子脚本，然后再手动投递这些任务。至于如何将一个完整的执行脚本拆分为多个，你既可以自己写程序完成，也可以使用 **ilus** 提供的
[yhbatch_slurm_jobs](https://github.com/ShujiaHuang/ilus/blob/master/scripts/yhbatch_slurm_jobs.py) 来完成，但要注意，我这里提供的程序是基于 slurm 集群的，不一定符合你的需要，但它的作用和意义我刚刚也做了说明，如果你觉得这个程序不能直接满足你的需求，你可以对其进行修改。虽然 **ilus** 没有自动化管理任务的投递和执行，但是通过这个方法却也能够增加流程控制的灵活性和操控性。

实现对任务完成状态的监控是非常重要的，特别是有成千上万的样本需要分析时，就显得更加重要。这个过程如果通过手动来检查，那么一定是低效且易出错的。**ilus**
充分考虑到了这一点，因此在生成流程的时候会为每个任务都添加一个可识别的完成标记，我们只需要查看每个任务是否有该标记就行了（具体情况参考下文**WGS** 的例子）。

虽然如此，一旦我们的任务数量庞大，如果每次都要手动打开文件检查对应的任务是否已经顺利结束的话，那么未免太过于麻烦。因此，我在 **ilus** 中实现了一个可以专门用于检查任务作业完成状态的程序，具体用法参考下文 **WGS** 的例子。

如何安装
--------

**ilus** 是基于 Python 编写的，同时支持 Python2.7+ 和 Python3.7+，稳定版本的代码已经发布至 PyPI。因此要使用 **ilus**, 直接通过 `pip` 这个 Python 包管理工具进行安装：

```bash
$ pip install ilus
```

该命令除了主程序 `ilus` 之外，还会自动将 **ilus** 所依赖的其它 Python 包自动装上。安装完成之后，在命令行中执行 `ilus`，如果能看到类似如下的内容，那么就说明安装成功了。

```bash
$ ilus
usage: ilus [-h] {WGS,genotype-joint-calling,VQSR} ...
ilus: error: too few arguments
```

快速开始
--------

通过执行 `ilus --help` 可以看到 `WGS`, `genotype-joint-calling` 和
`VQSR` 这三个功能模块。

```bash
$ ilus --help
usage: ilus [-h] {WGS,genotype-joint-calling,VQSR} ...

ilus: A WGS analysis pipeline.

optional arguments:
    -h, --help            show this help message and exit

ilus commands:
{WGS,genotype-joint-calling,VQSR}
    WGS                 Creating pipeline for WGS(from fastq to genotype VCF)
    genotype-joint-calling Genotype from GVCFs.
    VQSR                VQSR
```

下面，通过例子逐一介绍如何使用这三个模块。

### 全基因组数据分析

全基因组数据分析流程（WGS）的运行脚本通过 `ilus WGS` 来生成，用法如下：

```bash
$ ilus WGS --help
usage: ilus WGS [-h] -C SYSCONF -L FASTQLIST [-P WGS_PROCESSES]
            [-n PROJECT_NAME] [-f] [-c] -O OUTDIR

optional arguments:
  -h, --help            show this help message and exit
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
  -L FASTQLIST, --fastqlist FASTQLIST
                        Alignment FASTQ Index File.
  -O OUTDIR, --outdir OUTDIR
                        A directory for output results.

  -n PROJECT_NAME, --name PROJECT_NAME
                        Name of the project. Default value: test
  -P WGS_PROCESSES, --Process WGS_PROCESSES
                        Specific one or more processes (separated by comma) of
                        WGS pipeline. Defualt value:
                        align,markdup,BQSR,gvcf,genotype,VQSR. Possible
                        values: {align,markdup,BQSR,gvcf,genotype,VQSR}
  -f, --force_overwrite
                        Force overwrite existing shell scripts and folders.
  -c, --cram            Covert BAM to CRAM after BQSR and save alignment file storage.
```

其中，`-C`, `-L` 和 `-O` 这三个是 **必须参数**，其余的按照实际需要做选择。`-O` 参数比较简单，是输出目录，该目录如果不存在，**ilus**将自动创建。最重要的是 `-C` 和 `-L` 参数，前者是 **ilus** 的配置文件，如果没有这个文件 **ilus** 就无法正确生成分析流程，因此它十分重要；后者是输入文件，**这个文件的格式有固定要求**，一共 5 列，每一列都是流程所必须的信息。

下面，我分别对这两个文件的格式展开说明：

首先是配置文件，我们需要在文件中编写好分析流程所用的程序路径、`GATK bundle` 文件路径、参考序列的路径以及各个关键步骤所对应的参数。

需要注意的是 `bwa mem` 的比对索引文件前缀要与配置文件的 `{resources}{reference}` 的前缀相同，并放在同一个文件夹里。如下：

```bash
/path/human_reference/GRCh38/
|-- human_GRCh38.fa
|-- human_GRCh38.dict
|-- human_GRCh38.fa.amb
|-- human_GRCh38.fa.ann
|-- human_GRCh38.fa.bwt
|-- human_GRCh38.fa.fai
|-- human_GRCh38.fa.pac
`-- human_GRCh38.fa.sa
```

配置文件要使用 [Yaml 语法](https://zh.wikipedia.org/wiki/YAML) 进行编写，这里我提供一份 [配置文件的模板](https://github.com/ShujiaHuang/ilus/blob/master/tests/ilus_sys.yaml)，参考如下：

```yaml
# Configuration file specifying system details for running an analysis pipeline
aligner:
  bwa: /path/to/BioSoftware/local/bin/bwa
  bwamem_options: [-Y -M -t 8]

samtools:
    samtools: /path/to/BioSoftware/local/bin/samtools
    sort_options: ["-@ 8"]
    merge_options: ["-@ 8 -f"]
    stats_options: ["-@ 8"]

bcftools:
    bcftools: /path/to/BioSoftware/local/bin/bcftools
    concat_options: ["-a --rm-dups all"]

bedtools:
    bedtools: /path/to/BioSoftware/local/bin/bedtools
    genomecov_options: ["-bga -split"]

sambamba:
  sambamba: /path/to/BioSoftware/local/bin/sambamba
  sort_options: ["-t 8"]
  merge_options: ["-t 8"]
  markdup_options: []


verifyBamID2:
    # This is the VerifyBamID2: https://github.com/Griffan/VerifyBamID
    verifyBamID2: /path/to/BioSoftware/local/bin/verifyBamID2
    options: [
        # download from: https://github.com/Griffan/VerifyBamID/tree/master/resource 
        "--SVDPrefix /path/to/BioSoftware/verifyBamID2/1.0.6/resource/1000g.phase3.10k.b38.vcf.gz.dat"
    ]


bgzip: /path/to/BioSoftware/local/bin/bgzip
tabix: /path/to/BioSoftware/local/bin/tabix

gatk:
  gatk: /path/to/BioSoftware/gatk/4.1.4.1/gatk
  markdup_java_options: ["-Xmx10G", "-Djava.io.tmpdir=/your_path/cache"]
  bqsr_java_options: ["-Xmx8G", "-Djava.io.tmpdir=/your_path/cache"]
  hc_gvcf_java_options: ["-Xmx4G"]
  genotype_java_options: ["-Xmx8G"]
  vqsr_java_options: ["-Xmx10G"]

  CollectAlignmentSummaryMetrics_jave_options: ["-Xmx10G"]

  # Default adapter sequence is BGISEQ-500/MGISEQ/DNBSEQ in ilus. If you use illumina (or other) sequencing system 
  # you should change the value of this parameter. The most widely used adapter of Illumina is TruSeq adapters. If 
  # your data is from the TruSeq library, you can replace the parameter with the two following sequences: 
  # "--ADAPTER_SEQUENCE AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  # "--ADAPTER_SEQUENCE AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

  CollectAlignmentSummaryMetrics_options: [
    "--ADAPTER_SEQUENCE AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    "--ADAPTER_SEQUENCE AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"
  ]

  hc_gvcf_options: [""]
  genotypeGVCFs_options: [""]
  genomicsDBImport_options: ["--reader-threads 12"]
  use_genomicsDBImport: false  # Do not use genomicsDBImport to combine GVCFs by default

  vqsr_options: [
    "-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff",
    "-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0",
    "--max-gaussians 6"
  ]

  # Fro creating gvcf. The value could be a interval region file in bed format or 
  # could be chromosomes list here. I suggest you to use chromosome list here.
  interval: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
             "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
  

  # Specific variant calling intervals. 
  # The value could be a file in bed format (I show you a example bellow) or a interval of list.
  # Bed format of interval file only contain three columns: ``Sequencing ID``, ``region start`` and ``region end``,e.g.:
  #         chr1    10001   207666
  #         chr1    257667  297968

  # These invertals could be any regions alone the genome as you wish or just set the same as ``interval`` parameter above.
  # variant_calling_interval: ["./wgs_calling_regions.GRCh38.5M.interval.bed"]
  variant_calling_interval: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                             "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                             "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
  

  # GATK bundle
  bundle:
    hapmap: /path/to/BioDatahub/gatk/bundle/hg38/hapmap_3.3.hg38.vcf.gz
    omni: /path/to/BioDatahub/gatk/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
    1000G: /path/to/BioDatahub/gatk/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    mills: /path/to/BioDatahub/gatk/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    1000G_known_indel: /path/to/BioDatahub/gatk/bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz
    dbsnp: /path/to/BioDatahub/gatk/bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz


# Define resources to be used for individual programs on multicore machines.
# These can be defined specifically for memory and processor availability.
resources:
  reference: /path/to/BioDatahub/human_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

```

在配置文件中，`bwa`、`samtools`、`bcftools`、`bedtools`、`gatk`、`bgzip` 和 `tabix` 都是必须的生信软件，需要自行安装，再将路径填入到对应的参数里（如模板所示）。[verifyBamID2](https://github.com/Griffan/VerifyBamID) 仅用于计算样本是否存在污染，**它并不是一个必填的参数**，如果你的配置文件中没有这个参数，则代表流程不对样本的污染情况进行推算，如果有那么你要自行安装并下载与之配套的 `resource` 数据，模板里我也告诉你该去哪里下载相关的数据了。

要注意的是，配置文件中的 `variant_calling_interval` 参数。这是一个专门用来指定变异检测区间的参数，比如以上配置文件的例子，我给出了从 `chr1` 到 `chrM` 这 25 条染色体，意思就是告诉流程要对这 25 条染色体做变异检测。如果你在这个参数里只列出一条染色体，或者仅仅给出一个染色体区间，比如 `chr1:1-10000`，那么 **ilus** 也将只在你给定的这个区间里完成变异检测。

这是一个非常灵活有用的参数，`variant_calling_interval` 区间是可以任意指定的，除了可以按照我例子给出的赋值方式之外，还可以将区间 **文件的路径** 赋给这个参数。**我们知道 WGS 和 WES 有很多步骤是完全相同的，只在变异检测的区间上存在差别------WES数据没有必要也不能** 在全染色体上做变异检测，只在外显子捕获区域里进行就可以了。

这个时候你只需要将外显子捕获区域的文件------注意是文件，这个文件的内容是 `.bed` 文件格式，例子如下所示：

```bash
chr1    63697   63697
chr1    101158  101158
chr1    103241 103241
chr1    104108  104108
chr1    185336 185336
chr1    261495  261495
chr1    598862 598862
chr1    601606  601606
chr1    700596 700596
chr1    725086  725086
```

也可以参照 [GATK的说明](<https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists>)，这里你不需要手动拆分成一个个的区间，只需将文件的路径赋给这个参数就可以了，这时流程就成了 WES 分析流程。这也是为何称 `ilus` 是一个WGS和WES分析流程生成器的原因。

另外，**ilus** 必需的公用数据集是：`gatk bundle` 和基因组参考序列。

> 【注意】如果你项目的样本量少于 10 那么 GATK 将不计算 `InbreedingCoeff` 的值，此时配置文件中 `vqsr_options` 不需要设置 `-an InbreedingCoeff`，可以将其去掉。

接下来是由 `-L` 参数指定的输入文件，文件里包含了 `WGS/WES` 分析流程所必需的一切测序数据信息，**这个文件需要你自己来准备**，文件各列的格式信息如下：

- [1] SAMPLE，样本名

- [2] RGID，Read Group，使用 `bwa mem` 时通过 -R 参数指定的 `read group`

- [3] FASTQ1，Fastq1 文件的路径

- [4] FASTQ2，Fastq2 文件路径，如果是`Single End`测序，没有`fastq2`，此时该列用空格代替

- [5] LANE，fastq 的 `lane` 编号

    > 这五个信息中 `RGID` 最容易出错，`RGID`一定要设置正确（正确的编写方式参考以下例子），否则分析流程会出错。

另外，假如某个样本的测序量比较大，导致一个样本有多个 `lane` 的测序数据，或者同一个 `lane` 的数据被拆分成了多个子文件，这个时候不需要人工对这些 `fastq` 数据进行合并，只需要依照测序信息编写好这个输入文件即可。

那些属于同一个样本的数据，即使输入的 `fastq` 已被拆分成了成千上万份，流程中也会在各个子数据跑完比对和排序之后自动进行合并。
下面我给出一个输入文件的例子，其中就有样本的数据分拆输出的情况：

```bash
#SAMPLE RGID    FASTQ1  FASTQ2  LANE
HG002   "@RG\tID:CL100076190_L01\tPL:COMPLETE\tPU:CL100076190_L01_HG002\tLB:CL100076190_L01\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_2.clean.fq.gz  CL100076190_L01
HG002   "@RG\tID:CL100076190_L02\tPL:COMPLETE\tPU:CL100076190_L02_HG002\tLB:CL100076190_L02\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_2.clean.fq.gz  CL100076190_L02
HG003   "@RG\tID:CL100076246_L01\tPL:COMPLETE\tPU:CL100076246_L01_HG003\tLB:CL100076246_L01\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_2.clean.fq.gz   CL100076246_L01
HG003   "@RG\tID:CL100076246_L02\tPL:COMPLETE\tPU:CL100076246_L02_HG003\tLB:CL100076246_L02\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_2.clean.fq.gz   CL100076246_L02
HG004   "@RG\tID:CL100076266_L01\tPL:COMPLETE\tPU:CL100076266_L01_HG004\tLB:CL100076266_L01\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_2.clean.fq.gz   CL100076266_L01
HG004   "@RG\tID:CL100076266_L02\tPL:COMPLETE\tPU:CL100076266_L02_HG004\tLB:CL100076266_L02\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_2.clean.fq.gz   CL100076266_L02
HG005   "@RG\tID:CL100076244_L01\tPL:COMPLETE\tPU:CL100076244_L01_HG005\tLB:CL100076244_L01\tSM:HG005"  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_1.clean.fq.gz  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_2.clean.fq.gz  CL100076244_L01
```

接下来举例说明 **ilus WGS** 的使用和流程结构情况。

**例子1：从头开始生成 WGS 分析流程**

```bash
$ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -O output/
```

这个命令的意思是，项目 (-n) `my_wgs` 依据配置文件 (-C) `ilus_sys.yaml` 和输入数据(-L)`input.list` 在输出目录 `output` 中生成一个 WGS 分析流程。同时流程在完成分析之后将 `BAM` 自动转为 (-c)`CRAM` 格式。`CRAM` 比 `BAM` 更加节省空间，如果不设置 `-c` 参数，则保留原来的 `BAM` 文件。

以上命令顺利执行之后，在输出目录 `output` 中一共有 4 个文件夹（如下）：

```bash
00.shell/
01.alignment/
02.gvcf/
03.genotype/
```

它们分别用于存放流程产生的各类不同数据，其中：

-   `00.shell` 流程 `shell` 脚本的汇集目录；
-   `01.alignment` 以样本为单位存放比对结果；
-   `02.gvcf` 存放各个样本的 `gvcf` 结果；
-   `03.genotype` 存放最后变异检测的结果。

`00.shell` 目录里有分析流程的各个执行脚本和日志目录：

```bash
/00.shell
├── loginfo
│   ├── 01.alignment
│   ├── 01.alignment.e.log.list
│   ├── 01.alignment.o.log.list
│   ├── 02.markdup
│   ├── 02.markdup.e.log.list
│   ├── 02.markdup.o.log.list
│   ├── 03.BQSR
│   ├── 03.BQSR.e.log.list
│   ├── 03.BQSR.o.log.list
│   ├── 04.gvcf
│   ├── 04.gvcf.e.log.list
│   ├── 04.gvcf.o.log.list
│   ├── 05.genotype
│   ├── 05.genotype.e.log.list
│   ├── 05.genotype.o.log.list
│   ├── 06.VQSR
│   ├── 06.VQSR.e.log.list
│   └── 06.VQSR.o.log.list
├── my_wgs.step1.bwa.sh
├── my_wgs.step2.markdup.sh
├── my_wgs.step3.bqsr.sh
├── my_wgs.step4.gvcf.sh
├── my_wgs.step5.genotype.sh
└── my_wgs.step6.VQSR.sh
```

投递任务运行流程时，按顺序从 `step1` 依次执行到 `step6` 即可。`loginfo/` 文件夹下记录了各个样本所有步骤的运行状态，你可以通过检查各个任务的 `.o.log.list` 日志文件，获得每个样本是否都成功结束的标记。

如果成功了，可以在日志文件的末尾看到一个类似于 `[xxxx] xxxx done` 的标记。通过我在 **ilus** 中提供的程序 [check_jobs_status](https://github.com/ShujiaHuang/ilus/blob/master/scripts/check_jobs_status.py) 你可以很方便地知道哪些样本（步骤）已经顺利完成，哪些还没有。这个脚本会帮你将那些未完成的任务全部输出，方便检查问题和重新执行这部分未完成的任务。`check_jobs_status`
用法如下：

```bash
$ python check_jobs_status.py loginfo/01.alignment.o.log.list > bwa.unfinish.list
```

如果这个 list 文件为空，并输出了 `** All Jobs done **`，那么代表所有任务都成功结束了。

### 如何并行投递任务

**ilus** 所生成的流程脚本有一个特点，那就是它们都完全独立于 **ilus**，不会再依赖 **ilus** 的任何功能，这个时候你就算将 **ilus** 卸载删除掉都没关系。**而且执行脚本(shell)中的每一行都是可以独立运行的，彼此之间互不影响**。

这样做的好处是，你可以按照所用集群的特点将脚本拆分成若干个子脚本（如果你的样本不多或者集群资源不充足，也可以不拆分），然后再分别独立投递任务，这样可以大大提升任务完成速度，**这个也是目前并行投递ilus任务的唯一方式**。

至于具体要拆成多少个子脚本取决于你自己的需要，比如你一共有 10 个样本，第一步的 `xxx.step1.bwa.sh` 比对脚本中就一共有10个比对命令，每一行都是一个样本的 `bwa`。由于这10个命令都彼此独立互不依赖的，因此你可以将该步骤步拆分为 10（或者小于10）个子脚本，然后再手动投递这 10 个任务。至于如何将一个完整的执行脚本拆分为多个，你既可以自己写程序完成，也可以使用 **ilus** 提供的 [yhbatch_slurm_jobs程序](https://github.com/ShujiaHuang/ilus/blob/master/scripts/yhbatch_slurm_jobs.py) 来完成，但要注意，我这里提供的这个程序是基于 slurm 集群的，不一定符合你的需要，但它的作用和意义我刚刚也说了，如果你觉得这个程序不能满足你的需求，你可以进行修改。

**例子2：只生成 `WGS` 流程中的某个/某些步骤**

有时，我们并不打算（或者没有必要）从头到尾完整地将 `WGS` 流程执行下去，比如我们只想执行从 `fastq` 比对到生成 `gvcf` 这个步骤，暂时不想执行 `genotype` 和 `VQSR`，该怎么办呢？ilus 的 `-P` 参数就可以实现这个目的：

```bash
$ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P align,markdup,BQSR,gvcf -O ./output
```

这样就只生成从 `bwa` 到 `gvcf` 的执行脚本，这对于需要分批次完成分析的项目来说是很有用的。而且 `ilus` 所输出的结果是以样本为单位作区分的，因此在相同的输出目录下，只要样本编号是不同的，那么不同批次的数据就不会存在相互覆盖的问题。

除此之外，`-P` 参数还有一个用途，那就是假如某个 `WGS` 步骤跑错了，需要调整，之后再重新更新对应的步骤，那你就可以用 `-P` 重跑特定的步骤。比如我需要重新生成 `BQSR` 这个步骤的运行脚本，那么就可以这样做：

```bash
$ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P BQSR -O ./output
```

不过，要注意的是，**ilus** 为了节省项目对存储空间的消耗，只会为每一个样本保留 `BQSR` 之后的总 `BAM/CRAM` 文件。因此，如果你想重新跑 `BQSR` 那就需要先确保 `BQSR`
前一步（即，`markdup`）的 `BAM` 文件没有被被删除。

如果你一直使用的是 **ilus** 那么是不用担心这个问题的，因为 **ilus** 执行任务时具有 "原子属性"，也就是说只有当步骤中所有过程都成功结束了才会将那些完全不需要的文件删除掉。所以，如果 `BQSR` 这一步没有正常结束，那么前一步 `markdup` 的 `BAM` 文件是会被保留着的。

> -P
> 参数用来指定的分析模块必须属于「align,markdup,BQSR,gvcf,genotype,VQSR」中的一个或多个，并用英文逗号隔开。

### genotype-joint-calling

如果我们已经有了各个样本的 `gvcf` 数据，现在要用这些 `gvcf` 完成多样本的联合变异检测 `Joint-calling`，那么就可以使用
`genotype-joint-calling` 来实现。具体用法如下：

```bash
$ ilus genotype-joint-calling --help
usage: ilus genotype-joint-calling [-h] -C SYSCONF -L GVCFLIST
                                   [-n PROJECT_NAME] [--as_pipe_shell_order]
                                   [-f] -O OUTDIR

optional arguments:
  -h, --help            show this help message and exit
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
  -L GVCFLIST, --gvcflist GVCFLIST
                        GVCFs file list. One gvcf_file per-row and the format
                        should looks like: [interval gvcf_file_path]. Column
                        [1] is a symbol which could represent the genome
                        region of the gvcf_file and column [2] should be the
                        path.
  -O OUTDIR, --outdir OUTDIR
                        A directory for output results.
  -n PROJECT_NAME, --name PROJECT_NAME
                        Name of the project. [test]
  --as_pipe_shell_order
                        Keep the shell name as the order of `WGS`.
  -f, --force           Force overwrite existing shell scripts and folders.
```

`-L` 是 **ilus genotype-joint-calling** 的输入参数，它接受的是一个 `gvcf list` 文件，这个文件由两列构成，第一列是每个 `gvcf` 文件所对应的区间或者染色体编号，第二列是 `gvcf` 文件的路径。目前 **ilus** 要求各个样本的 `gvcf` 都按照主要染色体（1-22、X、Y、M）分开，举个例子：

```bash
$ ilus genotype-joint-calling -n my_project -C ilus_sys.yaml -L gvcf.list -O genotype --as_pipe_shell_order
```

其中 `gvcf.list` 的格式如下：

```bash
chr1    /path/sample1.chr1.g.vcf.gz
chr1    /paht/sample2.chr1.g.vcf.gz
chr2    /path/sample1.chr2.g.vcf.gz
chr2    /path/sample2.chr2.g.vcf.gz
...
chrM    /path/sample1.chrM.g.vcf.gz
chrM    /path/sample2.chrM.g.vcf.gz
```

这个例子里 `gvcf.list` 只有两个样本。参数 `--as_pipe_shell_order` 可加也可不加（默认是不加），它唯一的作用就是按照 **ilus WGS**流程的方式输出执行脚本的名字，维持和 `WGS` 流程一样的次序和相同的输出目录结构。

### VQSR

该功能仅用于生成基于 `GATK VQSR` 的执行脚本。我们如果已经有了最终的变异检测（VCF格式）结果，现在只想借助 `GATK VQSR` 对这个变异数据做质控，那么就可以使用这个模块了，用法与 `genotype-joint-calling` 大同小异，如下：

```bash
$ ilus VQSR --help
usage: ilus VQSR [-h] -C SYSCONF -L VCFLIST [-n PROJECT_NAME]
                 [--as_pipe_shell_order] [-f] -O OUTDIR

optional arguments:
  -h, --help            show this help message and exit
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
  -L VCFLIST, --vcflist VCFLIST
                        VCFs file list. One vcf_file per-row and the format
                        should looks like: [interval vcf_file_path]. Column
                        [1] is a symbol which could represent the genome
                        region of the vcf_file and column [2] should be the
                        path.
  -O OUTDIR, --outdir OUTDIR
                        A directory for output results.
  -n PROJECT_NAME, --name PROJECT_NAME
                        Name of the project. [test]
  --as_pipe_shell_order
                        Keep the shell name as the order of `WGS`.
  -f, --force           Force overwrite existing shell scripts and folders.
```

跟 `genotype-joint-calling` 相比不同的是，**ilus VQSR** 的输入文件是 `VCF` 文件列表，并且 **每行就是一个VCF文件的路径**，举个例子，如下：

```bash
/path/chr1.vcf.gz
/path/chr2.vcf.gz
...
/path/chrM.vcf.gz
```

其它参数与 `genotype-joint-calling` 相同。还有文件列表中的vcf不需要事先进行手动合并，`ilus VQSR` 会帮你合并。以下提供一个完整的例子：

```bash
$ ilus VQSR -C ilus_sys.yaml -L vcf.list -O genotype --as_pipe_shell_order
```

以上，文档结束。








