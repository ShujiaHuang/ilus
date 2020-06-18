Ilus
====

**ilus** 是一个轻量的、可拓展的、简单易用的生成 **半自动化** 全基因组测序数据（Whole genome sequencing, WGS）分析流程的Python包.

简介
----

**ilus** 的作用是生成数据分析流程，而不是执行具体的流程。具体的执行将在生成流程脚本之后，由用户自己手动完成，该过程独立于 **ilus**，这也是为何称之为半自动化的原因. 目前 **Ilus** 中含有三个功能模块，分别是：

- 第一、``WGS`` 全基因组数据分析流程模块，该功能模块基于 `GATK 的最佳实践 <https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows>`_，使用 ``bwa-mem + GATK`` 的构建分析流程，其中包含了比对、排序、同一个样本多lane数据合并、重复序列标记（Markduplicates）、样本 gvcf 生成、多样本联合变异检测（Joint-calling）和 变异质控（Variant quality score recalibrator, VQSR） 这5个过程。该模块同样适用于 ``WES`` 数据的分析，只需要将配置文件的 ``variant_calling_interval`` 设置为 WES 测序的外显子捕获区间即可。
- 第二、``genotype-joint-calling`` 多样本联合变异检测模块。该模块是从 ** ilus WGS** 中分离出来的。这样可以满足对已有 gvcf 样本单独进行变异检测的需要。或者碰到分多批次完成的WGS分析时，如果需要合并做变异检测，也只需要整理一个 gvcf 列表文件，并用该功能就可以了，没有必要从测序数据 fastq 开始。
- 第三、``VQSR`` 变异质控模块（VQSR），同样是从 **ilus WGS** 中分离出来的，这样可以方便对单独的变异检测结果（VCF格式）进行变异质控。

    需要注意的是 ilus 不包含原始 fastq 数据的质控，ilus WGS 流程默认你所输入的测序数据都是 clean data， 即已经经过了严格的质控。

**ilus** 不具体执行任务的原因有两点：

- 第一，避免在 ilus 的核心功能之外做过多的优化。不同的计算集群（本地和云上），任务作业被调度的方式或大或小都会存在差异，如果将这些情况都一一考虑进去，那么 ilus 会不可避免地变得臃肿复杂，而且还不一定能够覆盖所有的类型，同时还容易丢失 **ilus** 的重点。因此，作为一个轻量级的工具，``ilus`` 在设计时就不将自动投递和运行任务的功能考虑在内，而是严格依据你的输入数据和配置文件信息，生成分析流程的各个 Bash 执行脚本。

    如果要实现作业的自动投递和监控，**ilus** 希望将来和开发者协作通过设计外部插件的形式来实现，但目前你需要手动将这些执行脚本，分步骤投递和运行。

- 第二，增加灵活性和操控性。**ilus** 所生成的流程运行脚本有一个特点，它们都完全独立于 ilus，不依赖 ilus 的任何功能，这个时候就算是将 ilus 完全删除掉都没关系。而且脚本中的每一行都是可以独立运行的，彼此之间互不影响。这样做的好处是，你可以按照计算集群的特点将这个脚本拆分成若干个子脚本（当然如果任务不多或者集群资源不充足，也可以不拆分），然后分别独立投递任务。但具体要拆成多少个取决于你要并行成多少个任务。比如你一共有 10 个样本，第一步的 ``bwa`` 比对脚本中就一共有10个比对命令，每行代表一个样本的 ``bwa``，由于每行命令都是彼此独立互不依赖的。因此你可以将该步骤步的执行脚本拆分为 10 （或者小于10）个子任务，然后分别手动投递这10个任务。至于如何将一个完整的执行脚本拆分为多个，你既可以自己写程序完成，也可以使用 **ilus** 提供的 `yhbatch_slurm_jobs程序 <https://github.com/ShujiaHuang/ilus/blob/master/scripts/yhbatch_slurm_jobs.py>`_ 来完成，如果该程序不能直接满足你的需求，你可以对其进行修改。虽然 **ilus** 没有自动化管理任务的投递和执行，但是通过这个方法却可以增加灵活性和操控性。

任务完成状态的监控是非常重要的，特别是当有成千上万的样本需要完成分析时，就显得更加重要，这个过程如果是人工手动来检查，那么过程一定是低效且易出错的。虽然 **ilus** 目前没有提供实时的任务监控工具，但是也已经实现了一个程序专门用于检查用户的作业是否都已成功完成，具体用法见下文例子。

如何安装
-------

Ilus 是基于 Python 编写的，稳定版本的代码已经发布至 PyPI。因此要使用 Ilus, 你可以通过 ``pip`` 这个Python包管理工具进行安装：

.. code:: bash

    pip install ilus

该命令除了主程序 ilus 之外，还会自动地将 ``ilus`` 所依赖的 Python 包也会被自动装上。安装完成之后，在命令行中执行 ``ilus`` ，如果能看到类似如下的内容，那么就说明你已经安装成功了。


.. code:: bash

    $ ilus
    usage: ilus [-h] {WGS,genotype-joint-calling,VQSR} ...
    ilus: error: too few arguments


快速开始
-------

通过执行 ``ilus --help`` 可以查看三个功能模块，分别是 ``WGS``, ``genotype-joint-calling`` 和 ``VQSR``。

.. code:: bash

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


下面，通过例子分别对这三个功能的使用进行说明。

全基因组数据分析
--------------

全基因组数据分析流程的运行脚本通过 ``ilus WGS`` 来生成，用法如下：

.. code:: bash

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
      


其中，``-C``, ``-L`` 和 ``-O`` 这三个参数是 **必须参数** ，其它的参数按照我们的实际需要做选择即可。``-O`` 参数比较简单，为输出目录，该目录如果不存在，**ilus** 将会自动创建。最重要的是 ``-C`` 和 ``-L`` 参数，前者是 **ilus** 的配置文件，没有这个文件，**ilus** 就无法生成正确的流程，因此十分重要；后者是输入文件的列表文件，该列表文件一共有 5 列，每一列都是必须的信息。

以下分别对这两个参数的格式展开说明：

首先是配置文件，我们需要在其中指定 ``WGS`` 流程各个步骤中所用的程序的路径以及所使用到 ``GATK bundle`` 文件和参考序列的路径。

需要注意的是 ``BWA MEM`` 的索引文件前缀需要与配置文件的 {resources}{reference} 相同，并存放在同一个目录中。如下：

.. code:: bash

    /path/human_reference/GRCh38/
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.dict
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.amb
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.ann
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.bwt
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.pac
    `-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.sa


该配置文件使用 Yaml 语法进行编写，在这里我提供一份 `配置文件的例子 <https://github.com/ShujiaHuang/ilus/blob/master/tests/ilus_sys.yaml>`_，参考如下：

.. code:: yaml

    aligner:
      bwa: /path_to/bwa
      bwamem_options: [-Y -M -t 8]

    samtools:
        samtools: /path_to/samtools
        sort_options: ["-@ 8"]
        merge_options: ["-@ 8 -f"]
        stats_options: ["-@ 8"]

    bcftools:
        bcftools: /path_to/bcftools
        options: []

    bedtools:
        bedtools: /path_to/bedtools
        concat_options: []
        genomecov_options: ["-bga -split"]

    # https://github.com/Griffan/VerifyBamID
    verifyBamID2:
        verifyBamID2: /path_to/verifyBamID2
        options: [
            "--SVDPrefix /path_to/verifyBamID2_resource/1000g.phase3.10k.b38.vcf.gz.dat"
        ]


    bgzip: /path_to/bgzip
    tabix: /path_to/tabix

    gatk:
      gatk: /path_to/gatk
      markdup_java_options: ["-Xmx10G", "-Djava.io.tmpdir=/your_path/cache"]
      bqsr_java_options: ["-Xmx8G", "-Djava.io.tmpdir=/your_path/cache"]
      hc_gvcf_java_options: ["-Xmx4G"]
      genotype_java_options: ["-Xmx8G"]
      vqsr_java_options: ["-Xmx10G"]

      CollectAlignmentSummaryMetrics_jave_options: ["-Xmx10G"]

      # Adapter sequencing of BGISEQ-500. If you use illumina (or others) sequencing system you should
      # change the value of this parameter.
      CollectAlignmentSummaryMetrics_options: [
        "--ADAPTER_SEQUENCE AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        "--ADAPTER_SEQUENCE AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"
      ]

      genomicsDBImport_options: ["--reader-threads 12"]
      use_genomicsDBImport: false  # Do not use genomicsDBImport to combine GVCFs by default

      vqsr_options: [
        "-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum",
        "-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0",
        "--max-gaussians 6"
      ]

      # interval value could be a file which contain all interval regions in it or could be a list here
      interval: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                 "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
      
      # Specific variant calling intervals. The value could be a file in bed format (recommend) or a interval list,
      # and the value could be as the same as ``interval`` parameter above.
      # The first three columns in interval regions file must be ``Sequencing ID``, ``region start`` and ``region end``,e.g.:
      #         chr1    10001   207666
      #         chr1    257667  297968

      variant_calling_interval: ["./wgs_calling_regions.GRCh38.interval.bed"]
      # variant_calling_interval: [
      #  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
      #  "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
      #  "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
      #  "chrX", "chrY", "chrM"
      #]

      bundle:
        hapmap: /path_to/gatk/bundle/hg38/hapmap_3.3.hg38.vcf.gz
        omni: /path_to/gatk/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
        1000G: /path_to/gatk/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        mills: /path_to/gatk/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        1000G_known_indel: /path_to/gatk/bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz
        dbsnp: /path_to/gatk/bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz


    # Define resources to be used for individual programs on multicore machines.
    # These can be defined specifically for memory and processor availability.
    resources:
      reference: /path_to/human_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa


在配置文件中， ``bwa``、``samtools``、``bcftools``、``bedtools``、``gatk``、``bgzip`` 和 ``tabix`` 都是必须的生信软件，需要自行安装，并将路径填入到对应的参数中，`verifyBamID2 <https://github.com/Griffan/VerifyBamID>`_ 仅用于计算样本是否存在污染，并不是必填的参数，如果配置文件中没有这个参数，则代表流程不会对样本的污染情况进行计算。另外，所必须的数据则是：``gatk bundle`` 和基因组参考序列。


``-L`` 参数是输入文件，文件中包含了WGS/WES分析流程所必须的所有测序数据信息，各列的信息如下：

- [1] SAMPLE，样本名
- [2] RGID，Read Group，使用bwa mem时通过 -R 参数指定的 read group
- [3] FASTQ1，Fastq1 文件的路径
- [4] FASTQ2，Fastq2 文件路径，如果是Single End测序，没有fastq2，此时该列用空格代替
- [5] LANE，fastq 的 lane 编号

如果某个样本的测序量比较大，导致一个样本有多个 lane 数据，或者同一个 lane 的数据被拆分成了多个，这个时候不需要人工对这些 fastq 数据进行合并，只需要依照如上信息编写好即可。同一个样本的数据在流程中会在各个子数据跑完比对并完成排序之后自动对进行合并。下面给出这个输入文件的例子：

.. code:: bash

    #SAMPLE RGID    FASTQ1  FASTQ2  LANE
    HG002   "@RG\tID:CL100076190_L01\tPL:COMPLETE\tPU:CL100076190_L01_HG002\tLB:CL100076190_L01\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_2.clean.fq.gz  CL100076190_L01
    HG002   "@RG\tID:CL100076190_L02\tPL:COMPLETE\tPU:CL100076190_L02_HG002\tLB:CL100076190_L02\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_2.clean.fq.gz  CL100076190_L02
    HG003   "@RG\tID:CL100076246_L01\tPL:COMPLETE\tPU:CL100076246_L01_HG003\tLB:CL100076246_L01\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_2.clean.fq.gz   CL100076246_L01
    HG003   "@RG\tID:CL100076246_L02\tPL:COMPLETE\tPU:CL100076246_L02_HG003\tLB:CL100076246_L02\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_2.clean.fq.gz   CL100076246_L02
    HG004   "@RG\tID:CL100076266_L01\tPL:COMPLETE\tPU:CL100076266_L01_HG004\tLB:CL100076266_L01\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_2.clean.fq.gz   CL100076266_L01
    HG004   "@RG\tID:CL100076266_L02\tPL:COMPLETE\tPU:CL100076266_L02_HG004\tLB:CL100076266_L02\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_2.clean.fq.gz   CL100076266_L02
    HG005   "@RG\tID:CL100076244_L01\tPL:COMPLETE\tPU:CL100076244_L01_HG005\tLB:CL100076244_L01\tSM:HG005"  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_1.clean.fq.gz  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_2.clean.fq.gz  CL100076244_L01

以下提供使用 **ilus** 生成 WGS 分析流程的例子。


**例子1：从头开始生成 WGS 分析流程**

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -O output/

这个命令的意思是，项目 ``my_wgs``（-n）依据 ``ilus_sys.yaml`` 和 ``input.list`` 在输出目录 ``output`` 中生成 WGS 分析流程，并将最后的比对数据从 BAM 转为 CRAM (-c)。

在输出目录 ``output`` 一共有 4 个文件夹（如下），分别用于存放分析流程产生的各类数据。其中：

.. code:: bash
    
    00.shell/
    01.alignment/
    02.gvcf/
    03.genotype/

- ``00.shell`` 该目录是分析流程的汇集目录，其中，生成了分步骤生成了流程各个步骤的执行脚本，同时还包含一个日志文件目录： 

.. code:: bash

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


投递任务运行流程时，我们按顺序从 step1 执行到 step6 即可。``loginfo`` 目录记录了各个步骤各个样本的运行状态，我们可以检查各个步骤的 ``.o.log.list`` 日志文件，获得该样本是否成功结束的标记。如果成功了，那么在该日志文件的末尾会有一个 ``[xxxx] xxxx done`` 的标记。可以通过使用 **ilus** 提供的程序 `check_jobs_status <https://github.com/ShujiaHuang/ilus/blob/master/scripts/check_jobs_status.py>`_ 检查各个步骤是否已经全部顺利完成，如果有错那么该脚本会将未完成的任务输出，方便我们进行检测和重新执行。用法为：

.. code:: bash

    $ python check_jobs_status.py loginfo/01.alignment.o.log.list > bwa.unfinish.list

如果任务都是成功结束的，那么该 list 文件为空，并输出 ``** All Jobs done **``。

- ``01.alignment`` 用于存放各个样本的比对结果
- ``02.gvcf`` 用于存放各个样本的 ``gvcf`` 结果
- ``03.genotype`` 用于存放最后变异检测的结果

**例子2：只生成 WGS 流程中的某个/某些步骤**

有时，我们并打算（或者没有必要）从头到尾完整地将 WGS 流程执行下去，比如我们只想执行从 ``fastq`` 比对到生成 ``gvcf`` 这个步骤，暂时不想执行 ``genotype`` 和 ``VQSR``，那么这个时候我们可以通过 ``-P`` 参数指定特定的步骤来实现：

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P align,markdup,BQSR,gvcf -O ./output


这样就只生成从 ``bwa`` 到 ``gvcf`` 的执行脚本。

除此之外，如果某个 WGS 步骤跑错了，调整之后，需要重新更新时，你也可以用 ``-P`` 指定重跑特定的步骤。比如我想重生成 BQSR 这个步骤的运行脚本，那么就可以这样做：

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P BQSR -O ./output

需要注意的是，**ilus** 为了节省项目的空间，只会为每一个样本保留 BQSR 之后的 BAM/CRAM 文件，因此，如果你想重新跑 BQSR 需要确定在 BQSR 前一步（即，markdup）的 BAM 文件是否已经被删除了，如果原先 **ilus** 在BQSR这一步没有正常结束的话，那么该 markdup 的 BAM 文件应该还会被保留着的，**ilus** 执行任务时具有“原子属性”，也就是说只有当所有步骤都成功结束时才会删除在之后的分析中完全不需要的文件。

    -P 参数能够指定的步骤必须属于「align,markdup,BQSR,gvcf,genotype,VQSR」中的一个或多个.


genotype-joint-calling
----------------------

如果我们已经有了各个样本的 gvcf 数据，现在要用这些 gvcf 完成多样本的联合变异检测（Joint-calling），那么就可以使用 ``genotype-joint-calling`` 来实现。具体用法如下：

.. code:: bash

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


``-L`` 是 **ilus genotype-joint-calling** 的输入参数，它接受的是一个 ``gvcf list`` 文件，这个文件由两列构成，第一列是 gvcf 文件所对应的区间或者染色体编号，第二列是 gvcf 文件的路径，目前 **ilus** 要求各个样本的 gvcf 都按照主要染色体（1-22、X、Y、M）分开，举个例子：

.. code:: bash

    $ ilus genotype-joint-calling -n my_project -C ilus_sys.yaml -L gvcf.list -O genotype --as_pipe_shell_order

其中 ``gvcf.list`` 的格式如下：

.. code:: bash

    chr1    /path/sample1.chr1.g.vcf.gz
    chr1    /paht/sample2.chr1.g.vcf.gz
    chr2    /path/sample1.chr2.g.vcf.gz
    chr2    /path/sample2.chr2.g.vcf.gz
    ...
    chrM    /path/sample1.chrM.g.vcf.gz
    chrM    /path/sample2.chrM.g.vcf.gz

以上 ``gvcf.list`` 中只有两个样本。

参数 ``--as_pipe_shell_order`` 可加也可不加（默认是不加），它唯一的作用就是按照 **ilus WGS** 流程的方式输出执行脚本的名字，维持和 ``WGS`` 流程一样的次序和相同的目录结构。


VQSR
----

该功能仅用于生成基于 ``GATK VQSR`` 的执行脚本。我们如果已经有了最终的变异检测（VCF格式）结果，现在只想借助 ``GATK VQSR`` 完成这个变异数据的质控，那么就可以使用这个模块了，用法与 ``genotype-joint-calling`` 大同小异，如下：

.. code:: bash

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

跟 ``genotype-joint-calling`` 相比不同的是，**ilus VQSR** 中的输入文件是 VCF 文件列表，并且每行只有一列，为 vcf 文件的路径，举个例子，如下：

.. code:: bash

    /path/chr1.vcf.gz
    /path/chr2.vcf.gz
    ...
    /path/chrM.vcf.gz

**ilus VQSR** 的其它参数与 ``genotype-joint-calling`` 相同，以下提供一个完整的例子：

.. code:: bash

    $ ilus VQSR -C ilus_sys.yaml -L vcf.list -O genotype --as_pipe_shell_order



