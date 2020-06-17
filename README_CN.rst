Ilus
====

**Ilus** 是一个轻量级的、半自动的、可拓展的、简单易用的全基因组测序数据（Whole genome sequencing, WGS）分析流程包.

简介
----

**Ilus** 目前包含了三个独立功能模块，分别是：

* 第一、全基因组数据分析流程，该流程基于 `GATK 的最佳实践 <https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows>`_，即使用 bwa-mem + GATK 的方式进行构建，包含了比对、排序、同样本多lane数据合并、标记重复序列（Markduplicates）、单样本 gvcf 生成、多样本联合变异检测（Joint-calling）和 变异质控（Variant quality score recalibrator, VQSR） 这若干个过程。
* 第二、独立的多样本联合变异检测模块，之所以将该部分从原来的 WGS 流程中分离出来作为一个独立的功能，这是因为，在有些情况下，我们已经有了各个样本的 gvcf ，或者我们的数据是分多批次完成的，那么这个时候，只需要整理一个 gvcf 列表文件，并用该功能就可以接下去完成多样本的联合变异检测了，而没有必要从测序数据 fastq 开始从头分析。
* 第三、独立的变异质控模块（VQSR），该模块同样是为了方便对单独的变异检测结果（VCF格式）进行变异质控。


    需要注意的是 ilus 不包含原始 fastq 数据的质控，ilus WGS 流程默认你所输入的测序数据都是 clean data， 即已经经过了严格的质控。

另外，考虑到不同计算集群（本地和云上）的作业调度情况存在一定的区别，``ilus`` 并不会直接在你的计算机或者集群里自动投递和运行任务，而是严格依据你的输入数据（包括数据的顺序）和配置文件的信息，生成分析流程的各个 Bash 执行脚本。你需要手动将这些执行脚本，分步骤投递和运行，并且在投递下一个任务之前，要先确保上一步已经正确完成。这样做的目的是希望能够尽可能地满足不同用户使用不同计算环境时的需要，而且，这些步骤的分析脚本都会统一生成并汇集在一个固定的目录下（目录名通常是 ``00.shell``）。

除了可以按照顺序投递任务之外，假如你有充足的计算资源，那么你可以按照需要将一个汇总的运行脚本拆分成若干个子执行脚本，具体拆成多少个取决于你要并行成多少个任务。比如你一共有 10 个样本，第一步的 bwa 比对脚本中就一共有10行，每行代表一个样本的bwa，这里每一行的运行命令都是可以并行运行的，每个任务之间是彼此独立的。因此你可以将该步骤步的执行脚本拆分为 10 个（或者其它小于10的个数）子任务，然后分别手动投递这10个任务就行了。至于如何将一个完整的执行脚本拆分为多个，你既可以自己写程序完成，也可以使用 **ilus** 提供的程序 ``yhbatch_slurm_jobs.py`` 来完成。这样做虽然牺牲了自动化的任务投递，但却可以增加适应性和操控性。



如何安装
-------

Ilus 基于 Python 编写，并已经发布至 PyPI，因此要使用 Ilus, 你只需要执行以下命令即可：

.. code:: bash

    pip install ilus

该命令除了主程序 ilus 之外，还将自动地将 ilus 所依赖的其它 Python 包也一并装上。安装完成之后，你就可以在命令行直接执行 ``ilus`` 这个命令了，如果你能顺利执行并看到类似如下的内容，那么就说明你已经安装成功了。


.. code:: bash

    $ ilus
    usage: ilus [-h] {WGS,genotype-joint-calling,VQSR} ...
    ilus: error: too few arguments


快速开始
-------

目前，ilus 包含了如下三个功能模块：

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

以下，分别对这三个功能的使用进行说明。

全基因组数据分析
--------------

全基因组数据分析流程通过 ``ilus WGS`` 来生成，用法如下：

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
      -O OUTDIR, --outdir OUTDIR
                            A directory for output results.


在 WGS 功能中，只有 ``-C``, ``-L`` 和 ``-O`` 这三个参数是必须的，其它的参数按照需要进行选择即可。其中，``-O`` 参数比较简单，就是项目的输出目录，该目录如果不存在，那么 **ilus** 会为你自动新建一个目录。最重要的是 ``-C`` 和 ``-L`` 参数，前者是 **ilus** 的配置文件，没有这个文件，**ilus** 就无法生成正确的流程，因此十分重要，后者是输入文件的列表文件，该列表文件一共有 5 列，每一列都是必须的信息，以下分别对这两个参数的格式进行说明：

首先是配置文件，我们需要在其中指定 WGS 流程各个步骤中所用的程序的路径以及所使用到GATK bundle文件和参考序列的路径。

需要注意的是 BWA MEM 的索引文件前缀需要与配置文件的 {resources}{reference} 相同，并存放在同一个目录中。如下：

.. code:: bash

    /path/human_reference/GRCh38/
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.dict
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.amb
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.ann
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.bwt
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.pac
    |-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.sa


该配置文件使用 Yaml 语法进行编写，在这里我提供一份该配置文件的例子，如下：

.. code:: bash

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

      # Adapter sequencing of BGISEQ-500
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

      interval: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                 "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
      
      # Specific variant calling interval, this could be the same as ``interval`` above
      variant_calling_interval: ["./wgs_calling_regions.GRCh38.interval.bed"]

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


在该配置文件，对于 WGS 流程来说所必须的生信软件是：bwa、samtools，bcftools、bedtools、gatk、bgzip和tabix，文件中的 ``verifyBamID2`` 参数仅用于计算样本是否存在污染，并不是必填的参数，如果配置文件中没有这个参数，那么流程则不进行样本污染情况的计算。另外，所必须的数据则是：gatk bundle 和参考序列。


接着是 ``-L`` 输入文件，这是分析流程所必须的所有测序数据，各列的信息如下：

* [1] Sample ID 样本名
* [2] Read Group，使用bwa mem时通过 -R 参数指定的 read group)
* [3] Fastq1 路径
* [4] Fastq2 路径，如果是Single End测序，没有fastq2，则该列用空格代替
* [5] fastq 的 lane 编号

对于测序量比较大，含有多个 lane 数据的样本，或者同一个 lane 的数据被拆分成了多个，不需要人工合并，只需要依照如上信息按行编写好输入文件即可，同一个样本的数据在流程中会在适当的时候由程序自动对其进行合并。如下是一个输入文件的例子：

.. code:: bash

    #SAMPLE RGID    FASTQ1  FASTQ2  LANE
    HG002   "@RG\tID:CL100076190_L01\tPL:COMPLETE\tPU:CL100076190_L01_HG002\tLB:CL100076190_L01\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L01_read_2.clean.fq.gz  CL100076190_L01
    HG002   "@RG\tID:CL100076190_L02\tPL:COMPLETE\tPU:CL100076190_L02_HG002\tLB:CL100076190_L02\tSM:HG002"  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_1.clean.fq.gz  /path/HG002_NA24385_son/BGISEQ500/BGISEQ500_PCRfree_NA24385_CL100076190_L02_read_2.clean.fq.gz  CL100076190_L02
    HG003   "@RG\tID:CL100076246_L01\tPL:COMPLETE\tPU:CL100076246_L01_HG003\tLB:CL100076246_L01\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L01_read_2.clean.fq.gz   CL100076246_L01
    HG003   "@RG\tID:CL100076246_L02\tPL:COMPLETE\tPU:CL100076246_L02_HG003\tLB:CL100076246_L02\tSM:HG003"  /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_1.clean.fq.gz   /path/HG003_NA24149_father/BGISEQ500/BGISEQ500_PCRfree_NA24149_CL100076246_L02_read_2.clean.fq.gz   CL100076246_L02
    HG004   "@RG\tID:CL100076266_L01\tPL:COMPLETE\tPU:CL100076266_L01_HG004\tLB:CL100076266_L01\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L01_read_2.clean.fq.gz   CL100076266_L01
    HG004   "@RG\tID:CL100076266_L02\tPL:COMPLETE\tPU:CL100076266_L02_HG004\tLB:CL100076266_L02\tSM:HG004"  /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_1.clean.fq.gz   /path/HG004_NA24143_mother/BGISEQ500/BGISEQ500_PCRfree_NA24143_CL100076266_L02_read_2.clean.fq.gz   CL100076266_L02
    HG005   "@RG\tID:CL100076244_L01\tPL:COMPLETE\tPU:CL100076244_L01_HG005\tLB:CL100076244_L01\tSM:HG005"  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_1.clean.fq.gz  /path/HG005_NA24631_son/BGISEQ500/BGISEQ500_PCRfree_NA24631_CL100076244_L01_read_2.clean.fq.gz  CL100076244_L01

以下是一些使用 **ilus** 生成 WGS 分析流程的例子。


**例子1：从头开始执行 WGS 流程**

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -O ./output

这个命令的意思是，依据 ``ilus_sys.yaml`` 和 ``input.list`` 在输出目录 ``output`` 生成名为 （-n）``my_wgs`` 的 WGS 分析流程，并将最后的比对数据从 BAM 转为 CRAM (-c)。输出目录 ``output`` 有 4 个文件夹（如下），用于存放由 WGS 分析流程产生的各类数据。

.. code:: bash
    
    00.shell/
    01.alignment/
    02.gvcf/
    03.genotype/

从文件夹的名字，我们也可以了解到各个目录的具体作用。比如 ``00.shell`` 目录存放的是流程各个步骤的执行脚本和日志文件的目录：

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


我们依照从 step1 到 step6执行流程即可。loginfo目录记录了各个步骤各个样本的运行状态，我们可以检查各个步骤的 .o.log.list 日志文件，获得该样本是否成功结束的标记。如果成功结束了，那么在该日志文件的末尾会有一个 ``**[xx] xxxx done**`` 的标记。可以通过使用 **ilus** 提供的脚本 ``check_jobs_status.py`` 检查各个步骤是否已经全部顺利完成，如果有错那么该脚本会将未完成的任务输出，方便我们重新执行。用法为：

.. code:: bash

    $ python check_jobs_status.py loginfo/01.alignment.o.log.list > bwa.unfinish.list

如果任务都是成功结束的，那么该 list 文件为空，并输出 ``** All Jobs done **``。

**例子2：只执行 WGS 流程中某个/某些步骤**

有时候，我们并打算从头到尾完整地将 WGS 流程执行下去，比如我们只想执行从 fastq 比对到生成 gvcf 这个步骤，暂时不想执行 genotype joint-calling 和 VQSR，那么这个时候我们可以通过 ``-P`` 参数指定特定的步骤：

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P align,markdup,BQSR,gvcf -O ./output


这样就只会生成从 bwa 到 gvcf 的 shell 脚本。

除此之外，当你发现 WGS 的某个步骤跑错了，需要重新更新时，你也可以用 ``-P`` 指定重跑特定的步骤。比如我想重生成 BQSR 这个步骤的运行脚本，那么就可以这样做：

.. code:: bash

    $ ilus WGS -c -n my_wgs -C ilus_sys.yaml -L input.list -P BQSR -O ./output

需要注意的是，**ilus** 为了节省项目的空间，只会为每一个样本保留 BQSR 之后的 BAM/CRAM 文件，因此，如果你想重新跑 BQSR 需要确定在 BQSR 前一步（即，markdup）的 BAM 文件是否已经被删除了，如果原先 **ilus** 在BQSR这一步没有正常结束的话，那么该 markdup 的 BAM 文件应该还会被保留着的，**ilus** 执行任务时具有“原子属性”，也就是说只有当所有步骤都成功结束时才会删除在之后的分析中完全不需要的文件。


genotype-joint-calling
----------------------

如果我们已经有了各个样本的 gvcf 需要从这些 gvcf 开始完成多样本的联合变异检测（Joint-calling），那么就可以使用 ``genotype-joint-calling`` 来实现。具体用法如下：

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


在 **ilus genotype-joint-calling** 中输入的 gvcf list 文件，由两列构成，第一列是该 gvcf 所在的区间或者染色体编号，第二列是该 gvcf 文件的路径，举个例子：

.. code:: bash

    $ ilus genotype-joint-calling -n my_project -C ilus_sys.yaml -L gvcf.list -O 03.genotype --as_pipe_shell_order

其中 ``gvcf.list`` 的格式类似如下：

.. code:: bash

    chr1    /path/sample1.chr1.g.vcf.gz
    chr1    /paht/sample2.chr1.g.vcf.gz
    chr2    /path/sample1.chr2.g.vcf.gz
    chr2    /path/sample2.chr2.g.vcf.gz
    ...
    chrM    /path/sample1.chrM.g.vcf.gz
    chrM    /path/sample2.chrM.g.vcf.gz

以上假设 gvcf.list 中只有两个样本。

参数 ``--as_pipe_shell_order`` 可加也可不加（默认是不加），它唯一的作用就是按照 **ilus WGS** 流程的方式输出执行脚本的名字，维持和 WGS 流程一样的次序。


VQSR
----

该功能仅用于生成基于 ``GATK VQSR`` 的执行脚本。我们如果已经有了 VCF 结果，现在只想单独对这个变异数据跑 VQSR 进行初步的质控，那么就可以使用这个模块，用法与 ``genotype-joint-calling`` 大同小异，如下：

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

**ilus VQSR** 的其它参数与 ``genotype-joint-calling`` 相同，以下为一个完整的例子：

.. code:: bash

    $ ilus VQSR -C ilus_sys.yaml -L vcf.list -O 03.genotype --as_pipe_shell_order


