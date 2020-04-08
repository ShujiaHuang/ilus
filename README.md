# Ilus WGS数据分析流程

## 比对: bwamem
该步骤同时完成了bwamem和结果排序。命令行参数为：

```bash
usage: python ilus_process.py bwamem [-h] [-R RG_INFO] [-L LANE] [-O OUTDIR]
                                 [-C SYSCONF]
                                 fastq1 fastq2

positional arguments:
  fastq1                input fastq file of read1 for SE or PE.
  fastq2                input fastq file of read2 for PE (optional).

optional arguments:
  -h, --help            show this help message and exit
  -R RG_INFO, --readgroup RG_INFO
                        read group header line such as '@RG\tID:foo\tSM:bar'
  -L LANE, --lane LANE  Sequencing lane ID.
  -O OUTDIR, --outdir OUTDIR
                        A directory for output results.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```
给出例子：

```bash
time python ../../scripts/ilus_process.py bwamem -C ../../config/ilus_sys.yaml -R "@RG\tID:H0164ALXX140820.2\tPL:illumina\tPU:H0164ALXX140820.2\tLB:Solexa-272222\tPI:0\tSM:NA12878\tCN:BI" -L H0164ALXX140820.2 -O outdir ../data/germline/fastqs/NA12878_1.fastq.gz ../data/germline/fastqs/NA12878_2.fastq.gz
```

这里需要注意的是 -R 参数和bwa自带的 -R 参数一模一样，都需要完整的Read Group描述信息，另一个注意点是 -O 参数要求的是输出目录，而不是文件，原因是输出比对文件的时候，程序会自己按照 -L 为文件命名，而且由于在输出比对结果时，同时对结果进行了排序，因此，在文件名中会自动加上 `.sorted`。

## Markduplicates

调用sambamba对比对结果进行重复序列标记。命令行参数：
```bash
usage: ilus_process.py markdups [-h] [-I INBAM]
                                   [--REMOVE_DUPLICATES REMOVEDUPS]
                                   [-O OUTFILE] [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -I INBAM, --input INBAM
                        Input one or more BAM files to analyze. This argument
                        must be specified at least once. Reuqired
  --REMOVE_DUPLICATES REMOVEDUPS
                        If true do not write duplicates to the output file
                        instead of writing them with appropriate flags set.
                        [False]
  -O OUTFILE, --outfile OUTFILE
                        A BAM file's name for output results.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

操作例子：

```bash
time python ../../scripts/ilus_process.py markdups -C ../../config/ilus_sys.yaml -I out/H0164ALXX140820.2.sorted.bam -O out/H0164ALXX140820.2.sorted.markdups.bam
```

如果`--REMOVE_DUPLICATES`参数设置为`True`，则表示将删除重复标记的序列。


## BQSR
调用GATK4完成BQSR，命令行参数：

```bash
usage: ilus_process.py BQSR [-h] [-I INBAM] [-O OUTFILE] [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -I INBAM, --input INBAM
                        Input one or more BAM files to analyze. This argument
                        must be specified at least once. Reuqired
  -O OUTFILE, --outfile OUTFILE
                        A BAM file's name for output results.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

BQSR例子：
```bash
time python ../../scripts/ilus_process.py BQSR -C ../../config/ilus_sys.yaml -I out/H0164ALXX140820.2.sorted.markdups.bam -O out/H0164ALXX140820.2.sorted.markdups.BQSR.bam
```

这里会同时输出BQSR之后的BAM文件和该文件的索引文件，但是索引文件的后缀是替换了.bam的格式，如：`H0164ALXX140820.2.sorted.markdups.BQSR.bai` 而不是 `H0164ALXX140820.2.sorted.markdups.BQSR.bam.bai`这和sambamba以及samtools在给出BAM index的时候是不同的。


## HaplotypeCaller
通过调用GATK4 HaplotypeCaller进行变异检测，这里一般先完成gvcf，之后再进行genotype，参数如下：

```bash
usage: ilus_process.py HaplotypeCaller [-h] [-I INBAM] [-O OUTFILE]
                                          [-L INTERVAL] [-E EMIT_GVCF]
                                          [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -I INBAM, --input INBAM
                        A VCF file containing variants. Reuqired
  -O OUTFILE, --outfile OUTFILE
                        File to which variants should be written.
  -L INTERVAL, --intervals INTERVAL
                        One genomic intervals over which to operate.
  -E EMIT_GVCF, --emit-ref-confidence EMIT_GVCF
                        Emitting gvcf.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

HC为每个样本生成gvcf的例子：
```bash
time python ../../scripts/ilus_process.py HaplotypeCaller -C ../../config/ilus_sys.yaml -I out/H0164ALXX140820.2.sorted.markdups.BQSR.bam -E True -O out/NA12878.g.vcf.gz
```

`-E`参数的作用就是HaplotypeCaller的`--emit-ref-confidence GVCF`，如果为`True`那么生成g.vcf，否则就直接完成变异检测和Genotype。如：

```bash
time python ../../scripts/ilus_process.py HaplotypeCaller -C ../../config/ilus_sys.yaml -I out/H0164ALXX140820.2.sorted.markdups.BQSR.bam -O out/NA12878.vcf.gz
```


## GenotypeGVCFs
这是生成g.vcf，然后再进行 **单样本**或者 **多样本**joint-calling时候才需要调用的功能，命令行参数如下：

```bash
usage: ilus_process.py GenotypeGVCFs [-h] [-V VARIANTS] [-L INTERVAL]
                                        [-O OUTFILE] [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -V VARIANTS, --variant VARIANTS
                        Input one or more BAM files to analyze. This argument
                        must be specified at least once. Reuqired
  -L INTERVAL, --intervals INTERVAL
                        One genomic intervals over which to operate.
  -O OUTFILE, --outfile OUTFILE
                        File to which variants should be written.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

例子：
```bash
time python ../../scripts/ilus_process.py GenotypeGVCFs -C ../../config/ilus_sys.yaml -V out/NA12878.g.vcf.gz -O out/NA12878.vcf.gz
```
除了生成最后的VCF之外，还会同时为该vcf构建索引，因此我们都不需要额外进行索引的生成。

## VQSR
使用GATK4进行VQSR，命令行参数：

```bash
usage: ilus_process.py VQSR [-h] [-V VARIANTS]
                               [--snp-ts-filter-level SNP_TS_FILTER_LEVEL]
                               [--indel-ts-filter-level INDEL_TS_FILTER_LEVEL]
                               [-O OUTFILE] [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -V VARIANTS, --variant VARIANTS
                        Input one or more BAM files to analyze. This argument
                        must be specified at least once. Reuqired
  --snp-ts-filter-level SNP_TS_FILTER_LEVEL
                        The truth sensitivity level for SNP at which to start
                        filtering. Default: 99.0
  --indel-ts-filter-level INDEL_TS_FILTER_LEVEL
                        The truth sensitivity level for Indel at which to
                        start filtering. Default: 95.0
  -O OUTFILE, --outfile OUTFILE
                        File to which variants should be written.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

例子：

```bash
time python ../../scripts/ilus_process.py VQSR -C ../../config/ilus_sys.yaml -V out/NA12878.vcf.gz -O NA12878.VQSR.vcf.gz
```

不过在本例子中我们的数据很小，而VQSR要求要有足够的数据量，才能完成质控模型的训练，所以这个步骤其实是不能成功的，需要用更大的数据量（或者全基因组的结果）才行。

## 其它功能

### mergebam

Ilus的mergebam模块可以利用sambamba合并多个不同的BAM文件，命令行参数：

```bash
usage: python ilus_process.py mergebam [-h] [-I INBAM] [-O OUTFILE] [-C SYSCONF]

optional arguments:
  -h, --help            show this help message and exit
  -I INBAM, --input INBAM
                        Input one or more BAM files to analyze. This argument
                        must be specified at least once. Reuqired
  -O OUTFILE, --outfile OUTFILE
                        A BAM file's name for output results.
  -C SYSCONF, --conf SYSCONF
                        YAML configuration file specifying details about
                        system.
```

例子：

```bash
ime python ../../scripts/ilus_process.py mergebam -C ../../config/ilus_sys.yaml -I out/HXXXL1.chr1.sorted.bam -I out/HXXXL1.chr2.sorted.markdups.bam -O out/merge.bam
```


### fastq-splitter

同时支持gz和非gz格式，split之后的文件命令规则为：`fqname+序号+后缀`。其中`$fqname`将根据输入的fastq文件名进行自动推断，后缀(`$ext`)默认是gz格式，除非用`--uncompress`指定不压缩，如下：

```perl
    my ($fqname, $ext);
    if ($infile =~ m/\.fastq.gz$/) {
        $fqname = basename($infile, ".fastq.gz");
        $ext = "fastq";
    } elsif ($infile =~ m/\.fastq$/) {
        $fqname = basename($infile, ".fastq");
        $ext = "fastq";
    } elsif ($infile =~ m/\.fq.gz$/) {
        $fqname = basename($infile, ".fq.gz");
        $ext = "fq";
    } elsif ($infile =~ m/\.fq$/) {
        $fqname = basename($infile, ".fq");
        $ext = "fq";
    } else {
        die "Error: Input files are not fastq files(*.fq.gz/*.fq/*.fastq.gz/*.fastq)\n";
    }
    $ext .= ".gz" if !$uncompress;
```

具体用法如下：

```bash
time perl fastq-splitter.pl --n-parts 8 --outdir ./ SRR1770413.1.fq
time perl fastq-splitter.pl --n-parts 8 --outdir ./ SRR1770413.1.fq.gz

time perl fastq-splitter.pl --n-parts 8 --outdir ./ SRR1770413.1.fastq
time perl fastq-splitter.pl --n-parts 8 --outdir ./ SRR1770413.1.fastq.gz
```

文件名自动根据设定的split文件数目确定，比如：

```bash
$ perl fastq-splitter.pl --n-parts 20 --outdir ./ SRR1770413.1.fastq
```

得到的20份输出文件列表如下：
```
SRR1770413.1.01.fq
SRR1770413.1.02.fq
SRR1770413.1.03.fq
SRR1770413.1.04.fq
SRR1770413.1.05.fq
SRR1770413.1.06.fq
SRR1770413.1.07.fq
SRR1770413.1.08.fq
SRR1770413.1.09.fq
SRR1770413.1.10.fq
SRR1770413.1.11.fq
SRR1770413.1.12.fq
SRR1770413.1.13.fq
SRR1770413.1.14.fq
SRR1770413.1.15.fq
SRR1770413.1.16.fq
SRR1770413.1.17.fq
SRR1770413.1.18.fq
SRR1770413.1.19.fq
SRR1770413.1.20.fq
```

fastq-splitter支持同时对多份fastq进行split，如下例子，同时对两份fastq文件进行split之后，分别按照文件名的规则输出到`subdir`目录下。

```
$ perl fastq-splitter.pl --n-parts 100 --outdir ./subdir SRR1770413_1.fastq.gz SRR1770413_2.fastq.gz
```




