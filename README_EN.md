Ilus
====

**Ilus** is a handy pipeline generater for whole genome re-sequencing
(WGS) and whole exom sequencing data (WES) analysis.

English \| [简体中文](./README.rst)

Installation
------------

**ilus** is written by Python and release in PyPI. You can install it
just by running the following command:

``` {.bash}
pip install ilus
```

You can type **ilus** in your commandline terminal if everything is
smooth.

Quick Start
-----------

There are several functions in **ilus**:

``` {.bash}
$ ilus --help
```

You will see all the functions:

``` {.bash}
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

Here is all the Chinese [简体中文帮助文档](./README_CN.rst) for how to
creat WGS/WES pipeline by ilus.

English version is coming.
