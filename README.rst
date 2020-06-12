Ilus
====

**Ilus** is a handy NGS pipeline for whole genome re-sequencing (WGS) data analysis.

Installation
------------
Ilus is writed by Python and maintain in PyPI. You can install it just by running the following command:

.. code:: bash

    pip install ilus

You can type **ilus** in your commandline terminal if everything is smooth.

Quick Start
-----------

There are several functions in **ilus**:

.. code:: bash

    $ ilus --help

You will see all the functions:

.. code:: bash

    usage: ilus [-h] {WGS,genotype-joint-calling,VQSR} ...

    ilus: A WGS analysis pipeline.

    optional arguments:
        -h, --help            show this help message and exit

    ilus commands:
    {WGS,genotype-joint-calling,VQSR}
        WGS                 Creating pipeline for WGS(from fastq to genotype VCF)
        genotype-joint-calling Genotype from GVCFs.
        VQSR                VQSR

Run GATK Best Practices Workflow
--------------------------------

**ilus** implelement the whole WGS analysis workflow base on the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).



