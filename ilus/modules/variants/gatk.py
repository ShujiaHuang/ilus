"""GATK modules

Author: Shujia Huang
Date: 2020-04-19 15:19:56
"""
import os


def markduplicates(config, input_bam, output_markdup_bam, out_metrics_fname):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % config["gatk"]["markdup_java_options"] \
        if "markdup_java_options" in config["gatk"] \
           and len(config["gatk"]["markdup_java_options"]) else ""

    return ("time {gatk} {java_options} MarkDuplicates "
            "-I {input_bam} "
            "-M {out_metrics_fname} "
            "-O {output_markdup_bam}").format(**locals())


def baserecalibrator(config, input_bam, output_bqsr_bam, out_bqsr_recal_table):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % config["gatk"]["bqsr_java_options"] \
        if "bqsr_java_options" in config["gatk"] \
           and len(config["gatk"]["bqsr_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta
    known_site_1000G_indel = config["gatk"]["bundle"]["1000G_known_indel"]
    known_site_mills_gold_indels = config["gatk"]["bundle"]["mills"]
    known_site_dbsnp = config["gatk"]["bundle"]["dbsnp"]

    # create recalibrate table file for BQSR
    recal_data_cmd = ("time {gatk} {java_options} BaseRecalibrator "
                      "-R {reference} "
                      "--known-sites {known_site_1000G_indel} "
                      "--known-sites {known_site_mills_gold_indels} "
                      "--known-sites {known_site_dbsnp} "
                      "-I {input_bam} "
                      "-O {out_bqsr_recal_table}").format(**locals())

    # ApplyBQSR
    apply_bqsr_cmd = ("time {gatk} {java_options} ApplyBQSR "
                      "-R {reference} "
                      "--bqsr-recal-file {out_bqsr_recal_table} "
                      "-I {input_bam} "
                      "-O {output_bqsr_bam}").format(**locals())

    return recal_data_cmd + " && " + apply_bqsr_cmd


def haplotypecaller_gvcf(config, input_bam, output_gvcf_fname, interval=None):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % config["gatk"]["hc_gvcf_java_options"] \
        if "hc_gvcf_java_options" in config["gatk"] \
           and len(config["gatk"]["hc_gvcf_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta
    cmd = ("time {gatk} {java_options} HaplotypeCaller "
           "-R {reference} "
           "--emit-ref-confidence GVCF "
           "-I {input_bam} "
           "-O {output_gvcf_fname}").format(**locals())

    if interval:
        cmd += " -L %s" % interval

    return cmd


def genotypegvcfs(config, input_sample_gvcfs, output_vcf_fname):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % config["gatk"]["genotype_java_options"] \
        if "genotype_java_options" in config["gatk"] \
           and len(config["gatk"]["genotype_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta

    # create a combine gvcf
    directory, fname = os.path.split(output_vcf_fname)
    combine_gvcf_fname = os.path.join(directory, fname.replace(".vcf", ".g.vcf"))

    genotype_cmd = []
    if len(input_sample_gvcfs) > 1:
        gvcfs = " ".join(["-V %s" % s for s in input_sample_gvcfs])
        combine_gvcf_cmd = ("time {gatk} {java_options} CombineGVCFs "
                            "-R {reference} {gvcfs} "
                            "-O {combine_gvcf_fname}").format(**locals())
        genotype_cmd = [combine_gvcf_cmd]
    else:
        combine_gvcf_fname = input_sample_gvcfs[0]

    genotype_cmd.append(("time {gatk} {java_options} GenotypeGVCFs "
                         "-R {reference} "
                         "-V {combine_gvcf_fname} "
                         "-O {output_vcf_fname}").format(**locals()))

    genotype_cmd.append("rm -rf %s" % combine_gvcf_fname)
    return " && ".join(genotype_cmd)


def variantrecalibrator(config, input_vcf, output_vcf_fname):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % config["gatk"]["vqsr_java_options"] \
        if "vqsr_java_options" in config["gatk"] \
           and len(config["gatk"]["vqsr_java_options"]) else ""

    vqsr_cmd = ""

    return
