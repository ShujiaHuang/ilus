"""GATK modules

Author: Shujia Huang
Date: 2020-04-19 15:19:56
"""


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

    samtools = config["samtools"]["samtools"]
    bam_index_cmd = "time {samtools} index {output_bqsr_bam}".format(**locals())

    return recal_data_cmd + " && " + apply_bqsr_cmd + " && " + bam_index_cmd


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
