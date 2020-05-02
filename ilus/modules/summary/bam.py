"""Summary statistic information for BAM/CRAM

Author: Shujia Huang
Date: 2020-04-29 14:50:25
"""


def stats(config, input_bamfile, output_bamstats_fname):
    samtools = config["samtools"]["samtools"]
    stats_options = "%s" % " ".join(config["samtools"]["stats_options"]) \
        if "stats_options" in config["samtools"] and \
           len(config["samtools"]["stats_options"]) else ""

    cmd = "time {samtools} stats {stats_options} {input_bamfile} > {output_bamstats_fname}".format(**locals())
    return cmd


def genomecoverage(config, input_bamfile, output_cvg_fname):
    tabix = config["tabix"]
    bgzip = config["bgzip"]
    bedtools = config["bedtools"]["bedtools"]
    genomecov_options = "%s" % " ".join(config["bedtools"]["genomecov_options"]) \
        if "genomecov_options" in config["bedtools"] and \
           len(config["bedtools"]["genomecov_options"]) else ""

    cmd = ("time {bedtools} genomecov {genomecov_options} "
           "-ibam {input_bamfile} "
           "| {bgzip} > {output_cvg_fname} "
           "&& {tabix} -f -p bed {output_cvg_fname}").format(**locals())

    return cmd
