"""A bedtools modules.

Author: Shujia Huang
Date: 2020-04-28 20:29:09
"""


def concat(config, input_vcfs, output_vcf):

    tabix = config["tabix"]
    bedtools = config["bedtools"]["bedtools"]
    concat_options = config["bedtools"]["concat_options"] \
        if "concat_options" in config["bedtools"] else []

    concat_options = " ".join(concat_options)

    cmd = [("time {bedtools} concat {concat_options} "
            "-O z "
            "-o {output_vcf}").format(**locals())] + input_vcfs

    return " ".join(cmd) + " && time {tabix} -p vcf {output_vcf}".format(**locals())
