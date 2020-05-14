"""A bcftools modules.

Author: Shujia Huang
Date: 2020-04-28 20:29:09
"""


def concat(config, input_vcfs, output_vcf):

    tabix = config["tabix"]
    bcftools = config["bcftools"]["bcftools"]
    concat_options = "%s" % " ".join(config["bcftools"]["concat_options"]).replace("-O z", "") \
        if "concat_options" in config["bcftools"] and \
           len(config["bcftools"]["concat_options"]) else ""

    # Always output gz compress
    if not output_vcf.endswith(".gz"):
        output_vcf += ".gz"

    cmd = [("time {bcftools} concat {concat_options} "
            "-O z "
            "-o {output_vcf}").format(**locals())] + input_vcfs

    return " ".join(cmd) + " && time {tabix} -p vcf {output_vcf}".format(**locals())
