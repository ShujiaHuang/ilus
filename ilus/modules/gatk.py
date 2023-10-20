"""GATK modules

Author: Shujia Huang
Date: 2020-04-19 15:19:56
"""


def markduplicates(config, input_bam, output_markdup_bam, out_metrics_fname):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["markdup_java_options"]) \
        if "markdup_java_options" in config["gatk"] \
           and len(config["gatk"]["markdup_java_options"]) else ""

    return (f"time {gatk} {java_options} MarkDuplicates "
            f"-I {input_bam} "
            f"-M {out_metrics_fname} "
            f"-O {output_markdup_bam}")


def baserecalibrator(config, input_bam, output_bqsr_bam, out_bqsr_recal_table):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["bqsr_java_options"]) \
        if "bqsr_java_options" in config["gatk"] \
           and len(config["gatk"]["bqsr_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta
    known_site_1000G_indel = config["gatk"]["bundle"]["1000G_known_indel"]
    known_site_mills_gold_indels = config["gatk"]["bundle"]["mills"]
    known_site_dbsnp = config["gatk"]["bundle"]["dbsnp"]

    # create recalibrate table file for BQSR
    recal_data_cmd = (f"time {gatk} {java_options} BaseRecalibrator "
                      f"-R {reference} "
                      f"--known-sites {known_site_1000G_indel} "
                      f"--known-sites {known_site_mills_gold_indels} "
                      f"--known-sites {known_site_dbsnp} "
                      f"-I {input_bam} "
                      f"-O {out_bqsr_recal_table}")

    # ApplyBQSR
    apply_bqsr_cmd = (f"time {gatk} {java_options} ApplyBQSR "
                      f"-R {reference} "
                      f"--bqsr-recal-file {out_bqsr_recal_table} "
                      f"-I {input_bam} "
                      f"-O {output_bqsr_bam}")

    return recal_data_cmd + " && " + apply_bqsr_cmd


def haplotypecaller_gvcf(config, input_bam, output_gvcf_fname, interval=None):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["hc_gvcf_java_options"]) \
        if "hc_gvcf_java_options" in config["gatk"] \
           and len(config["gatk"]["hc_gvcf_java_options"]) else ""
    hc_options = " ".join(config["gatk"]["hc_gvcf_options"]) if "hc_gvcf_options" in config["gatk"] else ""

    reference = config["resources"]["reference"]  # reference fasta
    cmd = (f"time {gatk} {java_options} HaplotypeCaller "
           f"-R {reference} {hc_options} "
           f"--emit-ref-confidence GVCF "
           f"-I {input_bam} "
           f"-O {output_gvcf_fname}")

    if interval:
        cmd += f" -L {interval}"

    return cmd


def combineGVCFs(config, input_sample_gvcfs, output_combineGVCF_fname, interval=None):
    """Combine GVCFs by GATK genomicsDBImport or CombineGVCFs module.
    """
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["combineGVCFs_java_options"]) \
        if "combineGVCFs_java_options" in config["gatk"] \
           and len(config["gatk"]["combineGVCFs_java_options"]) else ""

    # set overwite existing genomicsdb workspace by default
    if "genomicsDBImport_options" not in config["gatk"]:
        config["gatk"]["genomicsDBImport_options"] = ["--overwrite-existing-genomicsdb-workspace true"]

    elif (("--overwrite-existing-genomicsdb-workspace false" not in config["gatk"]["genomicsDBImport_options"]) and
          ("--overwrite-existing-genomicsdb-workspace true" not in config["gatk"]["genomicsDBImport_options"])):
        config["gatk"]["genomicsDBImport_options"].append("--overwrite-existing-genomicsdb-workspace true")

    genomicsDBImport_options = "%s" % " ".join(config["gatk"]["genomicsDBImport_options"])
    reference = config["resources"]["reference"]  # reference fasta
    use_gDBI = config["gatk"]["use_genomicsDBImport"] if "use_genomicsDBImport" in config["gatk"] else False

    # Create command line for GenomicsDBImport or CombineGVCFs
    if use_gDBI:
        # use GenomicsDBImport
        sample_name_map = input_sample_gvcfs[0]  # Only one file
        combine_gvcf_cmd = (f"time {gatk} {java_options} GenomicsDBImport {genomicsDBImport_options} "
                            f"-R {reference} "
                            f"--sample-name-map {sample_name_map} "
                            f"--genomicsdb-workspace-path {output_combineGVCF_fname}")
    else:
        gvcfs = " ".join(["-V %s" % s for s in input_sample_gvcfs])
        combine_gvcf_cmd = (f"time {gatk} {java_options} CombineGVCFs "
                            f"-R {reference} {gvcfs} "
                            f"-O {output_combineGVCF_fname}")
    if interval:
        combine_gvcf_cmd += f" -L {interval}"

    return combine_gvcf_cmd


def genotypeGVCFs(config, input_combine_gvcf_fname, output_vcf_fname, interval=None):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["genotype_java_options"]) \
        if "genotype_java_options" in config["gatk"] \
           and len(config["gatk"]["genotype_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta
    use_gDBI = config["gatk"]["use_genomicsDBImport"] if "use_genomicsDBImport" in config["gatk"] else False

    genotypeGVCFs_options = " ".join(config["gatk"]["genotypeGVCFs_options"]) \
        if "genotypeGVCFs_options" in config["gatk"] else ""

    # Creat command line for genotypeGVCF
    if use_gDBI:
        genotype_cmd = (f"time {gatk} {java_options} GenotypeGVCFs "
                        f"-R {reference} {genotypeGVCFs_options} "
                        f"-V gendb://{input_combine_gvcf_fname} "
                        f"-O {output_vcf_fname}")
    else:
        genotype_cmd = (f"time {gatk} {java_options} GenotypeGVCFs "
                        f"-R {reference} {genotypeGVCFs_options} "
                        f"-V {input_combine_gvcf_fname} "
                        f"-O {output_vcf_fname}")
    if interval:
        genotype_cmd += f" -L {interval}"

    return genotype_cmd


def variantrecalibrator(config, input_vcf, output_vcf_fname):
    gatk = config["gatk"]["gatk"]
    java_options = "--java-options \"%s\"" % " ".join(config["gatk"]["vqsr_java_options"]) \
        if "vqsr_java_options" in config["gatk"] \
           and len(config["gatk"]["vqsr_java_options"]) else ""

    reference = config["resources"]["reference"]  # reference fasta
    vqsr_options = " ".join(config["gatk"]["vqsr_options"]) if "vqsr_options" in config["gatk"] else ""

    # Set name
    out_prefix = output_vcf_fname.replace(".gz", "").replace(".vcf", "")  # delete .vcf.gz
    out_snp_vqsr_fname = out_prefix + ".SNPs.vcf.gz"

    resource_hapmap = config["gatk"]["bundle"]["hapmap"]
    resource_omni = config["gatk"]["bundle"]["omni"]
    resource_1000G = config["gatk"]["bundle"]["1000G"]
    resource_dbsnp = config["gatk"]["bundle"]["dbsnp"]
    resource_mills_gold_indels = config["gatk"]["bundle"]["mills"]

    # SNP VQSR
    snp_vqsr_cmd = (f"time {gatk} {java_options} VariantRecalibrator "
                    f"-R {reference} "
                    f"-V {input_vcf} "
                    f"--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {resource_hapmap} "
                    f"--resource:omini,known=false,training=true,truth=false,prior=12.0 {resource_omni} "
                    f"--resource:1000G,known=false,training=true,truth=false,prior=10.0 {resource_1000G} "
                    f"--resource:dbsnp,known=true,training=false,truth=false,prior=5.0 {resource_dbsnp} "
                    f"{vqsr_options} "
                    f"-mode SNP "
                    f"--tranches-file {out_prefix}.SNPs.tranches "
                    f"-O {out_prefix}.SNPs.recal")

    apply_snp_vqsr_cmd = (f"time {gatk} {java_options} ApplyVQSR "
                          f"-R {reference} "
                          f"-V {input_vcf} "
                          f"--tranches-file {out_prefix}.SNPs.tranches "
                          f"--recal-file {out_prefix}.SNPs.recal "
                          f"--truth-sensitivity-filter-level 99.0 "
                          f"-mode SNP "
                          f"-O {out_snp_vqsr_fname}")

    # Indel VQSR after SNP
    indel_vqsr_cmd = (f"time {gatk} {java_options} VariantRecalibrator "
                      f"-R {reference} "
                      f"-V {out_snp_vqsr_fname} "
                      f"--resource:mills,known=true,training=true,truth=true,prior=12.0 {resource_mills_gold_indels} "
                      f"--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {resource_dbsnp} "
                      f"{vqsr_options} "
                      f"--tranches-file {out_prefix}.INDELs.tranches "
                      f"-mode INDEL "
                      f"-O {out_prefix}.INDELs.recal")
    apply_indel_vqsr_cmd = (f"time {gatk} {java_options} ApplyVQSR "
                            f"-R {reference} "
                            f"-V {out_snp_vqsr_fname} "
                            f"--truth-sensitivity-filter-level 99.0 "
                            f"--tranches-file {out_prefix}.INDELs.tranches "
                            f"--recal-file {out_prefix}.INDELs.recal "
                            f"-mode INDEL "
                            f"-O {output_vcf_fname} && rm -f {out_snp_vqsr_fname}")

    return " && ".join([snp_vqsr_cmd, apply_snp_vqsr_cmd, indel_vqsr_cmd, apply_indel_vqsr_cmd])


def collect_alignment_summary_metrics(config, input_bam, output_fname):
    gatk = config["gatk"]["gatk"]
    reference = config["resources"]["reference"]  # reference fasta

    java_options = f'--java-options "{" ".join(config["gatk"]["CollectAlignmentSummaryMetrics_jave_options"])}"' \
        if "CollectAlignmentSummaryMetrics_jave_options" in config["gatk"] \
           and len(config["gatk"]["CollectAlignmentSummaryMetrics_jave_options"]) else ""

    metric_options = " ".join(config["gatk"]["CollectAlignmentSummaryMetrics_options"]) \
        if "CollectAlignmentSummaryMetrics_options" in config["gatk"] \
           and len(config["gatk"]["CollectAlignmentSummaryMetrics_options"]) else ""

    cmd = (f"time {gatk} {java_options} CollectAlignmentSummaryMetrics {metric_options} "
           f"-R {reference} "
           f"-I {input_bam} "
           f"-O {output_fname}")

    return cmd


def mergevcfs(config, input_vcfs, output_fname):
    # Sometimes bcftools concat is better than MergeVcfs
    gatk = config["gatk"]["gatk"]

    java_options = f'--java-options "{" ".join(config["gatk"]["mergevcfs_java_options"])}"' \
        if "mergevcfs_java_options" in config["gatk"] \
           and len(config["gatk"]["mergevcfs_java_options"]) else ""

    cmd = [f"time {gatk} {java_options} MergeVcfs -O {output_fname}"] + [f"-I {f}" for f in input_vcfs]

    return " ".join(cmd)
