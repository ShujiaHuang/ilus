"""
"""
from _gatk4 import *


def gatk4_mark_duplicates(gatk, align_bams, outbamfile, remove_dups=False):
    """Run GATK MarkDuplicates.

    Parameters:
        ``gatk``: GATKRunner
        ``align_bams``: A list like
            contain at least one BAM files
    """
    bamname, ext = os.path.splitext(outbamfile)
    dup_metrics = "%s.metrics.txt" % bamname
    options = [("-I", in_bam) for in_bam in align_bams] + \
              [("-M", dup_metrics), ("-O", outbamfile),
               ("--REMOVE_DUPLICATES", "true" if remove_dups else "false")]

    gatk.run("MarkDuplicates", options)
    return outbamfile, dup_metrics


def gatk4_baserecalibrator(gatk, align_bams, outbamfile, config):
    """Run GATK BaseRecalibrator.
    """
    bamname, ext = os.path.splitext(outbamfile)
    bam_bai_name = "%s.bai" % bamname
    recal_table = "%s.recal_data.table" % bamname

    known_sites = [config["gatk_bundle"][ks]
                   for ks in ["1000G_phase1_indel", "mills", "dbsnp"]
                   if ks in config["gatk_bundle"]]

    opt_bqsr = [("-R", gatk.reference), ("-O", recal_table)] + \
               [("-I", in_bam) for in_bam in align_bams] + \
               [("--known-sites", ks) for ks in known_sites]

    gatk.run("BaseRecalibrator", opt_bqsr)

    # ApplyBQSR
    opt_applybqsr = [("-R", gatk.reference), ("--bqsr-recal-file", recal_table)] + \
                    [("-I", in_bam) for in_bam in align_bams] + \
                    [("-O", outbamfile)]

    gatk.run("ApplyBQSR", opt_applybqsr)
    return outbamfile, bam_bai_name, recal_table


def gatk4_haplotypecaller(gatk, align_bams, intervals, outvcffile, emit_gvcf=False):
    """Run HaplotypeCaller
    """
    if (not outvcffile.endswith(".vcf.gz")) and (not outvcffile.endswith(".vcf")):
        sys.stderr.write("Error: %s is not end with .vcf.gz or .vcf in "
                         "gatk4_haplotypecaller.\n" % outvcffile)
        sys.exit(1)

    options = [("-R", gatk.reference), ("-O", outvcffile)] + \
              [("-I", in_bam) for in_bam in align_bams]
    if intervals:
        options += [("-L", reg) for reg in intervals]

    if emit_gvcf:
        options.append(("--emit-ref-confidence", "GVCF"))

    gatk.run("HaplotypeCaller", options)
    return outvcffile


def gatk4_IndexFeatureFile(gatk, in_file):
    """Creates an index for a feature file, e.g. VCF or BED file.
    And the index file will be in the same directory as the input file.
    """

    options = [("--feature-file", in_file)]
    if in_file.endswith(".gz") and not os.path.exists("%s.tbi" % in_file):
        gatk.run("IndexFeatureFile", options)

    elif (in_file.endswith(".vcf") or in_file.endswith(".bed")) and (not os.path.exists("%s.idx" % in_file)):
        gatk.run("IndexFeatureFile", options)

    return


def gatk4_CombineGVCFs(gatk, in_gvcfs, intervals, outgvcffile):
    """CombineGVCFs"""
    options = [("-R", gatk.reference), ("-O", outgvcffile)]
    options += [("-V", gvcf) for gvcf in in_gvcfs]

    if intervals:
        options += [("-L", reg) for reg in intervals]

    gatk.run("CombineGVCFs", options)


def gatk4_MergeVCFs(gatk, in_vcfs, outvcffile):
    """ Merge VCFs or g.vcfs into a single vcf/g.vcf file.
    All in_vcfs should be the same samples.
    """
    options = [("-I", in_vcf) for in_vcf in in_vcfs]
    options += [("-O", outvcffile)]
    gatk.run("MergeVcfs", options)


def gatk4_GatherVCFs(gatk):
    pass


def gatk4_GenotypeGVCFs(gatk, in_gvcfs, intervals, outvcffile):

    if (not outvcffile.endswith(".vcf.gz")) and (not outvcffile.endswith(".vcf")):
        sys.stderr.write("Error: %s is not end with .vcf.gz or .vcf in "
                         "gatk4_GenotypeGVCFs.\n" % outvcffile)
        sys.exit(1)

    dirname = os.path.dirname(outvcffile)
    basename= os.path.basename(outvcffile).split(".vcf")[0]

    combine_gvcf_file = in_gvcfs[0]
    if len(in_gvcfs) > 1:
        # reset combine_gvcf_file
        combine_gvcf_file = os.path.join(dirname, "%s.combine.g.vcf.gz" % basename)
        gatk4_CombineGVCFs(gatk, in_gvcfs, intervals, combine_gvcf_file)

    # For GenotypeGVCFs
    options = [("-R", gatk.reference),
               ("-V", combine_gvcf_file),
               ("-O", outvcffile)]

    if intervals:
        options += [("-L", reg) for reg in intervals]

    gatk.run("GenotypeGVCFs", options)
    return outvcffile


def gatk4_VariantRecalibrator(gatk, in_vcfs, outvcffile, snp_max_gaussian=6, indel_max_gaussian=4,
                              snp_ts_filter_level=99.0, indel_ts_filter_level=95.0):

    if (not outvcffile.endswith(".vcf.gz")) and (not outvcffile.endswith(".vcf")):
        sys.stderr.write("Error: %s is not end with .vcf.gz or .vcf in "
                         "gatk4_VariantRecalibrator.\n" % outvcffile)
        sys.exit(1)

    dirname = os.path.dirname(outvcffile)
    basename= os.path.basename(outvcffile).split(".vcf")[0]

    snp_recal_file = os.path.join(dirname, "%s.snps.recal" % basename)
    snp_vqsr_file = os.path.join(dirname, "%s.snps.VQSR.vcf.gz" % basename)
    snp_tranches_file = os.path.join(dirname, "%s.snps.tranches" % basename)
    snp_rscriptfile = os.path.join(dirname, "%s.snps.plots.R" % basename)

    indel_recal_file = os.path.join(dirname, "%s.indels.recal" % basename)
    indel_tranches_file = os.path.join(dirname, "%s.indels.tranches" % basename)
    indel_rscriptfile = os.path.join(dirname, "%s.indels.plots.R" % basename)

    tranche = "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0"
    annotation = "-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum"

    # SNP
    snp_vqsr_opt = [("-R", gatk.reference), ("-O", snp_recal_file), ("--max-gaussians", snp_max_gaussian),
                    ("-mode", "SNP")]
    snp_vqsr_opt += [("-V", vcf) for vcf in in_vcfs]
    snp_vqsr_opt += [("-resource", "%s,%s:%s"%(n, r, f)) for n, r, f in gatk.resources["SNP"]]
    snp_vqsr_opt += [("", annotation), ("", tranche), ("--rscript-file", snp_rscriptfile),
                     ("--tranches-file", snp_tranches_file)]
    gatk.run("VariantRecalibrator", snp_vqsr_opt)

    snp_applyvqsr_opt = [("-R", gatk.reference), ("-O", snp_vqsr_file), ("-mode", "SNP")]
    snp_applyvqsr_opt += [("-V", vcf) for vcf in in_vcfs]
    snp_applyvqsr_opt += [("-ts-filter-level", snp_ts_filter_level),
                          ("--tranches-file", snp_tranches_file),
                          ("--recal-file", snp_recal_file)]
    gatk.run("ApplyVQSR", snp_applyvqsr_opt)

    # Indel
    indel_vqsr_opt = [("-R", gatk.reference), ("-V", snp_vqsr_file), ("--max-gaussians", indel_max_gaussian),
                      ("-O", indel_recal_file), ("-mode", "INDEL")]
    indel_vqsr_opt += [("-resource", "%s,%s:%s"%(n, r, f)) for n, r, f in gatk.resources["INDEL"]]
    indel_vqsr_opt += [("", annotation), ("", tranche), ("--rscript-file", indel_rscriptfile),
                       ("--tranches-file", indel_tranches_file)]
    gatk.run("VariantRecalibrator", indel_vqsr_opt)

    indel_applyvqsr_opt = [("-R", gatk.reference), ("-V", snp_vqsr_file),
                           ("-O", outvcffile), ("-mode", "INDEL")]
    indel_applyvqsr_opt += [("-ts-filter-level", indel_ts_filter_level),
                            ("--tranches-file", indel_tranches_file),
                            ("--recal-file", indel_recal_file)]
    gatk.run("ApplyVQSR", indel_applyvqsr_opt)

    return outvcffile, snp_recal_file, snp_tranches_file, snp_rscriptfile, \
           indel_recal_file, indel_tranches_file, indel_rscriptfile



