"""Run functions by provided a name with arguments.
"""
import sys
import os
import stat
import gzip

from ilus import utils
from ilus.modules.ngsaligner import bwa
from ilus.modules.bam.merge import mergebamfiles
from ilus.tools.sambamba import SambambaRunner, sambamba_mark_duplicates
from ilus.tools.gatk import GATKRunner, gatk4_mark_duplicates, gatk4_baserecalibrator, \
    gatk4_haplotypecaller, gatk4_GenotypeGVCFs, gatk4_VariantRecalibrator,gatk4_MergeVCFs


def bwamem(kwargs, out_folder_name, aione):
    """Run bwamem aligment for fastq to BAM"""

    # kwargs["args"].outdir is already an absolute path
    output_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "shell")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    if not utils.file_exists(kwargs["args"].fastqlist):
        sys.stderr.write("[ERROR] %s is not a file.\n" % kwargs["args"].fastqlist)
        sys.exit(1)

    sample_bamfiles_by_lane = {}  # {sample_id: [bwa1, bwa2, ...]}
    samples = []
    with gzip.open(kwargs["args"].fastqlist) if kwargs["args"].fastqlist.endswith(".gz") else \
            open(kwargs["args"].fastqlist) as I:

        # SAMPLE_ID RGID  FASTQ1  FASTQ2  LANE  LIBRARY PLATFORM        CENTER
        for line in I:
            if line[0] == "#":  # ignore header
                continue

            sample_id, rgID, fq1, fq2, lane = line.strip().split()[:5]
            sample_outdir = os.path.join(output_dirtory, sample_id)

            if sample_id not in sample_bamfiles_by_lane:
                utils.safe_makedir(sample_outdir)
                sample_bamfiles_by_lane[sample_id] = []

                # record the samples' id and keep the order as the same as input.
                samples.append([sample_id, sample_outdir])

            out_prefix = os.path.join(sample_outdir, sample_id + "_" + lane)
            lane_bam_file, cmd = bwa.bwa_mem(aione["config"], aione["config"]["resources"]["reference"],
                                             out_prefix, rgID, fq1, fq2)
            sample_bamfiles_by_lane[sample_id].append([lane_bam_file, cmd])

    # bwa_shell_files_dict = {}
    bwa_shell_files_list = []
    aione["sample_final_sorted_bam"] = []
    for sample, sample_outdir in samples:
        sample_final_bamfile = os.path.join(sample_outdir, sample+".sorted.bam")
        aione["sample_final_sorted_bam"].append([sample, sample_final_bamfile])

        lane_bam_files = []
        if len(sample_bamfiles_by_lane[sample]) == 1:
            lane_bam_files.append(sample_bamfiles_by_lane[sample][0])
            cmd = [sample_bamfiles_by_lane[sample][1]]

            if sample_final_bamfile != sample_bamfiles_by_lane[sample][0]:
                # single lane does not need to merge bamfiles
                cmd.append("mv -f %s %s" % (sample_bamfiles_by_lane[sample][0], sample_final_bamfile))

        else:
            samtools = aione["config"]["samtools"]["samtools"]
            samtools_merge_options = " ".join([str(x) for x in aione["config"]["samtools"].get("merge_options", [])])

            lane_bam_files = " ".join([f for f, _ in sample_bamfiles_by_lane[sample]])
            cmd = [c for _, c in sample_bamfiles_by_lane[sample]]

            # merge lane_bam_files into one and rm lane_bam_files
            cmd.append("{samtools} merge {samtools_merge_options} {sample_final_bamfile} "
                       "{lane_bam_files} && rm -rf {lane_bam_files}".format(**locals()))

        echo_mark_done = "echo \"[bwa] %s done\"" % sample
        cmd.append(echo_mark_done)

        sample_bwa_shell_file = os.path.join(shell_dirtory, sample+".bwa.sh")
        with open(sample_bwa_shell_file, "w") as OUT:
            OUT.write("#!/bin/bash\n")
            OUT.write("%s\n" % " && ".join(cmd))

        os.chmod(sample_bwa_shell_file, stat.S_IRWXU)  # 0700
        bwa_shell_files_list.append([sample, sample_bwa_shell_file])

    return bwa_shell_files_list  # [[sample, bwa_shell_file], ...]


def mergebam(kwargs, aione):
    return mergebamfiles(kwargs, aione)


def markduplicates(kwargs, aione):

    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output BAM files when running "
                         "markduplicates.\n")
        sys.exit(1)

    sambamba = SambambaRunner(aione["config"])
    aione["markdups_bam"] = sambamba_mark_duplicates(sambamba,
                                                     kwargs["args"].inbam,
                                                     kwargs["args"].outfile,
                                                     kwargs["args"].removedups)
    return aione


def gatk_markduplicates(kwargs, aione):

    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output BAM files when running "
                         "markduplicates.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"], aione["config"]["resources"]["gatk_bundle"])
    aione["markdups_bam"] = gatk4_mark_duplicates(gatk,
                                                  kwargs["args"].inbam,
                                                  kwargs["args"].outfile,
                                                  kwargs["args"].removedups)
    return aione


def baserecalibrator(kwargs, aione):
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output BAM files when running "
                         "BQSR.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"], aione["config"]["resources"]["gatk_bundle"])
    aione["bqsr_bam"] = gatk4_baserecalibrator(gatk,
                                               kwargs["args"].inbam,
                                               kwargs["args"].outfile,
                                               aione["config"]["resources"])
    return aione


def haplotypecaller(kwargs, aione):
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output VCF(or g.vcf) files when running "
                         "HaplotypeCaller.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"])
    name = "gvcf" if "g.vcf" in kwargs["args"].outfile else "genotype_vcf"

    aione[name] = gatk4_haplotypecaller(gatk,
                                        kwargs["args"].inbam,
                                        kwargs["args"].interval,
                                        kwargs["args"].outfile,
                                        kwargs["args"].emit_gvcf)
    return aione


def GenotypeGVCFs(kwargs, aione):
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output VCF files when running "
                         "GenotypeGVCFs.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"])
    aione["genotype_vcf"] = gatk4_GenotypeGVCFs(gatk,
                                                kwargs["args"].variants,
                                                kwargs["args"].interval,
                                                kwargs["args"].outfile)
    return aione


def variantrecalibrator(kwargs, aione):
    """Run VQSR"""
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output VCF files when running "
                         "GenotypeGVCFs.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"], aione["config"]["resources"]["gatk_bundle"])
    aione["qc_vcf"] = gatk4_VariantRecalibrator(gatk,
                                                kwargs["args"].variants,
                                                kwargs["args"].outfile,
                                                kwargs["args"].snp_max_gaussion,
                                                kwargs["args"].indel_max_gaussion,
                                                kwargs["args"].snp_ts_filter_level,
                                                kwargs["args"].indel_ts_filter_level)
    return aione

def MergeVCFs(kwargs,aione):
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output VCF files when running MergeVCFs.\n")
        sys.exit(1)

    gatk = GATKRunner(aione["config"])
    aione["MergeVCFs"] = gatk4_MergeVCFs(gatk,
                                         kwargs["args"].variants,
                                         kwargs["args"].outfile)
    return aione