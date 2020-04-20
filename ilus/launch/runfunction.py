"""Run functions by provided a name with arguments.
"""
import sys
import os
import stat
import gzip

from ilus import utils
from ilus.modules.ngsaligner import bwa
from ilus.modules.variants import gatk

from ilus.tools.sambamba import SambambaRunner, sambamba_mark_duplicates

from ilus.tools.gatk import GATKRunner, gatk4_mark_duplicates, gatk4_baserecalibrator, \
    gatk4_haplotypecaller, gatk4_GenotypeGVCFs, gatk4_VariantRecalibrator, gatk4_MergeVCFs

IS_RM_SUBBAM = True


def _create_cmd_file(out_shell_file, cmd):
    with open(out_shell_file, "w") as OUT:
        OUT.write("#!/bin/bash\n")
        OUT.write("%s\n" % " && ".join(cmd))

    os.chmod(out_shell_file, stat.S_IRWXU)  # 0700


def bwamem(kwargs, out_folder_name, aione):
    """Run bwamem aligment for fastq to BAM"""
    output_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "shell", "bwa")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    if not utils.file_exists(kwargs["args"].fastqlist):
        sys.stderr.write("[ERROR] %s is not a file.\n" % kwargs["args"].fastqlist)
        sys.exit(1)

    sample_bamfiles_by_lane = {}  # {sample_id: [bwa1, bwa2, ...]}
    samples = []
    with gzip.open(kwargs["args"].fastqlist) if kwargs["args"].fastqlist.endswith(".gz") else \
            open(kwargs["args"].fastqlist) as I:

        # SAMPLE_ID RGID  FASTQ1  FASTQ2  LANE  LIBRARY PLATFORM   CENTER
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
            lane_bam_file, cmd = bwa.bwa_mem(aione["config"], out_prefix, rgID, fq1, fq2)
            sample_bamfiles_by_lane[sample_id].append([lane_bam_file, cmd])

    # bwa_shell_files_dict = {}
    bwa_shell_files_list = []
    aione["sample_final_sorted_bam"] = []
    for sample, sample_outdir in samples:
        sample_final_bamfile = os.path.join(sample_outdir, sample + ".sorted.bam")
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

        sample_shell_file = os.path.join(shell_dirtory, sample + ".bwa.sh")
        _create_cmd_file(sample_shell_file, cmd)
        bwa_shell_files_list.append([sample, sample_shell_file])

    return bwa_shell_files_list  # [[sample, bwa_shell_file], ...]


def gatk_markduplicates(kwargs, out_folder_name, aione):
    """Markduplicates by GATK4
    """
    shell_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "shell", "markdup")
    utils.safe_makedir(shell_dirtory)

    aione["sample_final_markdup_bam"] = []
    markdup_shell_files_list = []
    for sample, sample_sorted_bam in aione["sample_final_sorted_bam"]:
        dirname, f_name = os.path.split(sample_sorted_bam)

        # Setting Output path of markduplicate BAM file as the same as ``sample_sorted_bam``
        out_markdup_bam_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".markdup.bam")
        out_markdup_metrics_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".metrics.txt")
        cmd = [
            gatk.markduplicates(aione["config"], sample_sorted_bam, out_markdup_bam_fname, out_markdup_metrics_fname)]

        if IS_RM_SUBBAM:
            cmd.append("rm -rf %s" % sample_sorted_bam)  # save disk space

        echo_mark_done = "echo \"[MarkDuplicates] %s done\"" % sample
        cmd.append(echo_mark_done)

        sample_shell_file = os.path.join(shell_dirtory, sample + ".markdup.sh")
        _create_cmd_file(sample_shell_file, cmd)
        markdup_shell_files_list.append([sample, sample_shell_file])
        aione["sample_final_markdup_bam"].append([sample, out_markdup_bam_fname])

    return markdup_shell_files_list  # [[sample, sample_shell_file], ...]


def gatk_baserecalibrator(kwargs, out_folder_name, aione):
    shell_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "shell", "bqsr")
    utils.safe_makedir(shell_dirtory)

    aione["sample_final_bqsr_bam"] = []
    bqsr_shell_files_list = []
    for sample, sample_markdup_bam in aione["sample_final_markdup_bam"]:
        dirname, f_name = os.path.split(sample_markdup_bam)

        out_bqsr_bam_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.bam")
        out_bqsr_recal_table = os.path.join(dirname, os.path.splitext(f_name)[0] + ".recal.table")

        cmd = [gatk.baserecalibrator(aione["config"], sample_markdup_bam, out_bqsr_bam_fname, out_bqsr_recal_table)]
        if IS_RM_SUBBAM:
            cmd.append("rm -rf %s" % sample_markdup_bam)

        echo_mark_done = "echo \"[BQSR] %s done\"" % sample
        cmd.append(echo_mark_done)
        sample_shell_file = os.path.join(shell_dirtory, sample + ".bqsr.sh")
        _create_cmd_file(sample_shell_file, cmd)
        bqsr_shell_files_list.append([sample, sample_shell_file])
        aione["sample_final_bqsr_bam"].append([sample, out_bqsr_bam_fname])

    return bqsr_shell_files_list


def gatk_haplotypecaller_gvcf(kwargs, out_folder_name, aione):
    """Create gvcf shell."""

    def _create_sub_shell(sample, sample_shell_dir, sample_output_dir, raw_interval=None):

        interval = raw_interval if raw_interval else "all"

        # in case the raw_interval is a full path file.
        interval, _ = os.path.splitext(os.path.split(interval)[-1])
        sample_shell_fname = os.path.join(sample_shell_dir, sample + ".%s.gvcf.sh" % interval)
        out_gvcf_fname = os.path.join(sample_output_dir, sample + ".%s.g.vcf.gz" % interval)

        if raw_interval:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname, interval=raw_interval)]
        else:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname)]

        echo_mark_done = "echo \"[GVCF] %s %s done\"" % (sample, interval)
        cmd.append(echo_mark_done)
        _create_cmd_file(sample_shell_fname, cmd)

        return sample_shell_fname, out_gvcf_fname

    shell_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "shell")
    output_dirtory = os.path.join(kwargs["args"].outdir, out_folder_name, "output")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    gvcf_shell_files_list = []
    aione["gvcf"] = {}
    for sample, sample_bqsr_bam in aione["sample_final_bqsr_bam"]:
        sample_shell_dir = os.path.join(shell_dirtory, sample)
        sample_output_dir = os.path.join(output_dirtory, sample)
        utils.safe_makedir(sample_shell_dir)
        utils.safe_makedir(sample_output_dir)

        if "interval" in aione["config"]["gatk"]:
            for interval in aione["config"]["gatk"]["interval"]:

                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, raw_interval=interval)

                interval, _ = os.path.splitext(os.path.split(interval)[-1])
                gvcf_shell_files_list.append([sample + ".%s" % interval, sample_shell_fname])
                if interval not in aione["gvcf"]:
                    aione["gvcf"][interval] = []

                aione["gvcf"][interval].append(out_gvcf_fname)

        else:

            # The whole genome
            interval = "all"
            sample_shell_fname, out_gvcf_fname = _create_sub_shell(sample, sample_shell_dir, sample_output_dir)
            if interval not in aione["gvcf"]:
                aione["gvcf"][interval] = []

            aione["gvcf"][interval].append(out_gvcf_fname)

    return gvcf_shell_files_list


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
