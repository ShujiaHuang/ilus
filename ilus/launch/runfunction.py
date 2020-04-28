"""Run functions by provided a name with arguments.
"""
import sys
import os
import stat
import gzip

from ilus import utils
from ilus.modules.ngsaligner import bwa
from ilus.modules.variants import gatk
from ilus.modules.vcf import bedtools

IS_RM_SUBBAM = True


def _create_cmd_file(out_shell_file, cmd):
    with open(out_shell_file, "w") as OUT:
        OUT.write("#!/bin/bash\n")
        OUT.write("%s\n" % " && ".join(cmd))

    os.chmod(out_shell_file, stat.S_IRWXU)  # 0700


def bwamem(kwargs, out_folder_name, aione):
    """Run bwamem aligment for fastq to BAM"""
    output_dirtory = os.path.join(kwargs.outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "bwa")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    if not utils.file_exists(kwargs.fastqlist):
        sys.stderr.write("[ERROR] %s is not a file.\n" % kwargs.fastqlist)
        sys.exit(1)

    sample_bamfiles_by_lane = {}  # {sample_id: [bwa1, bwa2, ...]}
    samples = []
    with gzip.open(kwargs.fastqlist) if kwargs.fastqlist.endswith(".gz") else \
            open(kwargs.fastqlist) as I:

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

        sample_shell_fname = os.path.join(shell_dirtory, sample + ".bwa.sh")
        if not os.path.exists(sample_shell_fname) or kwargs.overwrite:
            _create_cmd_file(sample_shell_fname, cmd)

        bwa_shell_files_list.append([sample, sample_shell_fname])

    return bwa_shell_files_list  # [[sample, bwa_shell_file], ...]


def gatk_markduplicates(kwargs, out_folder_name, aione):
    """Markduplicates by GATK4
    """
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "markdup")
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

        sample_shell_fname = os.path.join(shell_dirtory, sample + ".markdup.sh")
        if not os.path.exists(sample_shell_fname) or kwargs.overwrite:
            _create_cmd_file(sample_shell_fname, cmd)

        markdup_shell_files_list.append([sample, sample_shell_fname])
        aione["sample_final_markdup_bam"].append([sample, out_markdup_bam_fname])

    return markdup_shell_files_list  # [[sample, sample_shell_fname], ...]


def gatk_baserecalibrator(kwargs, out_folder_name, aione):
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "bqsr")
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
        sample_shell_fname = os.path.join(shell_dirtory, sample + ".bqsr.sh")
        if not os.path.exists(sample_shell_fname) or kwargs.overwrite:
            _create_cmd_file(sample_shell_fname, cmd)

        bqsr_shell_files_list.append([sample, sample_shell_fname])
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
        if not os.path.exists(sample_shell_fname) or kwargs.overwrite:
            _create_cmd_file(sample_shell_fname, cmd)

        return sample_shell_fname, out_gvcf_fname

    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell")
    output_dirtory = os.path.join(kwargs.outdir, out_folder_name, "output")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    gvcf_shell_files_list = []
    aione["gvcf"] = {}

    if "interval" not in aione["config"]["gatk"]:
        aione["config"]["gatk"]["interval"] = ["all"]

    aione["intervals"] = []
    for sample, sample_bqsr_bam in aione["sample_final_bqsr_bam"]:
        sample_shell_dir = os.path.join(shell_dirtory, sample)
        sample_output_dir = os.path.join(output_dirtory, sample)
        utils.safe_makedir(sample_shell_dir)
        utils.safe_makedir(sample_output_dir)

        for interval in aione["config"]["gatk"]["interval"]:

            if interval == "all":
                # The whole genome
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(sample, sample_shell_dir, sample_output_dir)
            else:
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, raw_interval=interval)

            # ``interval`` and ``aione["config"]["gatk"]["interval"]`` could be different.
            # The raw interval could be a file path.
            interval, _ = os.path.splitext(os.path.split(interval)[-1])
            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            gvcf_shell_files_list.append([sample + ".%s" % interval, sample_shell_fname])
            aione["gvcf"][interval].append(out_gvcf_fname)

    return gvcf_shell_files_list


def gatk_genotypeGVCFs(kwargs, out_folder_name, aione):
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell")
    output_dirtory = os.path.join(kwargs.outdir, out_folder_name, "output")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    genotype_vcf_shell_files_list = []
    aione["genotype"] = {}
    for interval in aione["intervals"]:
        genotype_vcf_fname = os.path.join(output_dirtory, "%s.%s.vcf.gz" % (kwargs.project_name, interval))
        sub_shell_fname = os.path.join(shell_dirtory, "%s.%s.genotype.sh" % (kwargs.project_name, interval))
        sample_gvcf_list = aione["gvcf"][interval]
        cmd = [gatk.genotypegvcfs(aione["config"], sample_gvcf_list, genotype_vcf_fname)]

        echo_mark_done = "echo \"[Genotype] %s done\"" % interval
        cmd.append(echo_mark_done)
        if not os.path.exists(sub_shell_fname) or kwargs.overwrite:
            _create_cmd_file(sub_shell_fname, cmd)

        genotype_vcf_shell_files_list.append([interval, sub_shell_fname])
        aione["genotype"][interval] = genotype_vcf_fname

    return genotype_vcf_shell_files_list


def gatk_variantrecalibrator(kwargs, out_folder_name, aione):
    """Run VQSR"""
    output_dirtory = os.path.join(kwargs.outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell")
    utils.safe_makedir(output_dirtory)
    utils.safe_makedir(shell_dirtory)

    genotype_vqsr_fname = os.path.join(output_dirtory, "%s.VQSR.vcf.gz" % kwargs.project_name)
    combine_vcf_fname = os.path.join(output_dirtory, "%s.raw.vcf.gz" % kwargs.project_name)
    shell_fname = os.path.join(shell_dirtory, "%s.VQSR.sh" % kwargs.project_name)

    cmd = []
    if len(aione["intervals"]) > 1:
        # concat-vcf
        concat_vcf_cmd = bedtools.concat(aione["config"],
                                         [aione["genotype"][interval] for interval in aione["intervals"]],
                                         combine_vcf_fname)
        cmd.append(concat_vcf_cmd)
    else:
        combine_vcf_fname = aione["genotype"][aione["interval"][0]]

    # VQSR
    echo_mark_done = "echo \"[VQSR] %s done\"" % genotype_vqsr_fname
    cmd.append(echo_mark_done)
    cmd.append(gatk.variantrecalibrator(aione["config"], combine_vcf_fname, genotype_vqsr_fname))

    if not os.path.exists(shell_fname) or kwargs.overwrite:
        _create_cmd_file(shell_fname, cmd)

    # Only one VQSR result
    return [["%s.VQSR" % kwargs.project_name, shell_fname]]
