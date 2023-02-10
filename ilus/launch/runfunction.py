"""Run functions by provided a name with arguments.

Author:  Shujia Huang
Date: 2020-04-19
"""
import sys
import os
import stat
import gzip

from ilus import utils
from ilus.modules.ngsaligner import bwa
from ilus.modules.variants import gatk
from ilus.modules.vcf import bcftools
from ilus.modules.summary import bam

IS_RM_SUBBAM = True


def _md(outdir, out_folder_name, is_dry_run=False):
    """Create folders for output and shell scripts in the same format.
    """
    output_dirtory = os.path.join(outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(outdir, out_folder_name, "shell")
    if not is_dry_run:
        utils.safe_makedir(output_dirtory)
        utils.safe_makedir(shell_dirtory)

    return output_dirtory, shell_dirtory


def _create_cmd_file(out_shell_file, cmd=None):
    if cmd is None:
        return

    with open(out_shell_file, "w") as OUT:
        OUT.write("#!/bin/bash\n")
        OUT.write("%s\n" % " && ".join(cmd))

    os.chmod(out_shell_file, stat.S_IRWXU)  # 0700
    return


def bwamem(kwargs, out_folder_name, aione, is_dry_run=False):
    """Create bwamem aligment shell scripts for fastq to BAM.
    """
    output_dirtory = os.path.join(kwargs.outdir, out_folder_name, "output")
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "bwa")
    if not is_dry_run:
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
                sample_bamfiles_by_lane[sample_id] = []
                if not is_dry_run:
                    utils.safe_makedir(sample_outdir)

                # record the samples' id and keep the order as the same as input.
                samples.append([sample_id, sample_outdir])

            out_prefix = os.path.join(sample_outdir, sample_id + "_" + lane)
            lane_bam_file, cmd = bwa.bwa_mem(aione["config"], out_prefix, rgID, fq1, fq2)
            sample_bamfiles_by_lane[sample_id].append([lane_bam_file, cmd])

    bwa_shell_files_list = []
    aione["sample_final_sorted_bam"] = []
    for sample, sample_outdir in samples:
        sample_final_bamfile = os.path.join(sample_outdir, sample + ".sorted.bam")
        aione["sample_final_sorted_bam"].append([sample, sample_final_bamfile])
        sample_shell_fname = os.path.join(shell_dirtory, sample + ".bwa.sh")
        bwa_shell_files_list.append([sample, sample_shell_fname])

        # Create the commandline script.
        cmd = []
        if len(sample_bamfiles_by_lane[sample]) == 1:
            lane_bam_file, cmd = sample_bamfiles_by_lane[sample][0][0], [sample_bamfiles_by_lane[sample][0][1]]
            if sample_final_bamfile != lane_bam_file:
                # single lane does not need to merge bamfiles
                cmd.append("mv -f %s %s" % (lane_bam_file, sample_final_bamfile))

        else:
            samtools = aione["config"]["samtools"]["samtools"]
            samtools_merge_options = " ".join(
                [str(x) for x in aione["config"]["samtools"].get("merge_options", [])])

            lane_bam_files = " ".join([f for f, _ in sample_bamfiles_by_lane[sample]])
            cmd = [c for _, c in sample_bamfiles_by_lane[sample]]

            # merge lane_bam_files into one and rm lane_bam_files
            cmd.append("{samtools} merge {samtools_merge_options} {sample_final_bamfile} "
                       "{lane_bam_files} && rm -rf {lane_bam_files}".format(**locals()))

        echo_mark_done = "echo \"[bwa] %s done\"" % sample
        cmd.append(echo_mark_done)

        if not is_dry_run and (not os.path.exists(sample_shell_fname) or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return bwa_shell_files_list  # [[sample, bwa_shell_file], ...]


def gatk_markduplicates(kwargs, out_folder_name, aione, is_dry_run=False):
    """Create shell scripts for Markduplicates by GATK4.
    """
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "markdup")
    if not is_dry_run:
        utils.safe_makedir(shell_dirtory)

    aione["sample_final_markdup_bam"] = []
    markdup_shell_files_list = []
    for sample, sample_sorted_bam in aione["sample_final_sorted_bam"]:
        dirname, f_name = os.path.split(sample_sorted_bam)

        # Setting Output path of markduplicate BAM file as the same as ``sample_sorted_bam``
        out_markdup_bam_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".markdup.bam")
        out_markdup_metrics_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".metrics.txt")
        sample_shell_fname = os.path.join(shell_dirtory, sample + ".markdup.sh")

        markdup_shell_files_list.append([sample, sample_shell_fname])
        aione["sample_final_markdup_bam"].append([sample, out_markdup_bam_fname])

        cmd = [gatk.markduplicates(aione["config"], sample_sorted_bam, out_markdup_bam_fname,
                                   out_markdup_metrics_fname)]

        if IS_RM_SUBBAM:
            cmd.append("rm -rf %s" % sample_sorted_bam)  # save disk space

        echo_mark_done = "echo \"[MarkDuplicates] %s done\"" % sample
        cmd.append(echo_mark_done)

        if not is_dry_run and (not os.path.exists(sample_shell_fname) or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return markdup_shell_files_list  # [[sample, sample_shell_fname], ...]


def gatk_baserecalibrator(kwargs, out_folder_name, aione, is_calculate_summary=True, is_dry_run=False):
    shell_dirtory = os.path.join(kwargs.outdir, out_folder_name, "shell", "bqsr")
    if not is_dry_run:
        utils.safe_makedir(shell_dirtory)

    aione["sample_final_bqsr_bam"] = []
    bqsr_shell_files_list = []

    is_calculate_contamination = True if "verifyBamID2" in aione["config"] else False
    for sample, sample_markdup_bam in aione["sample_final_markdup_bam"]:
        dirname, f_name = os.path.split(sample_markdup_bam)

        out_bqsr_bam_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.bam")
        out_bqsr_bai_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.bai")
        out_bqsr_recal_table = os.path.join(dirname, os.path.splitext(f_name)[0] + ".recal.table")

        out_alignment_summary_metric = os.path.join(
            dirname, os.path.splitext(f_name)[0] + ".AlignmentSummaryMetrics.txt")

        out_bamstats_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.stats")
        genome_cvg_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.depth.bed.gz")
        # when convert to CRAM format
        out_cram_fname = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.cram")

        aione["sample_final_bqsr_bam"].append([sample, out_bqsr_bam_fname if not kwargs.cram else out_cram_fname])
        sample_shell_fname = os.path.join(shell_dirtory, sample + ".bqsr.sh")
        bqsr_shell_files_list.append([sample, sample_shell_fname])

        cmd = [gatk.baserecalibrator(aione["config"],
                                     sample_markdup_bam,
                                     out_bqsr_bam_fname,
                                     out_bqsr_recal_table)]
        if IS_RM_SUBBAM:
            cmd.append("rm -rf %s" % sample_markdup_bam)

        if is_calculate_summary:
            cmd.append(gatk.collect_alignment_summary_metrics(
                aione["config"], out_bqsr_bam_fname, out_alignment_summary_metric)
            )
            cmd.append(bam.stats(aione["config"], out_bqsr_bam_fname, out_bamstats_fname))
            cmd.append(bam.genomecoverage(aione["config"], out_bqsr_bam_fname, genome_cvg_fname))

        if is_calculate_contamination:
            out_verifybamid_stat_prefix = os.path.join(dirname, os.path.splitext(f_name)[0] + ".BQSR.verifyBamID2")
            cmd.append(bam.verifyBamID2(aione["config"], out_bqsr_bam_fname, out_verifybamid_stat_prefix))

        if kwargs.cram:
            cmd.append(bwa.bam_to_cram(aione["config"], out_bqsr_bam_fname, out_cram_fname))
            cmd.append("rm -rf %s" % out_bqsr_bam_fname)
            cmd.append("rm -rf %s" % out_bqsr_bai_fname)

        echo_mark_done = "echo \"[BQSR] %s done\"" % sample
        cmd.append(echo_mark_done)

        if not is_dry_run and (not os.path.exists(sample_shell_fname) or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return bqsr_shell_files_list


def gatk_haplotypecaller_gvcf(kwargs, out_folder_name, aione, is_dry_run=False):
    """Create gvcf shell.
    """

    def _create_sub_shell(sample, sample_shell_dir, sample_output_dir, raw_interval=None, is_dry_run=False):
        interval_ = raw_interval if raw_interval else "all"
        # in case the raw_interval is a full path file.
        interval_, _ = os.path.splitext(os.path.split(interval_)[-1])
        sample_shell_fname_ = os.path.join(sample_shell_dir, sample + ".%s.gvcf.sh" % interval_)
        out_gvcf_fname_ = os.path.join(sample_output_dir, sample + ".%s.g.vcf.gz" % interval_)

        if raw_interval:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname_,
                                             interval=raw_interval)]
        else:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname_)]
        echo_mark_done = "echo \"[GVCF] %s %s done\"" % (sample, interval_)
        cmd.append(echo_mark_done)

        if (not is_dry_run) and (not os.path.exists(sample_shell_fname_) or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname_, cmd)

        return sample_shell_fname_, out_gvcf_fname_

    output_dirtory, shell_dirtory = _md(kwargs.outdir, out_folder_name, is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    gvcf_shell_files_list = []
    aione["gvcf"] = {}

    if "interval" not in aione["config"]["gatk"]:
        aione["config"]["gatk"]["interval"] = ["all"]

    aione["intervals"] = []  # chromosome id
    sample_map = {}  # Record sample_map for genomicsDBImport
    for sample, sample_bqsr_bam in aione["sample_final_bqsr_bam"]:

        sample_output_dir = os.path.join(output_dirtory, sample)
        sample_shell_dir = os.path.join(shell_dirtory, sample)
        if not is_dry_run:
            utils.safe_makedir(sample_output_dir)
            utils.safe_makedir(sample_shell_dir)

        for interval in aione["config"]["gatk"]["interval"]:
            if interval == "all":
                # All the whole genome once, take a lot of time, not suggested
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, is_dry_run=is_dry_run)
            else:
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, raw_interval=interval, is_dry_run=is_dry_run)

            # ``interval`` and ``aione["config"]["gatk"]["interval"]`` may be different.
            # The raw interval could be a file path.
            interval, _ = os.path.splitext(os.path.split(interval)[-1])
            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            gvcf_shell_files_list.append([sample + ".%s" % interval, sample_shell_fname])
            aione["gvcf"][interval].append(out_gvcf_fname)

            if is_use_gDBI and (interval not in sample_map):
                sample_map[interval] = []
            elif is_use_gDBI:
                sample_map[interval].append("%s\t%s" % (sample, out_gvcf_fname))

    if is_use_gDBI:
        aione["sample_map"] = {}
        for interval, value in sample_map.items():
            # create sample_map file for next process
            out_sample_map_fname = os.path.join(output_dirtory, "%s.gvcf.sample_map" % interval)
            aione["sample_map"][interval] = out_sample_map_fname

            if (not is_dry_run) and ((not os.path.exists(out_sample_map_fname)) or kwargs.overwrite):
                with open(out_sample_map_fname, "w") as OUT:
                    OUT.write("\n".join(value))

    return gvcf_shell_files_list


def _yield_gatk_combineGVCFs(project_name, variant_calling_intervals, output_dirtory, shell_dirtory,
                             is_use_gDBI, aione=None):
    for interval in variant_calling_intervals:
        interval_id = interval[0] if type(interval) is list else interval  # get chromosome id
        if interval_id in aione["gvcf"]:
            # load gvcf '.sample_map' file if using genomicsDBImport module to combine all the gvcf
            sample_gvcf_list = [aione["sample_map"][interval_id]] if is_use_gDBI else aione["gvcf"][interval_id]
        else:
            continue

        if not sample_gvcf_list:
            sys.stderr.write("[Error] Interval parameters in configuration file may be different "
                             "from that input gvcf file list when joint-calling by genotypeGVCFs.")
            sys.exit(1)

        interval_n = "_".join(interval) if type(interval) is list else interval.replace(":", "_")  # chr:s -> chr_s
        sub_shell_fname = os.path.join(shell_dirtory, "%s.%s.combineGVCFs.sh" % (project_name, interval_n))
        combineGVCF_fname = os.path.join(output_dirtory, "%s.%s.genomics_db" % (project_name, interval_n)) \
            if is_use_gDBI else os.path.join(output_dirtory, "%s.%s.g.vcf.gz" % (project_name, interval_n))

        # Create commandline for combining GVCFs
        calling_interval = "%s:%s-%s" % (interval[0], interval[1], interval[2]) \
            if type(interval) is list else interval

        combineGVCFs_cmd = gatk.combineGVCFs(aione["config"],
                                             sample_gvcf_list,
                                             combineGVCF_fname,
                                             interval=calling_interval)

        yield combineGVCFs_cmd, combineGVCF_fname, sub_shell_fname, interval_n, calling_interval


def gatk_combineGVCFs(kwargs, out_folder_name, aione, is_dry_run=False):
    """Combine GVCFs by GATK genomicsDBImport or CombineGVCFs module.
    """
    output_dirtory, shell_dirtory = _md(kwargs.outdir, out_folder_name, is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    combineGVCFs_shell_files_list = []
    aione["combineGVCFs"] = {}

    variant_calling_intervals = aione["config"]["gatk"]["variant_calling_interval"]
    for (combineGVCFs_cmd,
         combineGVCF_fname,
         sub_shell_fname,
         interval_n,  # interval marker
         calling_interval) in _yield_gatk_combineGVCFs(kwargs.project_name,
                                                       variant_calling_intervals,
                                                       output_dirtory,
                                                       shell_dirtory,
                                                       is_use_gDBI,
                                                       aione):

        combineGVCFs_shell_files_list.append([kwargs.project_name + "." + interval_n, sub_shell_fname])
        aione["combineGVCFs"][interval_n] = combineGVCF_fname

        echo_mark_done = "echo \"[CombineGVCFs] %s done\"" % calling_interval
        cmd = [combineGVCFs_cmd, echo_mark_done]

        if (not is_dry_run) and (not os.path.exists(sub_shell_fname) or kwargs.overwrite):
            _create_cmd_file(sub_shell_fname, cmd)

    return combineGVCFs_shell_files_list


# def gatk_genotypeGVCFs(kwargs, out_folder_name, aione, is_dry_run=False):
#     """Create shell for genotypeGVCFs.
#     """
#     output_dirtory, shell_dirtory = _md(kwargs.outdir, out_folder_name, is_dry_run=is_dry_run)
#     is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
#         if "use_genomicsDBImport" in aione["config"]["gatk"] else False
#
#     genotype_vcf_shell_files_list = []
#     aione["genotype_vcf_list"] = []
#
#     # Create commandline for combineGVCFs and genotypeGVCFs process for specific calling intervals
#     variant_calling_intervals = aione["config"]["gatk"]["variant_calling_interval"]
#     for (combineGVCFs_cmd,
#          combineGVCF_fname,
#          _,  # no need the sub_shell_fname in _yield_gatk_CombineGVCFs
#          interval_n,  # interval marker
#          calling_interval) in _yield_gatk_combineGVCFs(kwargs.project_name,
#                                                        variant_calling_intervals,
#                                                        output_dirtory,
#                                                        shell_dirtory,
#                                                        is_use_gDBI,
#                                                        aione):
#
#         sub_shell_fname = os.path.join(shell_dirtory, "%s.%s.genotype.sh" % (kwargs.project_name, interval_n))
#         genotype_vcf_fname = os.path.join(output_dirtory, "%s.%s.vcf.gz" % (kwargs.project_name, interval_n))
#
#         genotype_vcf_shell_files_list.append([kwargs.project_name + "." + interval_n, sub_shell_fname])
#         aione["genotype_vcf_list"].append(genotype_vcf_fname)
#
#         # Generate the genotype joint-calling process according to the input file
#         # ``combineGVCF_fname`` which create in  ``_yield_gatk_CombineGVCFs``
#         genotype_cmd = gatk.genotypeGVCFs(aione["config"],
#                                           combineGVCF_fname,
#                                           genotype_vcf_fname,
#                                           interval=calling_interval)
#
#         # delete the genomicsdb workspace (or the combine GVCF file).
#         delete_cmd = "rm -rf %s" % combineGVCF_fname \
#             if is_use_gDBI else "rm -rf %s %s.tbi" % (combineGVCF_fname, combineGVCF_fname)
#
#         echo_mark_done = "echo \"[Genotype] %s done\"" % calling_interval
#
#         # combine all the commands for the final genotype calling.
#         cmd = [combineGVCFs_cmd, genotype_cmd, delete_cmd, echo_mark_done]
#
#         if (not is_dry_run) and (not os.path.exists(sub_shell_fname) or kwargs.overwrite):
#             # Create shell file
#             _create_cmd_file(sub_shell_fname, cmd)
#
#     return genotype_vcf_shell_files_list


def gatk_genotypeGVCFs(kwargs, out_folder_name, aione, is_dry_run=False):
    """Create shell for genotypeGVCFs.
        """
    output_dirtory, shell_dirtory = _md(kwargs.outdir, out_folder_name, is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    genotype_vcf_shell_files_list = []
    aione["genotype_vcf_list"] = []

    variant_calling_intervals = aione["config"]["gatk"]["variant_calling_interval"]
    for interval in variant_calling_intervals:
        interval_n = "_".join(interval) if type(interval) is list else interval.replace(":", "_")  # chr:s -> chr_s
        sub_shell_fname = os.path.join(shell_dirtory, "%s.%s.genotype.sh" % (kwargs.project_name, interval_n))
        genotype_vcf_fname = os.path.join(output_dirtory, "%s.%s.vcf.gz" % (kwargs.project_name, interval_n))
        combineGVCF_fname = aione["combineGVCFs"][interval_n]

        genotype_vcf_shell_files_list.append([kwargs.project_name + "." + interval_n, sub_shell_fname])
        aione["genotype_vcf_list"].append(genotype_vcf_fname)

        calling_interval = "%s:%s-%s" % (interval[0], interval[1], interval[2]) \
            if type(interval) is list else interval

        # Create commandline for genotype joint-calling process
        genotype_cmd = gatk.genotypeGVCFs(aione["config"],
                                          combineGVCF_fname,
                                          genotype_vcf_fname,
                                          interval=calling_interval)

        # delete the genomicsdb workspace (or the combine GVCF file).
        delete_cmd = "rm -rf %s" % combineGVCF_fname \
            if is_use_gDBI else "rm -rf %s %s.tbi" % (combineGVCF_fname, combineGVCF_fname)
        echo_mark_done = "echo \"[Genotype] %s done\"" % calling_interval
        cmd = [genotype_cmd, delete_cmd, echo_mark_done]

        if (not is_dry_run) and (not os.path.exists(sub_shell_fname) or kwargs.overwrite):
            _create_cmd_file(sub_shell_fname, cmd)

    return genotype_vcf_shell_files_list


def gatk_variantrecalibrator(kwargs, out_folder_name, aione, is_dry_run=False):
    """Create shell scripts for VQSR.
    """
    output_dirtory, shell_dirtory = _md(kwargs.outdir, out_folder_name, is_dry_run=is_dry_run)
    genotype_vqsr_fname = os.path.join(output_dirtory, "%s.VQSR.vcf.gz" % kwargs.project_name)
    combine_vcf_fname = os.path.join(output_dirtory, "%s.raw.vcf.gz" % kwargs.project_name)
    shell_fname = os.path.join(shell_dirtory, "%s.VQSR.sh" % kwargs.project_name)

    cmd = []
    if len(aione["genotype_vcf_list"]) > 1:
        # concat-vcf
        concat_vcf_cmd = bcftools.concat(aione["config"], aione["genotype_vcf_list"], combine_vcf_fname)
        cmd.append(concat_vcf_cmd)
    else:
        combine_vcf_fname = aione["genotype_vcf_list"][0]

    # VQSR process
    cmd.append(gatk.variantrecalibrator(aione["config"], combine_vcf_fname, genotype_vqsr_fname))
    cmd.append("echo \"[VQSR] %s done\"" % genotype_vqsr_fname)

    if not is_dry_run and (not os.path.exists(shell_fname) or kwargs.overwrite):
        _create_cmd_file(shell_fname, cmd)

    # Only one VQSR result
    return [["%s.VQSR" % kwargs.project_name, shell_fname]]
