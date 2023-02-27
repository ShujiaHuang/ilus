"""Run functions by provided a name with arguments.

Author:  Shujia Huang
Date: 2020-04-19
"""
import sys
import stat
import gzip
from pathlib import Path
from typing import List, Tuple

from ilus.modules.utils import safe_makedir
from ilus.modules.ngsaligner import bwa
from ilus.modules.variants import gatk
from ilus.modules.vcf import bcftools
from ilus.modules.summary import bam

IS_RM_SUBBAM = True


def _md(outdir: Path, is_dry_run: bool = False) -> Tuple[Path, Path]:
    """Create folders for output and shell scripts in the same format.
    """
    output_dirtory = Path(outdir).joinpath("output")
    shell_dirtory = Path(outdir).joinpath("shell")
    if not is_dry_run:
        safe_makedir(output_dirtory)
        safe_makedir(shell_dirtory)

    # Return two dirctory in Path type
    return output_dirtory, shell_dirtory


def _create_cmd_file(out_shell_file: Path, cmd: List[str] = None):
    if cmd is None:
        return

    with open(str(out_shell_file), "w") as OUT:
        OUT.write("#!/bin/bash\n")
        OUT.write(f"{' && '.join(cmd)}\n")

    Path(out_shell_file).chmod(stat.S_IRWXU)  # 0700
    return


def bwamem(kwargs, out_folder_name: str, aione: dict = None,
           is_dry_run: bool = False):
    """Create bwamem aligment shell scripts for fastq to BAM.
    """
    output_dirtory = Path(kwargs.outdir).joinpath(out_folder_name, "output")
    shell_dirtory = Path(kwargs.outdir).joinpath(out_folder_name, "shell", "bwa")
    if not is_dry_run:
        safe_makedir(output_dirtory)
        safe_makedir(shell_dirtory)

    sample_bamfiles_by_lane = {}  # {sample_id: [bwa1, bwa2, ...]}
    samples = []
    with gzip.open(aione["fastqlist"]) if aione["fastqlist"].endswith(".gz") else \
            open(aione["fastqlist"]) as f:

        # SAMPLE_ID RGID  FASTQ1  FASTQ2  LANE  LIBRARY PLATFORM   CENTER
        for line in f:
            if line[0] == "#":  # ignore header
                continue

            sample_id, read_group_id, fq1, fq2, lane = line.strip().split()[:5]
            sample_outdir = output_dirtory.joinpath(sample_id)
            if sample_id not in sample_bamfiles_by_lane:
                sample_bamfiles_by_lane[sample_id] = []
                if not is_dry_run:
                    safe_makedir(sample_outdir)

                # record the samples' id and keep the order as the same as input.
                samples.append([sample_id, sample_outdir])

            out_prefix = sample_outdir.joinpath(sample_id + "_" + lane)
            lane_bam_file, cmd = bwa.bwa_mem(aione["config"], out_prefix, read_group_id, fq1, fq2)
            sample_bamfiles_by_lane[sample_id].append([lane_bam_file, cmd])

    bwa_shell_files_list = []
    aione["sample_final_sorted_bam"] = []
    for sample, sample_outdir in samples:
        sample_final_bamfile = sample_outdir.joinpath(f"{sample}.sorted.bam")
        aione["sample_final_sorted_bam"].append([sample, sample_final_bamfile])
        sample_shell_fname = shell_dirtory.joinpath(f"{sample}.bwa.sh")
        bwa_shell_files_list.append([sample, sample_shell_fname])

        if len(sample_bamfiles_by_lane[sample]) == 1:
            # single lane no need to merge
            lane_bam_file = sample_bamfiles_by_lane[sample][0][0]
            cmd = [sample_bamfiles_by_lane[sample][0][1]]
            if sample_final_bamfile != lane_bam_file:
                cmd.append(f"mv -f {lane_bam_file} {sample_final_bamfile}")

        else:
            samtools = aione["config"]["samtools"]["samtools"]
            samtools_merge_options = " ".join(
                [str(x) for x in aione["config"]["samtools"].get("merge_options", [])]
            )

            lane_bam_files = " ".join([f for f, _ in sample_bamfiles_by_lane[sample]])
            cmd = [c for _, c in sample_bamfiles_by_lane[sample]]

            # merge lane_bam_files together to be one big BAM file and rm lane_bam_files
            cmd.append(f"{samtools} merge {samtools_merge_options} {sample_final_bamfile} "
                       f"{lane_bam_files} && rm -rf {lane_bam_files}")

        echo_mark_done = f"echo \"[bwa] {sample} done\""
        cmd.append(echo_mark_done)

        if not is_dry_run and (not sample_shell_fname.exists() or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return bwa_shell_files_list  # [[sample, bwa_shell_file], ...]


def gatk_markduplicates(kwargs, out_folder_name: str, aione: dict = None,
                        is_dry_run: bool = False):
    """Create shell scripts for Markduplicates by GATK4.
    """
    shell_dirtory = Path(kwargs.outdir).joinpath(out_folder_name, "shell", "markdup")
    if not is_dry_run:
        safe_makedir(shell_dirtory)

    aione["sample_final_markdup_bam"] = []
    markdup_shell_files_list = []
    for sample, sample_sorted_bam in aione["sample_final_sorted_bam"]:
        # ``f_name_stem``: The final path component, minus its last suffix.
        dirname, f_name_stem = Path(sample_sorted_bam).parent, Path(sample_sorted_bam).stem

        # Setting Output path of markduplicate BAM file as the same as ``sample_sorted_bam``
        out_markdup_bam_fname = dirname.joinpath(f"{f_name_stem}.markdup.bam")
        out_markdup_metrics_fname = dirname.joinpath(f"{f_name_stem}.metrics.txt")
        sample_shell_fname = shell_dirtory.joinpath(f"{sample}.markdup.sh")

        markdup_shell_files_list.append([sample, sample_shell_fname])
        aione["sample_final_markdup_bam"].append([sample, out_markdup_bam_fname])

        cmd = [gatk.markduplicates(aione["config"], sample_sorted_bam, out_markdup_bam_fname,
                                   out_markdup_metrics_fname)]

        if IS_RM_SUBBAM:
            cmd.append(f"rm -rf {sample_sorted_bam}")  # save disk space

        echo_mark_done = f"echo \"[MarkDuplicates] {sample} done\""
        cmd.append(echo_mark_done)

        if not is_dry_run and (not sample_shell_fname.exists() or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return markdup_shell_files_list  # [[sample, sample_shell_fname], ...]


def gatk_baserecalibrator(kwargs, out_folder_name: str, aione: dict = None,
                          is_calculate_summary: bool = True,
                          is_dry_run: bool = False):
    shell_dirtory = Path(kwargs.outdir).joinpath(out_folder_name, "shell", "bqsr")
    if not is_dry_run:
        safe_makedir(shell_dirtory)

    aione["sample_final_bqsr_bam"] = []
    bqsr_shell_files_list = []

    is_calculate_contamination = True if "verifyBamID2" in aione["config"] else False
    for sample, sample_markdup_bam in aione["sample_final_markdup_bam"]:
        dirname, f_name_stem = Path(sample_markdup_bam).parent, Path(sample_markdup_bam).stem

        out_bqsr_bam_fname = dirname.joinpath(f"{f_name_stem}.BQSR.bam")
        out_bqsr_bai_fname = dirname.joinpath(f"{f_name_stem}.BQSR.bai")
        out_bqsr_recal_table = dirname.joinpath(f"{f_name_stem}.recal.table")

        out_alignment_summary_metric = dirname.joinpath(f"{f_name_stem}.AlignmentSummaryMetrics.txt")

        out_bamstats_fname = dirname.joinpath(f"{f_name_stem}.BQSR.stats")
        genome_cvg_fname = dirname.joinpath(f"{f_name_stem}.BQSR.depth.bed.gz")
        # when convert to CRAM format
        out_cram_fname = dirname.joinpath(f"{f_name_stem}.BQSR.cram")

        aione["sample_final_bqsr_bam"].append(
            [sample, out_bqsr_bam_fname if not kwargs.cram else out_cram_fname]
        )
        sample_shell_fname = shell_dirtory.joinpath(f"{sample}.bqsr.sh")
        bqsr_shell_files_list.append([sample, sample_shell_fname])

        cmd = [gatk.baserecalibrator(aione["config"],
                                     sample_markdup_bam,
                                     out_bqsr_bam_fname,
                                     out_bqsr_recal_table)]
        if IS_RM_SUBBAM:
            cmd.append(f"rm -rf {sample_markdup_bam}")

        if is_calculate_summary:
            cmd.append(gatk.collect_alignment_summary_metrics(
                aione["config"], out_bqsr_bam_fname, out_alignment_summary_metric)
            )
            cmd.append(bam.stats(aione["config"], out_bqsr_bam_fname, out_bamstats_fname))
            cmd.append(bam.genomecoverage(aione["config"], out_bqsr_bam_fname, genome_cvg_fname))

        if is_calculate_contamination:
            out_verifybamid_stat_prefix = dirname.joinpath(f"{f_name_stem}.BQSR.verifyBamID2")
            cmd.append(bam.verifyBamID2(aione["config"], out_bqsr_bam_fname, out_verifybamid_stat_prefix))

        if kwargs.cram:
            cmd.append(bwa.bam_to_cram(aione["config"], out_bqsr_bam_fname, out_cram_fname))
            cmd.append(f"rm -rf {out_bqsr_bam_fname}")
            cmd.append(f"rm -rf {out_bqsr_bai_fname}")

        echo_mark_done = f"echo \"[BQSR] {sample} done\""
        cmd.append(echo_mark_done)

        if not is_dry_run and (not sample_shell_fname.exists() or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname, cmd)

    return bqsr_shell_files_list


def gatk_haplotypecaller_gvcf(kwargs, out_folder_name: str, aione: dict = None,
                              is_dry_run: bool = False):
    """Create gvcf shell.
    """

    def _create_sub_shell(sample, sample_shell_dir, sample_output_dir,
                          raw_interval=None, is_dry_run=False):

        interval_ = raw_interval if raw_interval else "all"
        # in case the raw_interval is a full path file.
        interval_ = Path(interval_).name
        sample_shell_fname_ = sample_shell_dir.joinpath(f"{sample}.{interval_}.gvcf.sh")
        out_gvcf_fname_ = sample_output_dir.joinpath(f"{sample}.{interval_}.g.vcf.gz")

        if raw_interval:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname_,
                                             interval=raw_interval)]
        else:
            cmd = [gatk.haplotypecaller_gvcf(aione["config"], sample_bqsr_bam, out_gvcf_fname_)]
        echo_mark_done = f"echo \"[GVCF] {sample} {interval_} done\""
        cmd.append(echo_mark_done)

        if (not is_dry_run) and (not sample_shell_fname_.exists() or kwargs.overwrite):
            _create_cmd_file(sample_shell_fname_, cmd)

        return sample_shell_fname_, out_gvcf_fname_

    output_dirtory, shell_dirtory = _md(Path(kwargs.outdir).joinpath(out_folder_name),
                                        is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    gvcf_shell_files_list = []
    aione["gvcf"] = {}

    if "interval" not in aione["config"]["gatk"]:
        aione["config"]["gatk"]["interval"] = ["all"]

    aione["intervals"] = []  # chromosome id
    sample_map = {}  # Record sample_map for genomicsDBImport
    for sample, sample_bqsr_bam in aione["sample_final_bqsr_bam"]:

        sample_output_dir = output_dirtory.joinpath(sample)
        sample_shell_dir = shell_dirtory.joinpath(sample)
        if not is_dry_run:
            safe_makedir(sample_output_dir)
            safe_makedir(sample_shell_dir)

        for interval in aione["config"]["gatk"]["interval"]:
            if interval == "all":
                # All the whole genome once, take a lot of time, not suggested
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, is_dry_run=is_dry_run)
            else:
                sample_shell_fname, out_gvcf_fname = _create_sub_shell(
                    sample, sample_shell_dir, sample_output_dir, raw_interval=interval,
                    is_dry_run=is_dry_run)

            # ``aione["config"]["gatk"]["interval"]`` may be a file path.
            # Make sure the interval id does not contain any path if the raw interval is a file path
            interval = Path(interval).name
            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            gvcf_shell_files_list.append([sample + ".%s" % interval, sample_shell_fname])
            aione["gvcf"][interval].append(out_gvcf_fname)

            if is_use_gDBI and (interval not in sample_map):
                sample_map[interval] = []

            if is_use_gDBI:
                sample_map[interval].append(f"{sample}\t{out_gvcf_fname}")

    if is_use_gDBI:
        aione["sample_map"] = {}
        for interval, value in sample_map.items():
            # create sample_map file for next process
            out_sample_map_fname = output_dirtory.joinpath(f"{interval}.gvcf.sample_map")
            aione["sample_map"][interval] = out_sample_map_fname

            if (not is_dry_run) and (not out_sample_map_fname.exists() or kwargs.overwrite):
                with open(str(out_sample_map_fname), "w") as OUT:
                    OUT.write("\n".join(value))

    return gvcf_shell_files_list


def _yield_gatk_combineGVCFs(project_name, variant_calling_intervals, output_dirtory: Path,
                             shell_dirtory: Path, is_use_gDBI: bool, aione: dict = None):
    for interval in variant_calling_intervals:
        interval_id = interval[0] if type(interval) is list else interval  # get chromosome id
        if interval_id in aione["gvcf"]:
            # load gvcf '.sample_map' file if using genomicsDBImport module to combine all the gvcf
            sample_gvcf_list = [aione["sample_map"][interval_id]] \
                if is_use_gDBI else aione["gvcf"][interval_id]
        else:
            continue

        if not sample_gvcf_list:
            sys.stderr.write("[Error] Interval parameters in configuration file may be different "
                             "from that input gvcf file list in ``gatk_combineGVCFs`` function.")
            sys.exit(1)

        interval_n = "_".join(interval) if type(interval) is list else interval.replace(":", "_")  # chr:s -> chr_s

        sub_shell_fname = shell_dirtory.joinpath(f"{project_name}.{interval_n}.combineGVCFs.sh")
        combineGVCF_fname = output_dirtory.joinpath(f"{project_name}.{interval_n}.genomics_db") \
            if is_use_gDBI else output_dirtory.joinpath(f"{project_name}.{interval_n}.g.vcf.gz")

        # Create commandline for combining GVCFs
        calling_interval = "%s:%s-%s" % (interval[0], interval[1], interval[2]) \
            if type(interval) is list else interval

        combineGVCFs_cmd = gatk.combineGVCFs(aione["config"],
                                             sample_gvcf_list,
                                             combineGVCF_fname,
                                             interval=calling_interval)

        yield combineGVCFs_cmd, combineGVCF_fname, sub_shell_fname, interval_n, calling_interval


def gatk_combineGVCFs(kwargs, out_folder_name: str, aione: dict = None, is_dry_run: bool = False):
    """Combine GVCFs by GATK genomicsDBImport or CombineGVCFs module.
    """
    output_dirtory, shell_dirtory = _md(Path(kwargs.outdir).joinpath(out_folder_name),
                                        is_dry_run=is_dry_run)
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

        echo_mark_done = f"echo \"[CombineGVCFs] {calling_interval} done\""
        cmd = [combineGVCFs_cmd, echo_mark_done]

        if (not is_dry_run) and (not sub_shell_fname.exists() or kwargs.overwrite):
            _create_cmd_file(sub_shell_fname, cmd)

    return combineGVCFs_shell_files_list


def gatk_genotypeGVCFs(kwargs, out_folder_name: str, aione: dict = None, is_dry_run: bool = False):
    """Create shell for genotypeGVCFs.
    """
    output_dirtory, shell_dirtory = _md(Path(kwargs.outdir).joinpath(out_folder_name),
                                        is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    genotype_vcf_shell_files_list = []
    aione["genotype_vcf_list"] = []

    variant_calling_intervals = aione["config"]["gatk"]["variant_calling_interval"]
    for interval in variant_calling_intervals:
        interval_n = "_".join(interval) if type(interval) is list else interval.replace(":", "_")  # chr:s -> chr_s

        if interval_n not in aione["combineGVCFs"]:
            continue

        sub_shell_fname = shell_dirtory.joinpath(f"{kwargs.project_name}.{interval_n}.genotype.sh")
        genotype_vcf_fname = str(output_dirtory.joinpath(f"{kwargs.project_name}.{interval_n}.vcf.gz"))
        genotype_vcf_shell_files_list.append([kwargs.project_name + "." + interval_n, sub_shell_fname])
        aione["genotype_vcf_list"].append(genotype_vcf_fname)

        calling_interval = "%s:%s-%s" % (interval[0], interval[1], interval[2]) \
            if type(interval) is list else interval

        # Create commandline for genotype joint-calling process
        combineGVCF_fname = aione["combineGVCFs"][interval_n]
        genotype_cmd = gatk.genotypeGVCFs(aione["config"],
                                          combineGVCF_fname,
                                          genotype_vcf_fname,
                                          interval=calling_interval)

        # delete the genomicsdb workspace (or the combine GVCF file).
        # delete_cmd = f"rm -rf {combineGVCF_fname}" \
        #     if is_use_gDBI else f"rm -rf {combineGVCF_fname} {combineGVCF_fname}.tbi"

        echo_mark_done = f"echo \"[Genotype] {calling_interval} done\""
        cmd = [genotype_cmd, echo_mark_done]  # [genotype_cmd, delete_cmd, echo_mark_done]

        if (not is_dry_run) and (not sub_shell_fname.exists() or kwargs.overwrite):
            _create_cmd_file(sub_shell_fname, cmd)

    if not genotype_vcf_shell_files_list:
        sys.stderr.write("[Error] Interval parameters in configuration file may be different "
                         "from that input gvcf file list when calls ``gatk_genotypeGVCFs``.")
        sys.exit(1)

    return genotype_vcf_shell_files_list


def gatk_genotype(kwargs, out_folder_name: str, aione: dict = None, is_dry_run: bool = False):
    """CombineGVCFs and genotypeGVCFs into one function (not use).
    """
    output_dirtory, shell_dirtory = _md(Path(kwargs.outdir).joinpath(out_folder_name),
                                        is_dry_run=is_dry_run)
    is_use_gDBI = aione["config"]["gatk"]["use_genomicsDBImport"] \
        if "use_genomicsDBImport" in aione["config"]["gatk"] else False

    genotype_vcf_shell_files_list = []
    aione["genotype_vcf_list"] = []

    # Create commandline for combineGVCFs and genotypeGVCFs process for specific calling intervals
    variant_calling_intervals = aione["config"]["gatk"]["variant_calling_interval"]
    for (combineGVCFs_cmd,
         combineGVCF_fname,
         _,  # no need the sub_shell_fname in _yield_gatk_combineGVCFs
         interval_n,  # interval marker
         calling_interval) in _yield_gatk_combineGVCFs(kwargs.project_name,
                                                       variant_calling_intervals,
                                                       output_dirtory,
                                                       shell_dirtory,
                                                       is_use_gDBI,
                                                       aione):

        sub_shell_fname = shell_dirtory.joinpath(f"{kwargs.project_name}.{interval_n}.genotype.sh")
        genotype_vcf_fname = output_dirtory.joinpath(f"{kwargs.project_name}.{interval_n}.vcf.gz")

        genotype_vcf_shell_files_list.append([kwargs.project_name + "." + interval_n, sub_shell_fname])
        aione["genotype_vcf_list"].append(genotype_vcf_fname)

        # Generate the genotype joint-calling process according to the input file
        # ``combineGVCF_fname`` which create in  ``_yield_gatk_CombineGVCFs``
        genotype_cmd = gatk.genotypeGVCFs(aione["config"],
                                          combineGVCF_fname,
                                          genotype_vcf_fname,
                                          interval=calling_interval)

        # delete the genomicsdb workspace (or the combine GVCF file).
        delete_cmd = f"rm -rf {combineGVCF_fname}" \
            if is_use_gDBI else f"rm -rf {combineGVCF_fname} {combineGVCF_fname}.tbi"

        echo_mark_done = f"echo \"[Genotype] {calling_interval} done\""

        # combine all the commands for the final genotype calling.
        cmd = [combineGVCFs_cmd, genotype_cmd, delete_cmd, echo_mark_done]
        if (not is_dry_run) and (not sub_shell_fname.exists() or kwargs.overwrite):
            # Create shell file
            _create_cmd_file(sub_shell_fname, cmd)

    return genotype_vcf_shell_files_list


def gatk_variantrecalibrator(kwargs, out_folder_name: str, aione: dict = None,
                             is_dry_run: bool = False):
    """Create shell scripts for VQSR.
    """
    output_dirtory, shell_dirtory = _md(Path(kwargs.outdir).joinpath(out_folder_name),
                                        is_dry_run=is_dry_run)
    genotype_vqsr_fname = str(output_dirtory.joinpath(f"{kwargs.project_name}.VQSR.vcf.gz"))
    combine_vcf_fname = str(output_dirtory.joinpath(f"{kwargs.project_name}.raw.vcf.gz"))
    shell_fname = shell_dirtory.joinpath(f"{kwargs.project_name}.VQSR.sh")

    cmd = []
    if len(aione["genotype_vcf_list"]) > 1:
        # concat-vcf
        concat_vcf_cmd = bcftools.concat(aione["config"],
                                         aione["genotype_vcf_list"],
                                         combine_vcf_fname)
        cmd.append(concat_vcf_cmd)
    else:
        combine_vcf_fname = str(aione["genotype_vcf_list"][0])

    # VQSR process
    cmd.append(gatk.variantrecalibrator(aione["config"], combine_vcf_fname, genotype_vqsr_fname))
    cmd.append(f"echo \"[VQSR] {genotype_vqsr_fname} done\"")

    if not is_dry_run and (not shell_fname.exists() or kwargs.overwrite):
        _create_cmd_file(shell_fname, cmd)

    # Only one VQSR result
    return [[f"{kwargs.project_name}.VQSR", shell_fname]]
