"""Pipeline functions.

Author: Shujia Huang
Date: 2020-04-19

"""
import stat
import sys
from pathlib import Path
from typing import List, Tuple

from ilus.modules.utils import safe_makedir, file_exists, check_samplesheet
from ilus.launch.runfunction import (
    bwamem,
    gatk_markduplicates,
    gatk_baserecalibrator,
    gatk_haplotypecaller_gvcf,
    gatk_combineGVCFs,
    gatk_genotypeGVCFs,
    gatk_variantrecalibrator
)


def _create_a_total_shell_file(shell_list: List[Tuple[str, str]],
                               out_shell_filename: Path,
                               sub_shell_log_dir: Path,
                               o_log_file: str,
                               e_log_file: str) -> None:
    """Creat all the executable in a big shell script which gather all scripts from ``shell_list``.

        ``shell_list`` is a 2-D array: [[mark, shell_file], ...].
    """
    sub_shell_log_dir = str(sub_shell_log_dir).rstrip("/")
    o_log_template = "{marker}\t{sub_shell_log_dir}/{marker}.o.log\t{sub_shell}\n"
    e_log_template = "{marker}\t{sub_shell_log_dir}/{marker}.e.log\t{sub_shell}\n"
    shell_template = "{sub_shell} 2> {sub_shell_log_dir}/{marker}.e.log > {sub_shell_log_dir}/{marker}.o.log\n"

    with open(str(out_shell_filename), "w") as out_file, open(o_log_file, "w") as o_log, open(e_log_file, "w") as e_log:
        out_file.write("#!/bin/bash\n")
        for marker, sub_shell in shell_list:
            # record all the path of log files into a single file
            o_log.write(o_log_template.format(marker=marker, sub_shell_log_dir=sub_shell_log_dir, sub_shell=sub_shell))
            e_log.write(e_log_template.format(marker=marker, sub_shell_log_dir=sub_shell_log_dir, sub_shell=sub_shell))

            out_file.write(
                shell_template.format(marker=marker, sub_shell=sub_shell, sub_shell_log_dir=sub_shell_log_dir))

    out_shell_filename.chmod(stat.S_IRWXU)  # 0700
    return


def _make_process_shell(output_shell_fname: Path,
                        shell_log_directory: Path,
                        process_shells: List[Tuple[str, str]] = None,
                        is_overwrite: bool = False,
                        is_dry_run: bool = False) -> None:
    if is_dry_run:
        return

    safe_makedir(shell_log_directory)
    o_log_file = str(shell_log_directory).rstrip("/") + ".o.log.list"
    e_log_file = str(shell_log_directory).rstrip("/") + ".e.log.list"

    if not is_overwrite and Path(output_shell_fname).exists():
        sys.stderr.write(f"{output_shell_fname} is already exist. Please set '-f' "
                         f"parameter if you want to overwrite.\n")
        return

    _create_a_total_shell_file(process_shells,
                               output_shell_fname,
                               shell_log_directory,
                               o_log_file,
                               e_log_file)
    return


def WGS(kwargs, aione: dict = None) -> dict:
    """Create the WGS data analysis pipeline."""

    runner_module = {
        # [func, shell_file_name, shell_log_folder_name, output_folder_name]
        # create bwa/sort/merge process
        "align": [bwamem, f"{kwargs.project_name}.step1.bwa.sh", "01.alignment", "01.alignment"],

        # Create Markduplicates shells.
        "markdup": [gatk_markduplicates, f"{kwargs.project_name}.step2.markdup.sh", "02.markdup",
                    "01.alignment"],

        # Create BQSR+ApplyBQSR shells.
        "BQSR": [gatk_baserecalibrator, f"{kwargs.project_name}.step3.bqsr.sh", "03.BQSR", "01.alignment"],

        # Create GVCF shells
        "gvcf": [gatk_haplotypecaller_gvcf, f"{kwargs.project_name}.step4.gvcf.sh", "04.gvcf", "02.gvcf"],

        # combineGVCFs
        "combineGVCFs": [gatk_combineGVCFs, f"{kwargs.project_name}.step5.combineGVCFs.sh", "05.combineGVCFs",
                         "03.genotype"],

        # GenotypeGVCF
        "genotype": [gatk_genotypeGVCFs, f"{kwargs.project_name}.step6.genotype.sh", "06.genotype",
                     "03.genotype"],

        # Variant recalibrator
        "VQSR": [gatk_variantrecalibrator, f"{kwargs.project_name}.step7.VQSR.sh", "07.VQSR",
                 "03.genotype"],

        # Todo: Integrate summary and status statistic information of ilus pipeline.
        "summary": []
    }

    # Todo: Need a program to validate whether the tools, arguments and processes are
    #  available or not for the pipeline.

    wgs_pipeline_processes = ["align", "markdup", "BQSR", "gvcf", "combineGVCFs", "genotype", "VQSR"]
    processes_set = set(kwargs.wgs_processes.split(","))
    for p in processes_set:
        if p not in wgs_pipeline_processes:
            sys.stderr.write(f"[ERROR] {p} is not one of the wgs processes: "
                             f"{','.join(wgs_pipeline_processes)}\n")
            sys.exit(1)

    if not file_exists(kwargs.fastqlist):
        print(f"[ERROR] {kwargs.fastqlist} is not a file or empty.\n", file=sys.stderr)
        sys.exit(1)

    if check_samplesheet(kwargs.fastqlist):
        print(f"\n[INFO] Input file format of '{kwargs.fastqlist}' has been verified.",
              file=sys.stderr)

    aione["fastqlist"] = kwargs.fastqlist  # record the input file path of fastqlist.

    # Create project directory and return the abspath.
    # [Important] abspath will remove the last '/' of the path. e.g.: '/a/' -> '/a'
    kwargs.outdir = Path(kwargs.outdir).resolve()
    shell_dirtory = kwargs.outdir.joinpath("00.shell")
    shell_log_dirtory = shell_dirtory.joinpath("loginfo")

    if not kwargs.dry_run:
        safe_makedir(kwargs.outdir)
        safe_makedir(shell_dirtory)
        safe_makedir(shell_log_dirtory)

    for p in wgs_pipeline_processes:
        is_dry_run = True if kwargs.dry_run or (p not in processes_set) else False

        func, shell_fname, shell_log_folder, output_folder = runner_module[p]
        _make_process_shell(output_shell_fname=shell_dirtory.joinpath(shell_fname),
                            shell_log_directory=shell_log_dirtory.joinpath(shell_log_folder),
                            process_shells=func(kwargs, output_folder, aione, is_dry_run=is_dry_run),
                            is_overwrite=kwargs.overwrite,
                            is_dry_run=is_dry_run)

    return aione


def _f(kwargs, aione, shell_fname, shell_log_folder, function_name):
    # [Important] abspath will remove the last '/' in the path. e.g.: '/a/' -> '/a'
    kwargs.outdir = safe_makedir(kwargs.outdir.resolve())  # return abspath
    root_path, output_folder_name = kwargs.outdir.parent, kwargs.outdir.name

    tmp_dir = kwargs.outdir
    kwargs.outdir = root_path

    shell_dirtory = root_path.joinpath("00.shell" if kwargs.as_pipe_shell_order else "shell")
    safe_makedir(shell_dirtory)

    _make_process_shell(output_shell_fname=shell_dirtory.joinpath(shell_fname),
                        shell_log_directory=shell_dirtory.joinpath("loginfo", shell_log_folder),
                        process_shells=function_name(kwargs, output_folder_name, aione),
                        is_overwrite=kwargs.overwrite)

    kwargs.outdir = tmp_dir  # set back
    return


def genotypeGVCFs(kwargs, aione: dict = None) -> dict:
    """GenotypeGVCFs by GATK."""

    aione["intervals"] = []
    aione["gvcf"] = {}  # will be called in ``gatk_genotypeGVCFs``
    sample_map = {}  # Record sample_map for genomicsDBImport
    with open(kwargs.gvcflist) as I:
        # Format in gvcfilist: [chromosome_id  sampleid  gvcf_file_path]
        for line in I:
            if line.startswith("#"):
                continue

            try:
                interval, sample, gvcf = line.strip().split()
            except ValueError:
                raise ValueError(f"Input format error in '{kwargs.gvcflist}'. Expected: 3 columns \n\n"
                                 f"-- column 1: chromosome ID \n"
                                 f"-- column 2: sample ID\n"
                                 f"-- column 3: gvcf file path)\n\n"
                                 f"got error here: '{line.strip()}').")

            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            aione["gvcf"][interval].append(gvcf)
            if interval not in sample_map:
                sample_map[interval] = []

            sample_map[interval].append(f"{sample}\t{gvcf}")

    if aione["config"]["gatk"]["use_genomicsDBImport"]:
        # [Important] return abspath and abspath will remove the last '/' in the path. e.g.: '/a/' -> '/a'
        kwargs.outdir = Path(kwargs.outdir).resolve()
        safe_makedir(kwargs.outdir)

        aione["sample_map"] = {}
        for interval, value in sample_map.items():
            out_sample_map_fname = kwargs.outdir.joinpath(f"{interval}.gvcf.sample_map")
            aione["sample_map"][interval] = out_sample_map_fname
            with open(out_sample_map_fname, "w") as OUT:
                OUT.write("\n".join(value))

    # create shell scripts for genomicsDBImport or CombineGVCFs before the genotype joint-calling process
    combineGVCF_shell_fname, combineGVCF_shell_log_folder = [
        f"{kwargs.project_name}.step5.combineGVCFs.sh", "05.combineGVCFs"] \
        if kwargs.as_pipe_shell_order else [f"{kwargs.project_name}.combineGVCFs.sh", "combineGVCFs"]

    _f(kwargs, aione, combineGVCF_shell_fname, combineGVCF_shell_log_folder, gatk_combineGVCFs)

    # create shell scripts for genotype joint-calling
    genotype_shell_fname, genotype_shell_log_folder = [
        f"{kwargs.project_name}.step6.genotype.sh", "06.genotype"
    ] if kwargs.as_pipe_shell_order else [f"{kwargs.project_name}.genotype.sh", "genotype"]

    _f(kwargs, aione, genotype_shell_fname, genotype_shell_log_folder, gatk_genotypeGVCFs)

    return aione


def variantrecalibrator(kwargs, aione: dict = None) -> dict:
    aione["genotype_vcf_list"] = []  # will be called in ``gatk_variantrecalibrator``
    with open(kwargs.vcflist) as f_in:
        # Format in vcfilist one file per row
        for line in f_in:
            if line.startswith("#"):
                continue

            aione["genotype_vcf_list"].append(line.strip().split()[0])

    shell_fname, shell_log_folder = [
        f"{kwargs.project_name}.step7.VQSR.sh", "07.VQSR"
    ] if kwargs.as_pipe_shell_order else [f"{kwargs.project_name}.vqsr.sh", "genotype"]

    _f(kwargs, aione, shell_fname, shell_log_folder, gatk_variantrecalibrator)
    return aione
