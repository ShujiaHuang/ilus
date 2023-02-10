"""Pipeline functions.

Author: Shujia Huang
Date: 2020-04-19

"""
import os
import stat
import sys

from ilus.utils import safe_makedir
from ilus.launch.runfunction import bwamem, gatk_markduplicates, gatk_baserecalibrator, \
    gatk_haplotypecaller_gvcf, gatk_combineGVCFs, gatk_genotypeGVCFs, gatk_variantrecalibrator


def _create_a_total_shell_file(shell_list, out_shell_filename, sub_shell_log_dir, o_log_file, e_log_file):
    """Creat all the executable in a big shell script which gather all scripts from ``shell_list``.

    ``shell_list`` is a 2-D array: [[mark, shell_file], ...].
    """
    sub_shell_log_dir = sub_shell_log_dir.rstrip("/")
    with open(out_shell_filename, "w") as OUT, open(o_log_file, "w") as O_LOG, open(e_log_file, "w") as E_LOG:
        OUT.write("#!/bin/bash\n")
        for marker, sub_shell in shell_list:
            OUT.write("{sub_shell} 2> {sub_shell_log_dir}/{marker}.e.log > "
                      "{sub_shell_log_dir}/{marker}.o.log\n".format(**locals()))

            # record all the path of log files into a single file
            O_LOG.write("{marker}\t{sub_shell_log_dir}/{marker}.o.log\t{sub_shell}\n".format(**locals()))
            E_LOG.write("{marker}\t{sub_shell_log_dir}/{marker}.e.log\t{sub_shell}\n".format(**locals()))

    os.chmod(out_shell_filename, stat.S_IRWXU)  # 0700
    return


def _make_process_shell(output_shell_fname,
                        shell_log_directory,
                        process_shells=None,
                        is_overwrite=False,
                        is_dry_run=False):
    if is_dry_run:
        return

    safe_makedir(shell_log_directory)
    o_log_file = shell_log_directory.rstrip("/") + ".o.log.list"
    e_log_file = shell_log_directory.rstrip("/") + ".e.log.list"

    if not is_overwrite and os.path.exists(output_shell_fname):
        sys.stderr.write("%s is already exist. Please set -f parameter if "
                         "you want to overwrite.\n" % output_shell_fname)
        return

    _create_a_total_shell_file(process_shells, output_shell_fname, shell_log_directory, o_log_file, e_log_file)
    return


def wgs(kwargs, aione):
    # All the WGS processes.
    runner_module = {
        # [func, shell_file, shell_log_folder, output_folder]
        # create bwa/sort/merge process
        "align": [bwamem, kwargs.project_name + ".step1.bwa.sh", "01.alignment", "01.alignment"],

        # Create Markduplicates shells.
        "markdup": [gatk_markduplicates, kwargs.project_name + ".step2.markdup.sh", "02.markdup", "01.alignment"],

        # Create BQSR+ApplyBQSR shells.
        "BQSR": [gatk_baserecalibrator, kwargs.project_name + ".step3.bqsr.sh", "03.BQSR", "01.alignment"],

        # Create GVCF shells
        "gvcf": [gatk_haplotypecaller_gvcf, kwargs.project_name + ".step4.gvcf.sh", "04.gvcf", "02.gvcf"],

        # combineGVCFs
        "combineGVCFs": [gatk_combineGVCFs, kwargs.project_name + ".step5.combineGVCFs.sh", "05.combineGVCFs", "03.genotype"],

        # GenotypeGVCF
        "genotype": [gatk_genotypeGVCFs, kwargs.project_name + ".step6.genotype.sh", "06.genotype", "03.genotype"],

        # Variant recalibrator
        "VQSR": [gatk_variantrecalibrator, kwargs.project_name + ".step7.VQSR.sh", "07.VQSR", "03.genotype"],

        # Todo: Integrate summary and status statistic information into ilus pipeline.
        "summary": []
    }

    # Todo: Need a program to validate whether the tools, arguments and processes are
    #  available or not for the pipeline.

    wgs_processes_order = ["align", "markdup", "BQSR", "gvcf", "combineGVCFs", "genotype", "VQSR"]
    processes_set = set(kwargs.wgs_processes.split(","))
    for p in processes_set:
        if p not in wgs_processes_order:
            sys.stderr.write("[ERROR] %s is not one of the wgs processes: %s\n" % (p, ",".join(wgs_processes_order)))
            sys.exit(1)

    # Create project directory and return the abspath.
    # [Important] abspath will remove the last '/' in the path. e.g.: '/a/' -> '/a'
    kwargs.outdir = os.path.abspath(kwargs.outdir)
    shell_dirtory = os.path.join(kwargs.outdir, "00.shell")
    shell_log_dirtory = os.path.join(shell_dirtory, "loginfo")

    if not kwargs.dry_run:
        safe_makedir(kwargs.outdir)
        safe_makedir(shell_dirtory)
        safe_makedir(shell_log_dirtory)

    for p in wgs_processes_order:
        is_dry_run = True if kwargs.dry_run or (p not in processes_set) else False

        func, shell_fname, shell_log_folder, output_result_folder = runner_module[p]
        _make_process_shell(output_shell_fname=os.path.join(shell_dirtory, shell_fname),
                            shell_log_directory=os.path.join(shell_log_dirtory, shell_log_folder),
                            process_shells=func(kwargs, output_result_folder, aione, is_dry_run=is_dry_run),
                            is_overwrite=kwargs.overwrite,
                            is_dry_run=is_dry_run)

    return aione


def _f(kwargs, aione, shell_fname, shell_log_folder, function_name):
    # [Important] abspath will remove the last '/' in the path. e.g.: '/a/' -> '/a'
    kwargs.outdir = safe_makedir(os.path.abspath(kwargs.outdir))  # return abspath
    root_path, output_folder_name = os.path.split(kwargs.outdir)

    tmp_dir = kwargs.outdir
    kwargs.outdir = root_path

    shell_dirtory = os.path.join(root_path, "00.shell" if kwargs.as_pipe_shell_order else "shell")
    if not os.path.exists(shell_dirtory):
        safe_makedir(shell_dirtory)

    _make_process_shell(output_shell_fname=os.path.join(shell_dirtory, shell_fname),
                        shell_log_directory=os.path.join(shell_dirtory, "loginfo", shell_log_folder),
                        process_shells=function_name(kwargs, output_folder_name, aione),
                        is_overwrite=kwargs.overwrite)

    kwargs.outdir = tmp_dir  # set back
    return


def genotypeGVCFs(kwargs, aione):
    """GenotypeGVCFs by GATK."""

    aione["intervals"] = []
    aione["gvcf"] = {}  # will be called in ``gatk_genotypeGVCFs``
    sample_map = {}     # Record sample_map for genomicsDBImport
    with open(kwargs.gvcflist) as I:
        # Format in gvcfilist: [chromosome_id  sampleid  gvcf_file_path]
        for line in I:
            if line.startswith("#"):
                continue

            try:
                interval, sample, gvcf = line.strip().split()
            except ValueError:
                raise ValueError("[ERROR] format error in '%s'. \nExpected: 3 columns "
                                 "(column 1: chromosome ID, column 2: sample ID, column 3: gvcf file path), "
                                 "got '%s')." % (kwargs.gvcflist, line.strip()))

            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            aione["gvcf"][interval].append(gvcf)
            if interval not in sample_map:
                sample_map[interval] = []

            sample_map[interval].append("%s\t%s" % (sample, gvcf))

    if aione["config"]["gatk"]["use_genomicsDBImport"]:
        # [Important] return abspath and abspath will remove the last '/' in the path. e.g.: '/a/' -> '/a'
        kwargs.outdir = safe_makedir(os.path.abspath(kwargs.outdir))
        aione["sample_map"] = {}
        for interval, value in sample_map.items():
            out_sample_map_fname = os.path.join(kwargs.outdir, "%s.gvcf.sample_map" % interval)
            aione["sample_map"][interval] = out_sample_map_fname
            with open(out_sample_map_fname, "w") as OUT:
                OUT.write("\n".join(value))

    # create shell scripts for genomicsDBImport or CombineGVCFs before the genotype joint-calling process
    combineGVCF_shell_fname, combineGVCF_shell_log_folder = [
        kwargs.project_name + ".step5.combineGVCFs.sh", "05.combineGVCFs"] \
        if kwargs.as_pipe_shell_order else [kwargs.project_name + ".combineGVCFs.sh", "combineGVCFs"]
    _f(kwargs, aione, combineGVCF_shell_fname, combineGVCF_shell_log_folder, gatk_combineGVCFs)

    # create shell scripts for genotype joint-calling
    genotype_shell_fname, genotype_shell_log_folder = [
        kwargs.project_name + ".step6.genotype.sh", "06.genotype"] \
        if kwargs.as_pipe_shell_order else [kwargs.project_name + ".genotype.sh", "genotype"]
    _f(kwargs, aione, genotype_shell_fname, genotype_shell_log_folder, gatk_genotypeGVCFs)

    return aione


def variantrecalibrator(kwargs, aione):
    aione["genotype_vcf_list"] = []  # will be called in ``gatk_variantrecalibrator``
    with open(kwargs.vcflist) as I:
        # Format in vcfilist one file per row
        for line in I:
            if line.startswith("#"):
                continue

            aione["genotype_vcf_list"].append(line.strip().split()[0])

    shell_fname, shell_log_folder = [kwargs.project_name + ".step6.VQSR.sh", "06.VQSR"] \
        if kwargs.as_pipe_shell_order else [kwargs.project_name + ".vqsr.sh", "genotype"]

    _f(kwargs, aione, shell_fname, shell_log_folder, gatk_variantrecalibrator)
    return aione
