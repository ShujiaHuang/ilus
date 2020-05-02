"""Pipeline functions."""
import os
import stat

from ilus.utils import safe_makedir
from ilus.launch.runfunction import bwamem, gatk_markduplicates, gatk_baserecalibrator, \
    gatk_haplotypecaller_gvcf, gatk_genotypeGVCFs, gatk_variantrecalibrator


def _create_a_total_shell_file(shell_list, out_shell_filename, sub_shell_log_dir, o_log_file, e_log_file):
    """Creat all the executable a big shell script which gather all scripts from ``shell_list``.

    ``shell_list`` is a 2-D array: [[mark, shell_file], ...].
    """
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


def _make_process_shell(output_shell_fname, shell_log_directory, process_shells=None, is_overwrite=False):
    safe_makedir(shell_log_directory)
    o_log_file = shell_log_directory + ".o.log.list"
    e_log_file = shell_log_directory + ".e.log.list"

    if not is_overwrite and os.path.exists(output_shell_fname):
        print("%s is already exist. Please set -f parameter if you want to overwrite." % output_shell_fname)
        return

    _create_a_total_shell_file(process_shells, output_shell_fname, shell_log_directory, o_log_file, e_log_file)
    return


def wgs(kwargs, aione):
    # All the WGS processes.
    runner_module = {
        # [func, shell_file, shell_log_folder, output_folder]

        # bwa/sort/merge process
        "align": [bwamem, "step1.bwa.sh", "01.alignment", "01.alignment"],

        # Create Markduplicates running shells.
        "markdup": [gatk_markduplicates, "step2.markdup.sh", "02.markdup", "01.alignment"],

        # Create BQSR+ApplyBQSR running shells.
        "BQSR": [gatk_baserecalibrator, "step3.bqsr.sh", "03.BQSR", "01.alignment"],

        # Create GVCF running shells
        "gvcf": [gatk_haplotypecaller_gvcf, "step4.gvcf.sh", "04.gvcf", "02.gvcf"],

        # GenotypeGVCF
        "genotype": [gatk_genotypeGVCFs, "step5.genotype.sh", "05.genotype", "03.genotype"],

        # Variant recalibrator
        "VQSR": [gatk_variantrecalibrator, "step6.VQSR.sh", "06.VQSR", "03.genotype"],

        # Todo: Integrate summary and status statistic information into ilus pipeline.
        "summary": []
    }

    # Create project directory
    kwargs.outdir = safe_makedir(os.path.abspath(kwargs.outdir))  # return abspath
    shell_dirtory = os.path.join(kwargs.outdir, "00.shell")
    shell_log_dirtory = os.path.join(shell_dirtory, "loginfo")
    safe_makedir(shell_dirtory)
    safe_makedir(shell_log_dirtory)

    for p in ["align", "markdup", "BQSR", "gvcf", "genotype", "VQSR"]:
        func, shell_fname, shell_log_folder, output_result_folder = runner_module[p]
        _make_process_shell(output_shell_fname=os.path.join(shell_dirtory, shell_fname),
                            shell_log_directory=os.path.join(shell_log_dirtory, shell_log_folder),
                            process_shells=func(kwargs, output_result_folder, aione),
                            is_overwrite=kwargs.overwrite)

    return aione


def _f(kwargs, aione, shell_fname, shell_log_folder, function_name):

    kwargs.outdir = safe_makedir(os.path.abspath(kwargs.outdir))  # return abspath
    root_path, output_result_folder = os.path.split(kwargs.outdir)

    shell_dirtory = os.path.join(root_path, "00.shell" if kwargs.as_pipe_shell_order else "shell")
    if not os.path.exists(shell_dirtory):
        safe_makedir(shell_dirtory)

    shell_log_dirtory = os.path.join(shell_dirtory, "loginfo")
    safe_makedir(shell_log_dirtory)

    _make_process_shell(output_shell_fname=os.path.join(shell_dirtory, shell_fname),
                        shell_log_directory=os.path.join(shell_log_dirtory, shell_log_folder),
                        process_shells=function_name(kwargs, output_result_folder, aione),
                        is_overwrite=kwargs.overwrite)

    return


def genotypeGVCFs(kwargs, aione):
    """GenotypeGVCFs by GATK"""

    aione["intervals"] = []
    aione["gvcf"] = {}  # will be called in ``gatk_genotypeGVCFs``
    with open(kwargs.gvcflist) as I:
        # Format in gvcfilist: [Interval  gvcf_file_path]
        for line in I:
            if line.startswith("#"):
                continue

            interval, gvcf = line.strip().split()
            if interval not in aione["gvcf"]:
                aione["intervals"].append(interval)
                aione["gvcf"][interval] = []

            aione["gvcf"][interval].append(gvcf)

    shell_fname, shell_log_folder = ["step5.genotype.sh", "05.genotype"] if kwargs.as_pipe_shell_order \
        else ["genotype.sh", "genotype"]

    _f(kwargs, aione, shell_fname, shell_log_folder, gatk_genotypeGVCFs)

    return aione


def variantrecalibrator(kwargs, aione):

    aione["intervals"] = []
    aione["genotype"] = {}  # will be called in ``gatk_variantrecalibrator``
    with open(kwargs.vcflist) as I:
        # Format in vcfilist: [Interval  vcf_file_path]
        for line in I:
            if line.startswith("#"):
                continue

            interval, vcf = line.strip().split()
            if interval not in aione["genotype"]:
                aione["intervals"].append(interval)

            aione["genotype"][interval] = vcf

    shell_fname, shell_log_folder = ["step6.VQSR.sh", "06.VQSR"] if kwargs.as_pipe_shell_order \
        else ["vqsr.sh", "genotype"]

    _f(kwargs, aione, shell_fname, shell_log_folder, gatk_variantrecalibrator)

    return aione
