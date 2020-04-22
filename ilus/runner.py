"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

"""
import argparse
import os
import stat
import sys

import yaml
from datetime import datetime

from ilus.utils import safe_makedir
from ilus.launch import runfunction
from ilus.tools.gatk import check_bundlefile_index


def parse_commandline_args(args):
    """Parse input commandline arguments, handling multiple cases.
    """
    desc = "Ilus: Toolbox for NGS data analysis."
    cmdparser = argparse.ArgumentParser(description=desc)

    command = cmdparser.add_subparsers(help="Ilus supplemental commands")
    pipeline_cmd = command.add_parser("pipeline", help="Creating pipeline for WGS(from fastq to VCF)")
    pipeline_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                              help="YAML configuration file specifying details about system.")
    pipeline_cmd.add_argument("-L", "--fastqlist", dest="fastqlist", type=str, required=True,
                              help="Alignment FASTQ Index File.")
    pipeline_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                              help="Name of the project. [test]")
    pipeline_cmd.add_argument("-f", "--force", dest="overwrite", type=bool, default=False,
                              help="Force overwrite existing shell scripts and folders.")
    pipeline_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                              help="A directory for output results.")

    kwargs = {"args": cmdparser.parse_args(args)}

    return kwargs


def checkconfig(config):
    """check GATK bundle is been indexed or not"""
    if config["resources"]["gatk_bundle"]:
        check_bundlefile_index(config["variantcaller"]["gatk"], config["resources"]["gatk_bundle"])

    return


def create_a_total_shell_file(shell_list, out_shell_filename, sub_shell_log_dir, o_log_file, e_log_file):
    """Creat all the executable shell into a big single shell files.
    ``shell_list`` is a 2-D array: [[mark, shell_file], ...].
    """
    with open(out_shell_filename, "w") as OUT, open(o_log_file, "w") as O_LOG, open(e_log_file, "w") as E_LOG:
        OUT.write("#!/bin/bash\n")
        for marker, sub_shell in shell_list:
            OUT.write("{sub_shell} 2> {sub_shell_log_dir}/{marker}.e.log > "
                      "{sub_shell_log_dir}/{marker}.o.log\n".format(**locals()))

            # record all the path of log files into a single file
            O_LOG.write("{sub_shell_log_dir}/{marker}.o.log\n".format(**locals()))
            E_LOG.write("{sub_shell_log_dir}/{marker}.e.log\n".format(**locals()))

    os.chmod(out_shell_filename, stat.S_IRWXU)  # 0700
    return


def make_process_shell(output_shell_fname, shell_log_directory, process_shells=None, is_overwrite=False):
    safe_makedir(shell_log_directory)
    o_log_file = shell_log_directory + ".o.log.list"
    e_log_file = shell_log_directory + ".e.log.list"

    if not is_overwrite and os.path.exists(output_shell_fname):
        return

    create_a_total_shell_file(process_shells, output_shell_fname, shell_log_directory, o_log_file, e_log_file)
    return


if __name__ == "__main__":
    START_TIME = datetime.now()
    kwargs = parse_commandline_args(sys.argv[1:])

    runner_module = {  # [func, shell_file, shell_log_folder, output_folder]
        # bwa/sort/merge process
        "align": [runfunction.bwamem, "step1.bwa.sh", "01.alignment", "01.alignment"],

        # Create Markduplicates running shells.
        "markdup": [runfunction.gatk_markduplicates, "step2.markdup.sh", "02.markdup", "01.alignment"],

        # Create BQSR+ApplyBQSR running shells.
        "BQSR": [runfunction.gatk_baserecalibrator, "step3.bqsr.sh", "03.BQSR", "01.alignment"],

        # Create GVCF running shells
        "gvcf": [runfunction.gatk_haplotypecaller_gvcf, "step4.gvcf.sh", "04.gvcf", "02.gvcf"],

        # GenotypeGVCF
        "genotype": [runfunction.gatk_genotypeGVCFs, "step5.genotype.sh", "05.genotype", "03.genotype"],

        # Summary and status statistic
        "summary": []
    }
    process = ["align", "markdup", "BQSR", "gvcf", "genotype"]

    # all information in one dict.
    aione = {}

    # loaded global configuration file
    with open(kwargs["args"].sysconf) as C:
        aione["config"] = yaml.safe_load(C)

    # check bundle data is been index or not
    # checkconfig(aione["config"])

    # Create project directory
    kwargs["args"].outdir = safe_makedir(os.path.abspath(kwargs["args"].outdir))  # return abspath
    shell_dirtory = os.path.join(kwargs["args"].outdir, "00.shell")
    shell_log_dirtory = os.path.join(shell_dirtory, "loginfo")
    safe_makedir(shell_dirtory)
    safe_makedir(shell_log_dirtory)

    for p in process:
        func, shell_fname, shell_log_folder, output_result_folder = runner_module[p]
        make_process_shell(output_shell_fname=os.path.join(shell_dirtory, shell_fname),
                           shell_log_directory=os.path.join(shell_log_dirtory, shell_log_folder),
                           process_shells=func(kwargs, output_result_folder, aione))

    elapsed_time = datetime.now() - START_TIME
    print ("\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds))
