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
    bwamem_cmd = command.add_parser("bwamem", help="Run bwamem aligning fastq to reference")
    bwamem_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                            help="YAML configuration file specifying details about system.")
    bwamem_cmd.add_argument("-L", "--fastqlist", dest="fastqlist", type=str, required=True,
                            help="Alignment FASTQ Index File.")
    bwamem_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                            help="A directory for output results.")

    kwargs = {"args": cmdparser.parse_args(args)}

    return kwargs


def checkconfig(config):
    """check GATK bundle is been indexed or not"""
    if config["resources"]["gatk_bundle"]:
        check_bundlefile_index(config["variantcaller"]["gatk"], config["resources"]["gatk_bundle"])

    return


def create_a_total_shell_file(sample_shell_files, out_shell_filename, sub_shell_log_dirtory):
    """Creat all the executable shell into a big single shell files.
    ``sample_shell_files`` is a 2-D array: [[sample, sample_shell_file], ...].
    """
    with open(out_shell_filename, "w") as OUT:
        OUT.write("#!/bin/bash\n")
        for sample, sub_shell in sample_shell_files:
            OUT.write("{sub_shell} 2> {sub_shell_log_dirtory}/{sample}.e.log > "
                      "{sub_shell_log_dirtory}/{sample}.o.log\n".format(**locals()))

    os.chmod(out_shell_filename, stat.S_IRWXU)  # 0700
    return


if __name__ == "__main__":
    START_TIME = datetime.now()
    kwargs = parse_commandline_args(sys.argv[1:])

    # all information in one dict.
    aione = {}

    # loaded global configuration file specifying details about program and
    # resource.
    with open(kwargs["args"].sysconf) as C:
        aione["config"] = yaml.safe_load(C)

    # check bundle data is been index or not
    # checkconfig(aione["config"])

    # Create project diretory
    kwargs["args"].outdir = safe_makedir(os.path.abspath(kwargs["args"].outdir))
    shell_dirtory = os.path.join(kwargs["args"].outdir, "00.shell")
    shell_log_dirtory = os.path.join(kwargs["args"].outdir, "00.shell", "loginfo")
    safe_makedir(shell_dirtory)
    safe_makedir(shell_log_dirtory)

    # bwa/sort/merge process
    bwa_shell_log_dirtory = os.path.join(shell_log_dirtory, "01.alignment")
    safe_makedir(bwa_shell_log_dirtory)
    bwa_shell_files_list = runfunction.bwamem(kwargs, "01.alignment", aione)
    create_a_total_shell_file(bwa_shell_files_list,
                              os.path.join(shell_dirtory, "all_samples_bwa.sh"),
                              bwa_shell_log_dirtory)

    # Markdupkicates
    markdup_shell_log_dirtory = os.path.join(shell_log_dirtory, "02.markdup")
    safe_makedir(markdup_shell_log_dirtory)

    # BQSR+ApplyBQSR

    # gvcf

    elapsed_time = datetime.now() - START_TIME
    print ("\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds))
