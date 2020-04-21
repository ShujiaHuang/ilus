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


def make_process_shell(output_shell_fname, shell_log_directory, process_shells=None):
    safe_makedir(shell_log_directory)
    o_log_file = shell_log_directory + ".o.log.list"
    e_log_file = shell_log_directory + ".e.log.list"
    create_a_total_shell_file(process_shells, output_shell_fname, shell_log_directory, o_log_file, e_log_file)
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

    # Create project directory
    kwargs["args"].outdir = safe_makedir(os.path.abspath(kwargs["args"].outdir))  # return abspath

    shell_dirtory = os.path.join(kwargs["args"].outdir, "00.shell")
    shell_log_dirtory = os.path.join(shell_dirtory, "loginfo")
    safe_makedir(shell_dirtory)
    safe_makedir(shell_log_dirtory)

    # bwa/sort/merge process
    make_process_shell(output_shell_fname=os.path.join(shell_dirtory, "step1.bwa.sh"),
                       shell_log_directory=os.path.join(shell_log_dirtory, "01.alignment"),
                       process_shells=runfunction.bwamem(kwargs, "01.alignment", aione))

    # Create Markdupkicates running shells.
    make_process_shell(output_shell_fname=os.path.join(shell_dirtory, "step2.markdup.sh"),
                       shell_log_directory=os.path.join(shell_log_dirtory, "02.markdup"),
                       process_shells=runfunction.gatk_markduplicates(kwargs, "01.alignment", aione))

    # Create BQSR+ApplyBQSR running shells.
    make_process_shell(output_shell_fname=os.path.join(shell_dirtory, "step3.bqsr.sh"),
                       shell_log_directory=os.path.join(shell_log_dirtory, "03.BQSR"),
                       process_shells=runfunction.gatk_baserecalibrator(kwargs, "01.alignment", aione))

    # Create GVCF running shells
    make_process_shell(output_shell_fname=os.path.join(shell_dirtory, "step4.gvcf.sh"),
                       shell_log_directory=os.path.join(shell_log_dirtory, "04.gvcf"),
                       process_shells=runfunction.gatk_haplotypecaller_gvcf(kwargs, "02.gvcf", aione))

    # GenotypeGVCF
    make_process_shell(output_shell_fname=os.path.join(shell_dirtory, "step5.genotype.sh"),
                       shell_log_directory=os.path.join(shell_log_dirtory, "05.genotype"),
                       process_shells=runfunction.gatk_genotypeGVCFs(kwargs, "03.genotype", aione))

    # Summary and status statistic
    elapsed_time = datetime.now() - START_TIME
    print ("\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds))
