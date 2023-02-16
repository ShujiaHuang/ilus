"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

Author: Shujia Huang
Date: 2020-04-19

"""
import argparse
import os
import sys
import yaml
from datetime import datetime

from ilus.pipeline.wgs import wgs, genotypeGVCFs, variantrecalibrator
from ilus.utils import split_jobs, check_jobs_status

VERSION = "1.3.0"


def create_wgs_pipeline_command(commands):
    """Add arguments to create WGS pipeline command."""
    wgs_cmd = commands.add_parser(
        "WGS",
        help="Create a pipeline for WGS (from FASTQ to genotype VCF)."
    )

    wgs_cmd.add_argument(
        "-n", "--name",
        dest="project_name",
        type=str,
        default="test",
        help="Name of the project. Default value: test"
    )

    wgs_cmd.add_argument(
        "-C", "--conf",
        dest="sysconf",
        required=True,
        help="YAML configuration file specifying system details."
    )

    wgs_cmd.add_argument(
        "-L", "--fastqlist",
        dest="fastqlist",
        type=str,
        required=True,
        help="List of alignment FASTQ files."
    )

    wgs_cmd.add_argument(
        "-O", "--outdir",
        dest="outdir",
        required=True,
        help="Output directory for results."
    )

    wgs_cmd.add_argument(
        "-c", "--cram",
        dest="cram",
        action="store_true",
        help="Convert BAM to CRAM after BQSR and save alignment file storage."
    )

    wgs_cmd.add_argument(
        "-P", "--process",
        dest="wgs_processes",
        type=str,
        default="align,markdup,BQSR,gvcf,combineGVCFs,genotype,VQSR",
        help="Specify one or more processes (separated by comma) of WGS pipeline. "
             "Possible values: {align,markdup,BQSR,gvcf,combineGVCFs,genotype,VQSR}"
    )

    wgs_cmd.add_argument(
        "-f", "--force-overwrite",
        dest="overwrite",
        action="store_true",
        help="Force overwrite existing shell scripts and folders."
    )

    # Todo: Write a dry run function for testing the pipeline without truely run the pipeline.
    wgs_cmd.add_argument(
        "-dr", "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Dry run the pipeline for testing."
    )

    return wgs_cmd


def create_genotype_joint_calling_command(commands):
    """Add arguments to create genotype joint calling command."""
    genotype_cmd = commands.add_parser(
        "genotype-joint-calling",
        help="Genotype from GVCFs."
    )

    genotype_cmd.add_argument(
        "-C", "--conf",
        dest="sysconf",
        required=True,
        help="YAML configuration file specifying system details."
    )

    genotype_cmd.add_argument(
        "-L", "--gvcflist",
        dest="gvcflist",
        type=str,
        required=True,
        help="List of GVCF files. One gvcf_file per-row and the format should looks like: "
             "[interval\tgvcf_file_path]. Column [1] is a symbol which could represent "
             "the genome region of the gvcf_file and column [2] should be the path."
    )

    genotype_cmd.add_argument(
        "-O", "--outdir",
        dest="outdir",
        required=True,
        help="A directory for output results."
    )

    genotype_cmd.add_argument(
        "-n", "--name",
        dest="project_name",
        type=str,
        default="test",
        help="Name of the project. [test]"
    )

    genotype_cmd.add_argument(
        "--as_pipe_shell_order",
        dest="as_pipe_shell_order",
        action="store_true",
        help="Keep the shell name as the order of `WGS`."
    )

    genotype_cmd.add_argument(
        "-f", "--force",
        dest="overwrite",
        action="store_true",
        help="Force overwrite existing shell scripts and folders."
    )

    return genotype_cmd


def create_vqsr_command(commands):
    """Add arguments to create VQSR command."""
    vqsr_cmd = commands.add_parser("VQSR", help="VQSR")
    vqsr_cmd.add_argument(
        "-C", "--conf",
        dest="sysconf",
        required=True,
        help="YAML configuration file specifying details about system."
    )

    vqsr_cmd.add_argument(
        "-L", "--vcflist",
        dest="vcflist",
        type=str,
        required=True,
        help="VCFs file list. One file per-row."
    )

    vqsr_cmd.add_argument(
        "-O", "--outdir",
        dest="outdir",
        required=True,
        help="A directory for output results."
    )

    vqsr_cmd.add_argument(
        "-n", "--name",
        dest="project_name",
        type=str,
        default="test",
        help="Name of the project. [test]"
    )

    vqsr_cmd.add_argument(
        "--as_pipe_shell_order",
        dest="as_pipe_shell_order",
        action="store_true",
        help="Keep the shell name as the order of `WGS`."
    )

    vqsr_cmd.add_argument(
        "-f", "--force",
        dest="overwrite",
        action="store_true",
        help="Force overwrite existing shell scripts and folders."
    )

    return vqsr_cmd


def create_utility_module_command(commands):
    split_job_cmd = commands.add_parser("split-jobs", help="Split the whole shell into multiple jobs.")

    split_job_cmd.add_argument(
        "-I", "--input",
        dest="input",
        required=True,
        help="Input shell file."
    )

    split_job_cmd.add_argument(
        "-p", "--prefix",
        dest="prefix",
        type=str,
        default="work",
        help="The prefix name of output sub-shell file. [work]"
    )

    split_job_cmd.add_argument(
        "-n", "--number",
        dest="number",
        type=int,
        required=True,
        help="The number of sub tasks (shells)."
    )

    split_job_cmd.add_argument(
        "-t", "--parallel",
        dest="t",
        type=int,
        required=True,
        help="The number of parallel jobs."
    )

    check_job_cmd = commands.add_parser("check-jobs", help="Check the jobs have finished or not.")
    check_job_cmd.add_argument(
        "-I", "--input",
        dest="input",
        required=True,
        help="Input the log file of task."
    )

    return


def parse_commandline_args():
    """Parse input commandline arguments, handling multiple cases.
    """
    desc = f"ilus (Version = {VERSION}): A WGS/WES analysis pipeline generator."
    cmdparser = argparse.ArgumentParser(description=desc)
    cmdparser.add_argument(
        "-v", "--version",
        action="store_true",
        help="show the version of ilus and exit."
    )

    commands = cmdparser.add_subparsers(dest="command", title="ilus commands")

    # The arguments for WGS.
    create_wgs_pipeline_command(commands)

    # The arguments for GVCFs
    create_genotype_joint_calling_command(commands)

    # The arguments for VQSR
    create_vqsr_command(commands)

    # Utility tools
    create_utility_module_command(commands)

    return cmdparser.parse_args()


def main():
    START_TIME = datetime.now()
    runner = {
        "WGS": wgs,
        "genotype-joint-calling": genotypeGVCFs,
        "VQSR": variantrecalibrator
    }

    kwargs = parse_commandline_args()
    if kwargs.version:
        print("ilus " + VERSION, file=sys.stderr)
        sys.exit(0)

    if kwargs.command is None:
        sys.stderr.write("Please type: ilus -h or ilus --help to show the help message.\n")
        sys.exit(1)

    elif kwargs.command in runner:
        aione = {}  # All information in one dict.
        with open(kwargs.sysconf) as C:
            # loaded global configuration file
            aione["config"] = yaml.safe_load(C)

        # Normalize interval regions
        if "variant_calling_interval" in aione["config"]["gatk"]:
            intervals = []  # A 2-D array
            if os.path.isfile(aione["config"]["gatk"]["variant_calling_interval"][0]):
                with open(aione["config"]["gatk"]["variant_calling_interval"][0]) as I:
                    """ Bed format:
                    chr1	10001	207666
                    chr1	257667	297968
                    """
                    for line in I:
                        if line.startswith("#"):
                            continue
                        else:
                            intervals.append(line.strip().split()[:3])

                # Update by regular regions information
                aione["config"]["gatk"]["variant_calling_interval"] = intervals
        else:
            sys.stderr.write("[Error] 'variant_calling_interval' parameter in configure "
                             "file %s is required.\n" % kwargs.sysconf)
            sys.exit(1)

        runner[kwargs.command](kwargs, aione)
    elif kwargs.command == "split-jobs":
        split_jobs(kwargs.input, kwargs.number, kwargs.t, prefix=kwargs.prefix)

    elif kwargs.command == "check-jobs":
        check_jobs_status(kwargs.input)

    else:
        raise ValueError("[ERROR] Invalid command: '%s'." % kwargs.command)

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** %s done, %d seconds elapsed **\n" %
                     (sys.argv[1], elapsed_time.seconds))
