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

from ilus.pipeline import WGS, genotypeGVCFs, variantrecalibrator, \
    create_wgs_pipeline_command, create_genotype_joint_calling_command, \
    create_vqsr_command
from ilus.utils import split_jobs, check_jobs_status

PROG_NAME = "ilus"
VERSION = "1.3.0"


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
    cmdparser = argparse.ArgumentParser(
        prog=PROG_NAME,
        description=f"{PROG_NAME} (Version = {VERSION}): A WGS/WES analysis pipeline generator.",
        epilog="That's how you could use %(prog)s"
    )

    cmdparser.add_argument(
        "-v", "--version",
        action="store_true",
        help=f"show the version of {PROG_NAME} and exit."
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
        "WGS": WGS,
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
