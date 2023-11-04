"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

Author: Shujia Huang
Date: 2020-04-19

"""
import argparse
import sys
import yaml
from pathlib import Path
from datetime import datetime

# Import specific functions of ilus
from ilus.pipeline import (
    create_wgs_pipeline_command, WGS,
    create_genotype_joint_calling_command, genotypeGVCFs,
    create_vqsr_command, variantrecalibrator
)
from ilus.modules.utils import split_jobs, check_jobs_status

PROG_NAME = "ilus"
VERSION = "1.3.2"


def create_split_job_command(commands):
    # Create subparser for the "split-jobs" command
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
        help="The prefix name for output sub-shell. (default: %(default)s)"
    )
    split_job_cmd.add_argument(
        "-n", "--number",
        dest="number",
        type=int,
        required=True,
        help="Number of sub job."
    )
    split_job_cmd.add_argument(
        "-t", "--parallel",
        dest="t",
        type=int,
        required=True,
        help="Parallel number for per sub job."
    )

    return


def create_check_job_command(commands):
    # Create subparser for the "check-jobs" command
    check_job_cmd = commands.add_parser("check-jobs", help="Check the jobs have finished or not.")
    check_job_cmd.add_argument(
        "-I", "--input",
        dest="input",
        required=True,
        help="Task log file with suffix '.o.log.list’. For ilus pipeline, which could always "
             "be found in the folder 'loginfo', e.g.: loginfo/01.alignment.o.log.list"
    )

    return


def parse_commandline_args():
    """Parse input commandline arguments, handling multiple cases.
    """
    cmdparser = argparse.ArgumentParser(
        prog=PROG_NAME,
        description=f"{PROG_NAME} (Version = {VERSION}): A WGS/WES analysis pipeline generator.",
        epilog="That's how you can use %(prog)s"
    )

    cmdparser.add_argument(
        "-v", "--version",
        action="store_true",
        help=f"show the version of {PROG_NAME} and exit."
    )

    commands = cmdparser.add_subparsers(dest="command", title=f"{PROG_NAME} commands")

    # The arguments for the whole pipeline of WGS.
    create_wgs_pipeline_command(commands)

    # The arguments for joint-calling process
    create_genotype_joint_calling_command(commands)

    # The arguments for VQSR process
    create_vqsr_command(commands)

    # Utility tools
    create_split_job_command(commands)
    create_check_job_command(commands)

    return cmdparser.parse_args()


def load_config(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


def get_variant_calling_intervals(calling_interval_parameter):
    """get variant calling intervals into a bed format"""

    intervals = calling_interval_parameter
    if (type(calling_interval_parameter) is str) \
            and (Path(calling_interval_parameter).is_file()):
        # A file for recording interval
        interval_file = calling_interval_parameter
        with open(interval_file) as f:
            """Bed format:
            chr1	10001	207666
            chr1	257667	297968
            """
            # return the value to be a list of interval regions
            intervals = [line.strip().split()[:3] for line in f if not line.startswith("#")]

    elif type(calling_interval_parameter) is not list:
        raise ValueError(f"'variant_calling_interval' parameter could only be a file path or "
                         f"a list of chromosome id in the configure file (.yaml).\n")

    return intervals


def run_command(args):
    if args.version:
        print(f"{PROG_NAME} {VERSION}", file=sys.stderr)
        sys.exit(0)

    if args.command is None:
        print(f"Please type: {PROG_NAME} -h or {PROG_NAME} --help to show the help message.\n",
              file=sys.stderr)
        sys.exit(1)

    if args.command == "split-jobs":
        split_jobs(args.input, args.number, args.t, prefix=args.prefix)
        return

    if args.command == "check-jobs":
        check_jobs_status(args.input)
        return

    runner = {
        "WGS": WGS,
        "genotype-joint-calling": genotypeGVCFs,
        "VQSR": variantrecalibrator,
    }

    if args.command not in runner:
        raise ValueError(f"Invalid command: {args.command}")

    # loaded global configuration file
    config = load_config(args.sysconf)
    if not args.use_sentieon and "variant_calling_interval" in config["gatk"]:
        config["gatk"]["variant_calling_interval"] = get_variant_calling_intervals(
            config["gatk"]["variant_calling_interval"])

    elif args.use_sentieon and "variant_calling_interval" in config["sentieon"]:
        config["sentieon"]["variant_calling_interval"] = get_variant_calling_intervals(
            config["sentieon"]["variant_calling_interval"])
    else:
        raise ValueError(f"'variant_calling_interval' parameter is required "
                         f"in the configure file: {args.sysconf}.\n")

    # recording all information into one single dict.
    aione = {"config": config}
    runner[args.command](args, aione)

    return


def main():
    START_TIME = datetime.now()

    args = parse_commandline_args()
    run_command(args)

    elapsed_time = datetime.now() - START_TIME
    print(f"\n{PROG_NAME} (version: {VERSION}) for '{args.command}' done, "
          f"{elapsed_time.seconds} seconds elapsed.", file=sys.stderr)

    return
