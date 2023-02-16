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


def load_config(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


def get_intervals(interval_file):
    if not os.path.isfile(interval_file):
        raise ValueError(f"Invalid interval file: {interval_file}")

    with open(interval_file) as f:
        """Bed format:
        chr1	10001	207666
        chr1	257667	297968
        """
        return [line.strip().split()[:3] for line in f if not line.startswith("#")]


def run_command(args):
    if args.version:
        print("ilus " + VERSION, file=sys.stderr)
        sys.exit(0)

    if args.command is None:
        print("Please type: ilus -h or ilus --help to show the help message.\n", file=sys.stderr)
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

    if "variant_calling_interval" in config["gatk"]:
        if type(config["gatk"]["variant_calling_interval"]) is str:
            # A file for recording interval
            interval_file = config["gatk"]["variant_calling_interval"]
            config["gatk"]["variant_calling_interval"] = get_intervals(interval_file)

        elif type(config["gatk"]["variant_calling_interval"]) is not list:
            raise ValueError(f"'variant_calling_interval' parameter could only be a file path or a "
                             f"list of chromosome id in the configure file: {args.sysconf}.\n")

    else:
        raise ValueError(f"'variant_calling_interval' parameter is required "
                         f"in the configure file: {args.sysconf}.\n")

    # recording all information into one single dict.
    aione = {"config": config}
    runner[args.command](args, aione)

    return


def main():
    START_TIME = datetime.now()

    run_command(parse_commandline_args())

    elapsed_time = datetime.now() - START_TIME
    print(f"\n** {sys.argv[1]} done, {elapsed_time.seconds} seconds elapsed. **\n")

    return
