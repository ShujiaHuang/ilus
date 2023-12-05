"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

Author: Shujia Huang
Date: 2020-04-19

"""
import argparse
import sys
import yaml
from datetime import datetime

# Import specific functions of ilus
from ilus.pipeline import (
    create_wgs_pipeline_command, create_capseq_pipeline_command, WGS,
    create_genotype_joint_calling_command, genotypeGVCFs,
    create_vqsr_command, variantrecalibrator
)
from ilus.modules.utils import split_jobs, check_jobs_status

PROG_NAME = "ilus"
VERSION = "2.0.0"


def create_split_job_command(commands):
    # Create sub-parser for the "split-jobs" command
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

    # The arguments for the whole pipeline of WES.
    create_capseq_pipeline_command(commands)

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


def run_command(args):
    """ Main function for pipeline.
    """
    if args.version:
        print(f"{PROG_NAME} {VERSION}", file=sys.stderr)
        sys.exit(1)

    if args.command is None:
        print(f"Please type: {PROG_NAME} -h or {PROG_NAME} --help "
              f"to show the help message.\n", file=sys.stderr)
        sys.exit(1)

    # Create WGS pipeline.
    elif args.command == "WGS":
        # loaded global configuration file and record all information into one single dict.
        aione = {"config": load_config(args.sysconf)}
        WGS(args, aione)

    elif args.command == "capseq":
        # Create pipeline for capture-sequencing. WES is one kind of capture sequencing.
        aione = {"config": load_config(args.sysconf)}
        WGS(args, aione, is_capture_seq=True)

    elif args.command == "genotype-joint-calling":
        aione = {"config": load_config(args.sysconf)}
        genotypeGVCFs(args, aione)

    elif args.command == "VQSR":
        aione = {"config": load_config(args.sysconf)}
        variantrecalibrator(args, aione)

    elif args.command == "split-jobs":
        split_jobs(args.input, args.number, args.t, prefix=args.prefix)

    elif args.command == "check-jobs":
        check_jobs_status(args.input)

    else:
        raise ValueError(f"Invalid command: {args.command}")

    return


def main():
    START_TIME = datetime.now()

    args = parse_commandline_args()
    run_command(args)

    elapsed_time = datetime.now() - START_TIME
    print(f"\n{PROG_NAME} (version: {VERSION}) for '{args.command}' done, "
          f"{elapsed_time.seconds} seconds elapsed.", file=sys.stderr)

    return
