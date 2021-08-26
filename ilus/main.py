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


def parse_commandline_args():
    """Parse input commandline arguments, handling multiple cases.
    """
    desc = "ilus: A WGS/WES analysis pipeline generator."
    cmdparser = argparse.ArgumentParser(description=desc)
    commands = cmdparser.add_subparsers(dest="command", title="ilus commands")

    # The standard pipeline for WGS.
    pipeline_cmd = commands.add_parser("WGS", help="Creating pipeline for WGS(from fastq to genotype VCF)")
    pipeline_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                              help="YAML configuration file specifying details about system.")
    pipeline_cmd.add_argument("-L", "--fastqlist", dest="fastqlist", type=str, required=True,
                              help="Alignment FASTQ Index File.")
    pipeline_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                              help="A directory for output results.")

    pipeline_cmd.add_argument("-P", "--Process", dest="wgs_processes", type=str,
                              help="Specific one or more processes (separated by comma) of WGS pipeline. "
                                   "Defualt value: align,markdup,BQSR,gvcf,genotype,VQSR. "
                                   "Possible values: {align,markdup,BQSR,gvcf,genotype,VQSR}",
                              default="align,markdup,BQSR,gvcf,genotype,VQSR")

    pipeline_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                              help="Name of the project. Default value: test")
    pipeline_cmd.add_argument("-f", "--force_overwrite", dest="overwrite", action="store_true",
                              help="Force overwrite existing shell scripts and folders.")
    pipeline_cmd.add_argument("-c", "--cram", dest="cram", action="store_true",
                              help="Covert BAM to CRAM after BQSR and save alignment file storage.")

    # Genotype from GVCFs
    genotype_cmd = commands.add_parser("genotype-joint-calling", help="Genotype from GVCFs.")
    genotype_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                              help="YAML configuration file specifying details about system.")
    genotype_cmd.add_argument("-L", "--gvcflist", dest="gvcflist", type=str, required=True,
                              help="GVCFs file list. One gvcf_file per-row and the format should looks like: "
                                   "[interval\tgvcf_file_path]. Column [1] is a symbol which could represent "
                                   "the genome region of the gvcf_file and column [2] should be the path.")
    genotype_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                              help="A directory for output results.")

    genotype_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                              help="Name of the project. [test]")
    genotype_cmd.add_argument("--as_pipe_shell_order", dest="as_pipe_shell_order", action="store_true",
                              help="Keep the shell name as the order of `WGS`.")
    genotype_cmd.add_argument("-f", "--force", dest="overwrite", action="store_true",
                              help="Force overwrite existing shell scripts and folders.")

    # Genotype from VQSR
    vqsr_cmd = commands.add_parser("VQSR", help="VQSR")
    vqsr_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                          help="YAML configuration file specifying details about system.")
    vqsr_cmd.add_argument("-L", "--vcflist", dest="vcflist", type=str, required=True,
                          help="VCFs file list. One file per-row.")
    vqsr_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                          help="A directory for output results.")

    vqsr_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                          help="Name of the project. [test]")
    vqsr_cmd.add_argument("--as_pipe_shell_order", dest="as_pipe_shell_order", action="store_true",
                          help="Keep the shell name as the order of `WGS`.")
    vqsr_cmd.add_argument("-f", "--force", dest="overwrite", action="store_true",
                          help="Force overwrite existing shell scripts and folders.")

    return cmdparser.parse_args()


def main():
    START_TIME = datetime.now()
    runner = {
        "WGS": wgs,
        "genotype-joint-calling": genotypeGVCFs,
        "VQSR": variantrecalibrator
    }

    kwargs = parse_commandline_args()
    if kwargs.command is None:
        print("Please type: ilus -h or ilus --help to show the help message.")
        sys.exit(1)

    # all information in one dict.
    aione = {}

    # loaded global configuration file
    with open(kwargs.sysconf) as C:
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
        sys.stderr.write("[Error] 'variant_calling_interval' parameter in [\"gatk\"][\"variant_calling_interval\"] "
                         "in configure file %s is required." % kwargs.sysconf)
        sys.exit(1)

    runner[kwargs.command](kwargs, aione)
    elapsed_time = datetime.now() - START_TIME
    print("\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds))
