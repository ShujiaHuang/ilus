"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

"""
import argparse
import sys
import yaml

from datetime import datetime

from ilus.pipeline import wgs, genotypeGVCFs, variantrecalibrator


# from ilus.tools.gatk import check_bundlefile_index
# def checkconfig(config):
#     """check GATK bundle is been indexed or not"""
#     if config["resources"]["gatk_bundle"]:
#         check_bundlefile_index(config["variantcaller"]["gatk"], config["resources"]["gatk_bundle"])
#
#     return


def parse_commandline_args():
    """Parse input commandline arguments, handling multiple cases.
    """
    desc = "ilus: A WGS analysis pipeline."
    cmdparser = argparse.ArgumentParser(description=desc)
    commands = cmdparser.add_subparsers(dest="command", title="ilus commands")

    # The standard pipeline for WGS.
    pipeline_cmd = commands.add_parser("WGS", help="Creating pipeline for WGS(from fastq to genotype VCF)")
    pipeline_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                              help="YAML configuration file specifying details about system.")
    pipeline_cmd.add_argument("-L", "--fastqlist", dest="fastqlist", type=str, required=True,
                              help="Alignment FASTQ Index File.")
    pipeline_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                              help="Name of the project. [test]")
    pipeline_cmd.add_argument("-f", "--force", dest="overwrite", action="store_true",
                              help="Force overwrite existing shell scripts and folders.")
    pipeline_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                              help="A directory for output results.")

    # Genotype from GVCFs
    genotype_cmd = commands.add_parser("genotype-joint-calling", help="Genotype from GVCFs.")
    genotype_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                              help="YAML configuration file specifying details about system.")
    genotype_cmd.add_argument("-L", "--gvcflist", dest="gvcflist", type=str, required=True,
                              help="GVCFs file list. One gvcf_file per-row and the format should looks like: "
                                   "[interval\tgvcf_file_path]. Column [1] is a symbol which could represent "
                                   "the genome region of the gvcf_file and column [2] should be the path.")
    genotype_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                              help="Name of the project. [test]")
    genotype_cmd.add_argument("--as_pipe_shell_order", dest="as_pipe_shell_order", action="store_true",
                              help="Keep the shell name as the order of `pipeline`.")
    genotype_cmd.add_argument("-f", "--force", dest="overwrite", action="store_true",
                              help="Force overwrite existing shell scripts and folders.")
    genotype_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                              help="A directory for output results.")

    # Genotype from GVCFs
    vqsr_cmd = commands.add_parser("VQSR", help="VQSR")
    vqsr_cmd.add_argument("-C", "--conf", dest="sysconf", required=True,
                          help="YAML configuration file specifying details about system.")
    vqsr_cmd.add_argument("-L", "--vcflist", dest="vcflist", type=str, required=True,
                          help="VCFs file list. One vcf_file per-row and the format should looks like: "
                               "[interval\tvcf_file_path]. Column [1] is a symbol which could represent "
                               "the genome region of the vcf_file and column [2] should be the path.")
    vqsr_cmd.add_argument("-n", "--name", dest="project_name", type=str, default="test",
                          help="Name of the project. [test]")
    vqsr_cmd.add_argument("--as_pipe_shell_order", dest="as_pipe_shell_order", action="store_true",
                          help="Keep the shell name as the order of `pipeline`.")
    vqsr_cmd.add_argument("-f", "--force", dest="overwrite", action="store_true",
                          help="Force overwrite existing shell scripts and folders.")
    vqsr_cmd.add_argument("-O", "--outdir", dest="outdir", required=True,
                          help="A directory for output results.")

    return cmdparser.parse_args()


if __name__ == "__main__":
    START_TIME = datetime.now()
    runner = {
        "WGS": wgs,
        "genotype-joint-calling": genotypeGVCFs,
        "VQSR": variantrecalibrator
    }

    kwargs = parse_commandline_args()

    # all information in one dict.
    aione = {}

    # loaded global configuration file
    with open(kwargs.sysconf) as C:
        aione["config"] = yaml.safe_load(C)

    runner[kwargs.command](kwargs, aione)
    elapsed_time = datetime.now() - START_TIME
    print ("\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds))
