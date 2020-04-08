"""Run analysis pipeline for NGS data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

"""
import argparse
import sys
import yaml

from datetime import datetime

from ilus.launch import runfunction
from ilus.launch.function_cmdline_args import bwamem_args, mergebam_args, markdup_args,\
    bqsr_args, haplotypecaller_args, genotypegvcfs_args, vqsr_args, mergevcfs_args

from ilus.tools.gatk import check_bundlefile_index


def parse_commandline_args(args):
    """Parse input commandline arguments, handling multiple cases.
    """
    functions ={"bwamem": bwamem_args,
                "mergebam": mergebam_args,
                "markdups": markdup_args,
                "BQSR": bqsr_args,
                "HaplotypeCaller": haplotypecaller_args,
                "GenotypeGVCFs": genotypegvcfs_args,
                "VQSR": vqsr_args,
                "MergeVcfs" : mergevcfs_args
                }

    desc = "Olympus: toolbox for NGS data analysis."
    parser = argparse.ArgumentParser(description=desc)
    subparser = parser.add_subparsers(help="Olympus supplemental commands")

    if len(args) > 0 and args[0] in functions: # the corresponding parser
        sp = functions[args[0]](subparser)#sub procedure parse

        # configuration for system
        sp.add_argument("-C", "--conf", dest="sysconf",
                        help="YAML configuration file specifying details about "
                             "system.")
    else:
        subparser.add_parser("bwamem", help=("Run bwamem aligning fastq to reference"))
        subparser.add_parser("mergebam", help=("Run merge bamfile by sambamba."))
        subparser.add_parser("markdups", help=("Run markduplicates by sambamba."))
        subparser.add_parser("BQSR", help=("Run BQSR by GATK4."))
        subparser.add_parser("HaplotypeCaller", help=("Run HaplotypeCaller by GATK4."))
        subparser.add_parser("GenotypeGVCFs", help=("Run GenotypeGVCFs by GATK4."))
        subparser.add_parser("VQSR", help=("Run VQSR by GATK4."))
        subparser.add_parser("MergeVcfs", help=("Run MergeVcfs by GATK4."))

    kwargs = {"args": parser.parse_args(args)}
    return kwargs


def checkconfig(config):
    """check GATK bundle is been indexed or not"""
    if config["resources"]["gatk_bundle"]:
        check_bundlefile_index(config["variantcaller"]["gatk"],
                              config["resources"]["gatk_bundle"])

    return


if __name__ == "__main__":
    START_TIME = datetime.now()
    kwargs = parse_commandline_args(sys.argv[1:])

    # all information in one.
    aione = {}
    # loaded global configuration file specifying details about program and
    # resource.

    with open(kwargs["args"].sysconf) as C:
        aione["config"] = yaml.load(C)

    # check bundle data is been index or not
    checkconfig(aione["config"])

    if "bwamem" in sys.argv[1:] and kwargs["args"]:
        runfunction.bwamem(kwargs, aione)

    if "mergebam" in sys.argv[1:] and kwargs:
        runfunction.mergebam(kwargs, aione)

    if "markdups" in sys.argv[1:] and kwargs["args"]:
        runfunction.markduplicates(kwargs, aione)

    if "BQSR" in sys.argv[1:] and kwargs["args"]:
        runfunction.baserecalibrator(kwargs, aione)

    if "HaplotypeCaller" in sys.argv[1:] and kwargs["args"]:
        runfunction.haplotypecaller(kwargs, aione)

    if "GenotypeGVCFs" in sys.argv[1:] and kwargs["args"]:
        runfunction.GenotypeGVCFs(kwargs, aione)

    if "VQSR" in sys.argv[1:] and kwargs["args"]:
        runfunction.variantrecalibrator(kwargs, aione)

    if "MergeVcfs" in sys.argv[1:] and kwargs["args"]:
        runfunction.MergeVCFs(kwargs, aione)

    elapsed_time = datetime.now() - START_TIME
    print "\n** %s done, %d seconds elapsed **\n" % (sys.argv[1], elapsed_time.seconds)

