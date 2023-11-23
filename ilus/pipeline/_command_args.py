"""Arguments for pipelines

Author: Shujia Huang
Date:   2023-02-16
"""
import argparse


def _parent_parser():
    """Serves as a parent parser to record all the common arguments for ilus.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-n", "--name",
        dest="project_name",
        type=str,
        default="test",
        help="Name of the project. (Default: %(default)s)"
    )

    parser.add_argument(
        "-C", "--conf",
        dest="sysconf",
        type=str,
        required=True,
        help="YAML configuration file specifying system details."
    )

    parser.add_argument(
        "-O", "--outdir",
        dest="outdir",
        type=str,
        required=True,
        help="Output directory for results."
    )

    parser.add_argument(
        "-f", "--force-overwrite",
        dest="overwrite",
        action="store_true",
        help="Force overwrite existing shell scripts and folders."
    )

    parser.add_argument(
        "--use-sentieon",
        dest="use_sentieon",
        action="store_true",
        help="Use sentieon (doc: https://support.sentieon.com/manual) "
             "to create analysis pipeline."
    )

    return parser


def _add_germline_short_variant_discovery_argument(command):
    """Add argument for variant discovery.
    """
    command.add_argument(
        "-I", "--fastqlist",
        dest="fastqlist",
        type=str,
        required=True,
        help="Input the list of FASTQ files."
    )

    command.add_argument(
        "--clean-raw-data",
        dest="clean_raw_data",
        action="store_true",
        help="Set this option to clean raw sequencing data (fastq)."
    )

    command.add_argument(
        "--delete-clean-fastq",
        dest="delete_clean_fastq",
        action="store_true",
        help="Set this option to delete clean fastq after aligning all reads to reference."
    )

    command.add_argument(
        "-c", "--cram",
        dest="cram",
        action="store_true",
        help="Convert BAM to CRAM after BQSR and save alignment file storage."
    )

    command.add_argument(
        "-P", "--process",
        dest="wgs_processes",
        type=str,
        default="align,markdup,BQSR,gvcf,combineGVCFs,genotype,VQSR",
        help="Specify one or more processes (separated by comma) of WGS pipeline. "
             "Possible values: %(default)s"
    )

    # Todo: Write a dry run function for testing the pipeline without truely run the pipeline.
    command.add_argument(
        "-dr", "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Dry run the pipeline for testing."
    )

    return command


def create_wgs_pipeline_command(commands):
    """Add arguments to create WGS pipeline command."""
    # Generate the commandline argument for WGS pipeline.
    wgs_cmd = _add_germline_short_variant_discovery_argument(
        commands.add_parser(
            "WGS",
            parents=[_parent_parser()],
            help="Create pipeline scripts for WGS (from FASTQ to genotype VCF).")
    )

    # Specific variant calling intervals.
    # The value could be a file in bed format (I show you a example bellow) or a interval of list.
    # Bed format of interval file only contain three columns: ``Sequencing ID``, ``region start``
    # and ``region end``, e.g.:
    #   chr1    10001   207666
    #   chr1    257667  297968
    #
    # These invertals could be any regions alone the genome as you wish or just set the same as
    # ``--interval`` argument above.
    wgs_cmd.add_argument(
        "--interval",
        dest="interval",
        type=str,
        required=False,
        help="Interval strings (separate by comma) or a file (BED/Picard,1-based) that "
             "will be used in variants calling. "
             "e.g: '--interval chr1:1-2,chr2,chr3:4-5' "
             "or '--interval interval_file.bed'."
    )

    return wgs_cmd


def create_wes_pipeline_command(commands):
    """All the arguments for creating WES pipeline."""
    wes_cmd = _add_germline_short_variant_discovery_argument(
        commands.add_parser(
            "WES",
            parents=[_parent_parser()],
            help="Create pipeline scripts for WES (from FASTQ to genotype VCF).")
    )

    wes_cmd.add_argument(
        "--capture-interval",
        dest="interval",  # WES capture region file (BED/Picard,1-based)
        type=str,
        required=True,
        help="WES capture intervals: string or file (BED/Picard)"
    )

    return wes_cmd


def create_genotype_joint_calling_command(commands):
    """Add arguments to create genotype joint calling command."""
    genotype_cmd = commands.add_parser(
        "genotype-joint-calling",
        parents=[_parent_parser()],
        help="Genotype from GVCFs."
    )

    genotype_cmd.add_argument(
        "-I", "--gvcflist",
        dest="gvcflist",
        type=str,
        required=True,
        help="List of GVCF files. One gvcf_file per-row and the format should looks like: "
             "[interval\tgvcf_file_path]. Column [1] is a symbol which could represent "
             "the genome region of the gvcf_file and column [2] should be the path."
    )

    genotype_cmd.add_argument(
        "--as_pipe_shell_order",
        dest="as_pipe_shell_order",
        action="store_true",
        help="Keep the shell name as the order of `WGS`."
    )

    return genotype_cmd


def create_vqsr_command(commands):
    """Add arguments to create VQSR command."""
    vqsr_cmd = commands.add_parser("VQSR",
                                   parents=[_parent_parser()],
                                   help="VQSR")

    vqsr_cmd.add_argument(
        "-I", "--vcflist",
        dest="vcflist",
        type=str,
        required=True,
        help="VCFs file list. One file per-row."
    )

    vqsr_cmd.add_argument(
        "--as_pipe_shell_order",
        dest="as_pipe_shell_order",
        action="store_true",
        help="Keep the shell name as the order of `WGS`."
    )

    return vqsr_cmd
