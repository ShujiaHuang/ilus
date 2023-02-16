"""Arguments for pipelines

Author: Shujia Huang
Date:   2023-02-16
"""


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

