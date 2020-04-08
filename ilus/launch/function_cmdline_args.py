"""Provide commandline arguments for subfunctions.
"""

def bwamem_args(cdl_parser):
    # commandline arguments for bwamem
    bwamem_parser = cdl_parser.add_parser("bwamem", help=("Run bwamem aligning fastq to reference"))

    # fastq metadata
    bwamem_parser.add_argument("-R", "--readgroup", dest="rg_info",
                               help="read group header line such as '@RG\\tID:foo\\tSM:bar'")
    bwamem_parser.add_argument("-L", "--lane", dest="lane", help="Sequencing lane ID.")

    bwamem_parser.add_argument("-O", "--outdir", dest="outdir", help="A directory for output results.")
    bwamem_parser.add_argument("fastq1", help="input fastq file of read1 for SE or PE.")
    bwamem_parser.add_argument("fastq2", default="",
                               help="input fastq file of read2 for PE (optional). ")

    return bwamem_parser


def mergebam_args(cdl_parser):
    """commandline arguments for markduplicates by ``sambamba markdup``"""
    mergebam_parser = cdl_parser.add_parser("mergebam", help=("Run merge bamfile by sambamba."))

    mergebam_parser.add_argument("-I", "--input", dest="inbam", action='append', default=[],
                                help="Input one or more BAM files to analyze. This "
                                     "argument must be specified at least once. Reuqired")
    mergebam_parser.add_argument("-O", "--outfile", dest="outfile",
                                help="A BAM file's name for output results.")

    return mergebam_parser


def markdup_args(cdl_parser):
    """commandline arguments for markduplicates by ``sambamba markdup``"""
    markdup_parser = cdl_parser.add_parser("markdups", help=("Run markduplicates by sambamba."))

    markdup_parser.add_argument("-I", "--input", dest="inbam", action='append', default=[],
                                help="Input one or more BAM files to analyze. This "
                                     "argument must be specified at least once. Reuqired")
    markdup_parser.add_argument("--REMOVE_DUPLICATES", dest="removedups", type=bool, default=False,
                                help="If true do not write duplicates to the output file "
                                     "instead of writing them with appropriate flags set. [False]")
    markdup_parser.add_argument("-O", "--outfile", dest="outfile",
                                help="A BAM file's name for output results.")

    return markdup_parser


def gatk_markdup_args(cdl_parser):
    """commandline arguments for markduplicates by GATK4 Markduplicates.
    """
    markdup_parser = cdl_parser.add_parser("markdups", help=("Run markduplicates by GATK4."))

    markdup_parser.add_argument("-I", "--input", dest="inbam", action='append', default=[],
                                help="Input one or more BAM files to analyze. This "
                                     "argument must be specified at least once. Reuqired")
    markdup_parser.add_argument("--REMOVE_DUPLICATES", dest="removedups", type=bool, default=False,
                                help="If true do not write duplicates to the output file "
                                     "instead of writing them with appropriate flags set. [False]")
    markdup_parser.add_argument("-O", "--outfile", dest="outfile",
                                help="A BAM file's name for output results.")

    return markdup_parser


def bqsr_args(cdl_parser):
    """commandline arguments for BQSR by GATK4 BaseRecalibrator and ApplyBQSR.
    """
    bqsr_parser = cdl_parser.add_parser("BQSR", help=("Run BQSR by GATK4."))

    bqsr_parser.add_argument("-I", "--input", dest="inbam", action='append', default=[],
                             help="Input one or more BAM files to analyze. This "
                                     "argument must be specified at least once. Reuqired")
    bqsr_parser.add_argument("-O", "--outfile", dest="outfile",
                             help="A BAM file's name for output results.")

    return bqsr_parser


def haplotypecaller_args(cdl_parser):
    """HaplotypeCaller"""
    hc_parser = cdl_parser.add_parser("HaplotypeCaller", help=("Run HaplotypeCaller by GATK4."))

    hc_parser.add_argument("-I", "--input", dest="inbam", action='append', default=[],
                           help="Input one or more BAM files to analyze. This "
                                "argument must be specified at least once. Reuqired")
    hc_parser.add_argument("-O", "--outfile", dest="outfile",
                             help="File to which variants should be written.")
    hc_parser.add_argument("-L", "--intervals", action="append", dest="interval", default=[],
                           help="One genomic intervals over which to operate.")
    hc_parser.add_argument("-E", "--emit-ref-confidence", dest="emit_gvcf", type=bool, default=False,
                           help="Emitting gvcf.")

    return hc_parser


def genotypegvcfs_args(cdl_parser):
    genotype_parser = cdl_parser.add_parser("GenotypeGVCFs", help=("Run GenotypeGVCFs by GATK4."))

    genotype_parser.add_argument("-V", "--variant", dest="variants", action='append', default=[],
                                 help="Input one or more BAM files to analyze. This "
                                      "argument must be specified at least once. Reuqired")
    genotype_parser.add_argument("-L", "--intervals", action="append", dest="interval", default=[],
                                 help="One genomic intervals over which to operate.")
    genotype_parser.add_argument("-O", "--outfile", dest="outfile",
                           help="File to which variants should be written.")

    return genotype_parser


def vqsr_args(cdl_parser):
    vqsr_parser = cdl_parser.add_parser("VQSR", help=("Run VQSR by GATK4."))

    vqsr_parser.add_argument("-V", "--variant", dest="variants", action='append', default=[],
                             help="Input one or more BAM files to analyze. This "
                             "argument must be specified at least once. Reuqired")
    vqsr_parser.add_argument("--snp-max-gaussians", dest="snp_max_gaussion", type=int,
                             default=6, help="Max number of Gaussians for the positive model "
                                             "for SNP. Default: 6")
    vqsr_parser.add_argument("--indel-max-gaussians", dest="indel_max_gaussion", type=int,
                             default=4, help="Max number of Gaussians for the positive model "
                                             "for INDEL. Default: 4")
    vqsr_parser.add_argument("--snp-ts-filter-level", dest="snp_ts_filter_level", type=float,
                             default=99.0, help="The truth sensitivity level for SNP at which "
                                                "to start filtering. Default: 99.0")
    vqsr_parser.add_argument("--indel-ts-filter-level", dest="indel_ts_filter_level", type=float,
                             default=95.0, help="The truth sensitivity level for Indel at which "
                                                "to start filtering. Default: 95.0")
    vqsr_parser.add_argument("-O", "--outfile", dest="outfile",
                             help="File to which variants should be written.")

    return vqsr_parser


def mergevcfs_args(cdl_parser):
    mergevcfs_parser = cdl_parser.add_parser("MergeVcfs", help=("Run MergeVcfs by GATK4."))

    mergevcfs_parser.add_argument("-I", "--INPUT", dest="variants", action='append', default=[],
                             help="Input one or more VCF files to analyze. This "
                                  "argument must be specified at least once. Reuqired")
    mergevcfs_parser.add_argument("-O","--OUTPUT", dest="outfile",
                            help="The merged VCF or BCF file")
    return mergevcfs_parser