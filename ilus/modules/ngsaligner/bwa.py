"""NGS alignments with BWA
"""


def bwa_mem(config, out_prefix, rgID, fastq1, fastq2=""):
    """Perform piped alignment of fastq input files, generating sorted output BAM

    Parameters:
        ``aione``: A dict like.
            Data all in one by {key: value}
    """
    if fastq2 == ".":
        fastq2 = ""

    ref_index = config["resources"]["reference"]
    return _get_bwa_mem_cmd(fastq1, fastq2, ref_index, rgID, out_prefix, config)


def _get_bwa_mem_cmd(fastq1, fastq2, ref_index, rg_info, outfile_prefix, config):

    bwa = config["aligner"]["bwa"]
    bwa_options = " ".join([str(x) for x in config["aligner"].get("bwamem_options", [])])
    samtools = config["samtools"]["samtools"]
    samtools_options = " ".join([str(x) for x in config["samtools"].get("sort_options", [])])

    # skip seeds with more than INT occurrences
    bwa_cmd = ("time {bwa} mem {bwa_options} -R {rg_info} {ref_index} "
               "{fastq1} {fastq2} | {samtools} view -bS > {outfile_prefix}.bam && "
               "{samtools} sort {samtools_options} -o {outfile_prefix}.sorted.bam {outfile_prefix}.bam "
               "&& rm -rf {outfile_prefix}.bam").format(**locals())

    return "{outfile_prefix}.sorted.bam".format(**locals()), bwa_cmd
