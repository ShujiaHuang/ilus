"""NGS alignments with BWA
"""
from pathlib import Path


def filter(config, outdir, input_fastq1, input_fastq2=None):
    """preprocessing normal Fastq files
    """
    if input_fastq2 == ".":
        input_fastq2 = None

    soapnuke = config["cleanfastq"]["soapnuke"]
    soapnuke_options = " ".join([str(x) for x in config["cleanfastq"].get("soapnuke_options", [])])

    # The name of clean fq1
    fq1_prefix = str(Path(input_fastq1).stem).replace('fastq', 'fq').split('.fq')[0]
    clean_fq1 = f"{fq1_prefix}.clean.fq.gz"  # Basename of clean fq1
    if input_fastq2:  # PE data
        # The name of clean fq2
        fq2_prefix = str(Path(input_fastq2).stem).replace('fastq', 'fq').split('.fq')[0]
        clean_fq2 = f"{fq2_prefix}.clean.fq.gz"  # basename of clean fq2
        soapnuke_filter_cmd = (f"time {soapnuke} filter {soapnuke_options} "
                               f"-1 {input_fastq1} -2 {input_fastq2} "
                               f"-C {clean_fq1} -D {clean_fq2} "
                               f"--outDir {outdir}")
    else:
        soapnuke_filter_cmd = (f"time {soapnuke} filter {soapnuke_options} "
                               f"-1 {input_fastq1}"
                               f"-C {clean_fq1} "
                               f"--outDir {outdir}")
        clean_fq2 = ""

    return (soapnuke_filter_cmd,
            str(Path(outdir).joinpath(clean_fq1)),  # set full path for clean fq
            str(Path(outdir).joinpath(clean_fq2)) if clean_fq2 else "")

