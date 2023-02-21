"""Check the content of input fastq list.

Author: Shujia Huang
Date: 2023-02-21
"""
import sys
import gzip


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> " \
                    f"{error}\n{context.strip()}: '{context_str.strip()}'"

    print(error_str, file=sys.stderr)
    sys.exit(1)


def check_samplesheet(input_fname):
    """This function checks that the samplesheet follows the following structure:

    #SAMPLE	RGID	FASTQ1	FASTQ2	LANE
    HL003   "@RG\tID:HL003_HL003_L01_65\tPL:COMPLETE\tPU:HL003_L01_65_HL003\tLB:HL003_L01\tSM:HL003"        HL003_L01_65_1.fq.gz     HL003_L01_65_2.fq.gz     HL003_L01_65
    HL003   "@RG\tID:HL003_HL003_L01_67\tPL:COMPLETE\tPU:HL003_L01_67_HL003\tLB:HL003_L01\tSM:HL003"        HL003_L01_67_1.fq.gz     .       HL003_L01_67
    HL003   "@RG\tID:HL003_HL003_L01_68\tPL:COMPLETE\tPU:HL003_L01_68_HL003\tLB:HL003_L01\tSM:HL003"        HL003_L01_68_1.fq.gz     HL003_L01_68_2.fq.gz     HL003_L01_68
    """

    with gzip.open(input_fname) if input_fname.endswith(".gz") else open(input_fname) as f:
        HEADER = ["#SAMPLE", "RGID", "FASTQ1", "FASTQ2", "LANE"]
        header = f.readline().strip().split()
        if header[:len(HEADER)] != HEADER:
            if header[: len(HEADER)] != HEADER:
                print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}",
                      file=sys.stderr)
                sys.exit(1)

        # Check sample
        for line in f:
            col = [x for x in line.strip().split() if x]

            # Check valid number of columns per row
            if len(col) < len(HEADER):
                print_error(f"Invalid number of columns (minimum = {len(HEADER)})!", "Line", line)

            # Check sample name
            sample_id, read_group_id, fq1, fq2, lane = col[:len(HEADER)]
            read_group_id = read_group_id.replace("'", '"').strip('"')

            if sample_id.find(" ") != -1:
                print_error(f"Spaces have been replaced by underscores for sample: '{sample_id}'")

            if not sample_id:
                print_error("Sample entry has not been specified!", "Line", line)

            # Check FastQ file extension
            if fq1 == ".":
                print_error("FASTQ1 should be a file", "Line", line)

            for fq in [fq1, fq2]:
                if fq == ".":
                    continue

                fq_fn = fq.replace(".gz", "")
                if not fq_fn.endswith(".fastq") and not fq_fn.endswith(".fq"):
                    print_error(f"{fq} file does not have extension '.fastq.gz', '.fq.gz', '.fastq' or '.fq'!",
                                "Line",
                                line)

            # Check Read group
            PLATFORM = ["ILLUMINA", "COMPLETE", "SOLID", "LS454", "HELICOS", "PACBIO"]
            if not read_group_id.startswith("@RG"):
                print_error("Read Groups must start with '@RG'", "Line", read_group_id)
            elif "ID:" not in read_group_id:
                print_error("Read Groups missing 'ID' field", "Line", line)
            elif "PL:" not in read_group_id:
                print_error("Read Groups missing 'PL' field", "Line", line)
            elif "PU:" not in read_group_id:
                print_error("Read Groups missing 'PU' field", "Line", line)
            elif "LB:" not in read_group_id:
                print_error("Read Groups missing 'LB' field", "Line", line)
            elif "SM:" not in read_group_id:
                print_error("Read Groups missing 'SM' field", "Line", line)

            platform = read_group_id.split("PL:")[-1].split("\t")[0].split("\\t")[0]
            if platform.upper() not in PLATFORM:
                print_error(f"Invalid value in 'PL' field of read group: '{platform}', which "
                            f"could only be one of: {','.join(PLATFORM)}", "Line", line)
    return True
