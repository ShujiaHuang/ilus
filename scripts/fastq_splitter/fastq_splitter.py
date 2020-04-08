"""
FASTQ Splitter  -  a script for partitioning a FASTQ file into pieces

Version 0.1.1 (Aug 16, 2018)
Copyright (c) 2018 Shujia Huang
"""
import argparse
import os
import sys

from datetime import datetime

import gzip

version_description = "FASTQ Splitter 0.1.1 Copyright (c) 2018 Shujia Huang"

LINE_COUNT = 0
START_TIME = datetime.now()


def split_file(infile, n_parts, outdir):

    if not os.path.exists(infile):
        sys.stderr.write("Error: Can't find file: %s\n" % infile)
        sys.exit(1)

    fqname, ext = '', ''
    if infile.endswith(".fastq.gz"):
        fqname = os.path.basename(infile).split(".fastq.gz")[0]
        ext = "fastq.gz"
    elif infile.endswith(".fq.gz"):
        fqname = os.path.basename(infile).split(".fq.gz")[0]
        ext = "fq.gz"
    elif infile.endswith(".fastq"):
        fqname = os.path.basename(infile).split(".fastq")[0]
        ext = "fastq"
    elif infile.endswith(".fq"):
        fqname = os.path.basename(infile).split(".fq")[0]
        ext = "fq"
    else:
        sys.stderr.write("Error: The input files are not fastq format(*.fq.gz/*.fq/*.fastq.gz/*.fastq)\n")

    total_read_num, total_base_num = get_file_size(infile)

    elapsed_time = datetime.now() - START_TIME    
    print "Loaded %s: %d sequences, %d bp, %d seconds elapsed" % (infile, total_read_num, total_base_num, elapsed_time.seconds)
    print "=> dividing into %d parts:" % n_parts

    read_num_per_file = total_read_num/n_parts if total_read_num % n_parts == 0 else int(total_read_num/n_parts)+1
    num_len = len(str(n_parts))
    with gzip.open(infile) if infile.endswith(".gz") else open(infile) as I:
        for part in range(1, n_parts+1):
            part_file = "%s.%0*d.%s" % (fqname, num_len, part, ext)
            out_sub_file = '/'.join([outdir, part_file])
            
            print out_sub_file

            written = 0
            with gzip.open(out_sub_file, "wb") if out_sub_file.endswith(".gz") else open(out_sub_file, "w") as OUT:
                is_done = False
                while not is_done and written < read_num_per_file:
                    written += 1
                    is_done, _, read = get_fastq_read(I)
                    OUT.write("%s\n" % read)


def get_file_size(infile):

    total_read_num, total_base_num = 0, 0
    with gzip.open(infile) if infile.endswith(".gz") else open(infile) as I:

        is_done = False
        while not is_done:
            is_done, read_len, read = get_fastq_read(I)
            if not is_done:
                total_read_num += 1
                total_base_num += read_len

    return total_read_num, total_base_num


def get_fastq_read(iter_fh):
    truncated = "Error: Incomplete FASTQ entry at the end -> " \
                "looks like truncated input!"

    global LINE_COUNT

    name_line = fetch_next(iter_fh).strip()
    is_done = False if name_line else True
    if is_done:
        return True, 0, ""

    LINE_COUNT += 1
    if name_line[0] != "@":
        sys.stderr.write("Error parsing line %d: FASTQ entry does not start "
                         "with '@':\n%s\n" % (LINE_COUNT, name_line))

    read_seq = fetch_next(iter_fh).strip()
    LINE_COUNT += 1
    is_done = False if read_seq else True

    plus_line = fetch_next(iter_fh).strip()
    LINE_COUNT += 1
    is_done = False if plus_line else True

    if plus_line[0] != "+":
        sys.stderr.write("Error parsing line %d: Expecting '+', "
                         "found '%s'" % (LINE_COUNT, plus_line))

    read_qual = fetch_next(iter_fh).strip()
    LINE_COUNT += 1
    is_done = False if read_qual else True
    if len(read_seq) != len(read_qual):
        sys.stderr.write("Error: Misformatted FASTQ entry in input line %d: "
                         "quality length (%d) differs from sequence length "
                         "(%d):\n%s\n%s\n" % (LINE_COUNT, len(read_seq),
                                              len(read_qual), read_seq,
                                              read_qual))

    read = "\n".join([name_line, read_seq, plus_line, read_qual])
    return is_done, len(read_seq), read


def fetch_next(iter_fh):

    if iter_fh == '': return ''

    try:
        line = iter_fh.next()
    except StopIteration:
        line = ''

    return line


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=version_description)
    parser.add_argument("-i", "--infastq", dest='fastqfile',
                        help="Input fastq(or .fastq.gz) file")
    parser.add_argument("-n", "--n-parts", dest='nparts', type=int,
                        help="Divid into <n> parts. required")
    parser.add_argument("-d", "--outdir", dest="outdir",
                        help="Output directory. required")

    args = parser.parse_args()

    if args.outdir is None:
        parser.error("Splitting directory is not specified ('--outdir' option)\n")

    if args.nparts is None:
        parser.error("Divide number is not specified ('--n-parts' option)\n")

    if args.nparts <= 0:
        parser.error("Non-positive number of parts.\n")

    if args.nparts <= 1:
        sys.stderr.write("WARNING: Divide into %d parts, means we don't have to divide anything.\n" % args.nparts)
        sys.exit(0)

    args.fastqfile = os.path.abspath(args.fastqfile)
    args.outdir = os.path.abspath(args.outdir)

    split_file(args.fastqfile, args.nparts, args.outdir)

    elapsed_time = datetime.now() -START_TIME 
    print "\n** Fastq splitter done, %d seconds elapsed **\n" % elapsed_time.seconds












