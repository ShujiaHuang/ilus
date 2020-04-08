"""
Fastq basecontent

Version 0.1.1 (Aug 16, 2018)
Copyright (c) 2018 Shujia Huang
"""
import argparse
import os
import sys
import gzip

from datetime import datetime


version_description = "FASTQ Basecontent 0.1.1 Copyright (c) 2018 Shujia Huang"

LINE_COUNT = 0
START_TIME = datetime.now()


def main(infile):

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

    total_read_num, total_base_num, fastqcontent = get_base_content(infile)
    for k, v in fastqcontent.items():

        print "#%s" % k
        for kk, bb in v.items():
            print "%s\t%s" % (kk, bb)
        print "\n"

    elapsed_time = datetime.now() - START_TIME    
    print "Loaded %s: %d sequences, %d bp, %d seconds elapsed" % (infile, total_read_num, total_base_num, elapsed_time.seconds)

    
def get_base_content(infile):


    fastqcontent = {"Base":{"A": 0, "C": 0, "G": 0, "T": 0},
                    "Qual":{},
                    "ReadID":{},
                    "Plus_line":{}}

    total_read_num, total_base_num = 0, 0
    with gzip.open(infile) if infile.endswith(".gz") else open(infile) as I:

        is_done = False
        while not is_done:
            is_done, name_line, read_seq, plus_line, read_qual = get_fastq_read(I)
            if not is_done and total_read_num < 1000000:
                total_read_num += 1
                total_base_num += len(read_seq)

                for name, line in zip(("ReadID", "Base", "Plus_line", "Qual"),
                                   (name_line, read_seq, plus_line, read_qual)):
                    for b in line:
                        fastqcontent[name][b] = fastqcontent[name].get(b, 0) + 1

    return total_read_num, total_base_num, fastqcontent


def get_fastq_read(iter_fh):

    global LINE_COUNT

    name_line = fetch_next(iter_fh).strip()
    is_done = False if name_line else True
    if is_done:
        return True, "", "", "", ""

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

    return is_done, name_line, read_seq, plus_line, read_qual


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

    args = parser.parse_args()

    args.fastqfile = os.path.abspath(args.fastqfile)

    main(args.fastqfile)

    elapsed_time = datetime.now() -START_TIME 
    print "\n** Fastq splitter done, %d seconds elapsed **\n" % elapsed_time.seconds












