"""
Get runnable regions by divide N region alone the reference.
Author: Shujia Huang
Date: 2018-09-07
"""
import sys
import argparse


def get_fasta_size_from_faifile(infile):
    """
    Get chromosome size from .fai file

    .fai format:
        chr1	248956422	112	70	71
        chr2	242193529	252513167	70	71
        chr3	198295559	498166716	70	71
        chr4	190214555	699295181	70	71
        chr5	181538259	892227221	70	71
        chr6	170805979	1076358996	70	71
    """

    fa = {}
    chrom_order_list = []
    with open(infile) as I:
        for r in I:
            col = r.strip().split()
            fa[col[0]] = int(col[1])
            chrom_order_list.append(col[0])

    return chrom_order_list, fa


def get_nonN_region(in_n_interval_file, thd_n_size, fa):
    """Read N region to get non-n region

    ``in_n_interval_file`` format:

        chr1	1	10000
        chr1	207667	257666
        chr1	297969	347968
        chr1	535989	585988
        chr1	2702782	2746290
    """
    # New robust implementation:
    # 1. Read N intervals per chromosome, normalize starts (treat start<=0 as 1).
    # 2. Filter N intervals by threshold: only keep those with size >= thd_n_size
    #    as "true N" regions that should be removed. Small N regions are
    #    ignored (treated as non-N).
    # 3. Merge overlapping/adjacent true N regions and compute the complement
    #    (non-N regions) for each chromosome.

    from collections import defaultdict

    ns = defaultdict(list)

    with open(in_n_interval_file) as I:
        for line in I:
            parts = line.strip().split()
            if not parts:
                continue
            chromid = parts[0]
            try:
                n_start = int(parts[1])
                n_end = int(parts[2])
            except (IndexError, ValueError):
                continue

            # Normalize: treat starts <= 0 as 1 (allow 0-based input)
            if n_start < 1:
                n_start = 1

            if n_end < n_start:
                # malformed line
                continue

            n_region_size = n_end - n_start + 1
            ns[chromid].append((n_start, n_end))

    reg = {}

    for chrom, chrom_size in fa.items():
        intervals = ns.get(chrom, [])
        if not intervals:
            reg[chrom] = [[1, chrom_size]]
            continue

        intervals.sort()

        # If there are N regions covering the chromosome start, advance cur
        # past their end so that non-N starts after them (this preserves the
        # original behavior where a leading N at the chromosome start shifts
        # the first callable region).
        cur = 1
        i = 0
        while i < len(intervals) and intervals[i][0] <= cur:
            # intervals starting at or before current pointer
            cur = max(cur, intervals[i][1] + 1)
            i += 1

        # Now collect "true" N intervals (size >= threshold) that occur at or
        # after the current pointer; these will split callable regions.
        true_ns = []
        for s, e in intervals:
            if e < cur:
                continue
            if e - s + 1 >= thd_n_size:
                # trim start to be at least cur
                s = max(s, cur)
                true_ns.append((s, e))

        # Merge true N intervals
        merged = []
        for s, e in sorted(true_ns):
            if not merged:
                merged.append([s, e])
            else:
                if s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])

        # Compute complement (non-N regions) starting from cur
        nonn = []
        if not merged:
            if cur <= chrom_size:
                nonn.append([cur, chrom_size])
        else:
            if cur <= merged[0][0] - 1:
                nonn.append([cur, merged[0][0] - 1])

            for j in range(len(merged) - 1):
                a_end = merged[j][1]
                b_start = merged[j + 1][0]
                if a_end + 1 <= b_start - 1:
                    nonn.append([a_end + 1, b_start - 1])

            if merged[-1][1] < chrom_size:
                nonn.append([merged[-1][1] + 1, chrom_size])

        reg[chrom] = nonn

    return reg


def overlap_region(start1, end1, start2, end2):
    """Get the intersection region
    """
    s, e = 0, 0
    if end1 <= end2 and start1 <= start2:
        s, e = start2, end1

    elif end1 <= end2 and start1 > start2:
        if end1 < start1:
            sys.stderr.write("[ERROR] Start1(%d) > end1(%d)" % (start1, end1))
            sys.exit(1)

        s, e = start1, end1

    elif start1 <= start2 and end1 > end2:
        if end2 < start2:
            sys.stderr.write("[ERROR] Start2(%d) > end2(%d)\n" % (start2, end2))
            sys.exit(1)

        s, e = start2, end2

    elif (start1 > start2) and (start1 <= end2 < end1):
        s, e = start1, end2

    else:
        sys.stderr.write("[ERROR] The input region is not overlap with each other: "
                         "[%d-%d, %d-%d].\n" % (start1, end1, start2, end2))

    return s, e


def main(opts):
    chrom_order_list, fa = get_fasta_size_from_faifile(opts.ref_fai_file)

    if opts.n_reg_file:
        call_region = get_nonN_region(opts.n_reg_file, opts.nsize, fa)
    else:
        # The whole genome region without removing N regions.
        call_region = {c: [[1, s]] for c, s in fa.items()}

    for chrom in chrom_order_list:
        for start_pos, end_pos in call_region[chrom]:
            for pos in range(start_pos, end_pos + 1, opts.win):
                s = pos
                e = end_pos if pos + opts.win >= end_pos else pos + opts.win - 1

                print("\t".join(map(str, [chrom, s, e])))


if __name__ == "__main__":
    usage = "\nUsage: python %prog [options] -f <fasta.fai file> [-n N region bed file] > Output"
    cmdparse = argparse.ArgumentParser(description=usage)
    cmdparse.add_argument("-f", "--fai", dest="ref_fai_file", required=True,
                          help="The reference FASTA .fai file.")
    cmdparse.add_argument("-n", "--nbed", dest="n_reg_file", required=False,
                          help="The N region .bed file.")
    cmdparse.add_argument("-N", "--nbase", dest="nsize", type=int, default=1000,
                          help="The biggest continue N size could be allowed "
                               "in one region. [1000]")
    cmdparse.add_argument("-w", "--win", dest="win", type=int, default=5000000,
                          help="The max window size. [5000000]")
    opt = cmdparse.parse_args()

    main(opt)
    sys.stderr.write("** For the flowers bloom in the desert **\n")
