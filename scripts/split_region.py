import sys
import argparse

def main(opts):

    with open(opts.in_bed_fname, "r") as I:
        for line in I:
            # chr3	91276995	93655574
            # chr4	47839017	49336924
            if line.startswith("#"):
                continue

            col = line.strip().split()
            chrom, start_pos, end_pos = col[0], int(col[1]), int(col[2])
            for pos in range(start_pos, end_pos + 1, opts.win):
                s = pos
                e = end_pos if pos + opts.win >= end_pos else pos + opts.win - 1

                print("\t".join(map(str, [chrom, s, e])))


if __name__ == '__main__':
    usage = '\nUsage: python %prog [options]  > Output'
    cmdparse = argparse.ArgumentParser(description=usage)

    cmdparse.add_argument('-I', '--bed', dest='in_bed_fname', metavar='STR', required=True,
                          help='The region bed file.')
    cmdparse.add_argument('-w', '--win', dest='win', metavar='int', type=int,
                          help='The max window size.', default=100000)

    opt = cmdparse.parse_args()

    main(opt)
    sys.stderr.write('** For the flowers bloom in the desert **\n')