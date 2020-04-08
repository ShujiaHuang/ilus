"""
Get runnable regions by divide N region alone the reference.
Author: Shujia Huang
Date: 2018-09-07
"""
import sys
import optparse

MAXSIZE = 4550000 # 4.6Mb
SLOWCHUNKSIZE = 500000  # 0.5 Mb

def main(opts):
    """
    """
    fa = get_fasta_size_from_faifile(opts.ref_fai_file)

    # get non N region
    call_region = get_nonN_region(opts.n_reg_file, opts.nsize, fa)
    for chrom, regions in sorted(call_region.items(), key = lambda x:(x[0])):
        # if chrom[:2] == 'GL': continue
        for start_pos, end_pos in regions:
            print "\t".join(map(str, [chrom, start_pos, end_pos]))


def get_fasta_size_from_faifile(infile):
    """
    Get chromosome size from .fai file

    .fai format:
        1       249250621       52      60      61
        2       243199373       253404903       60      61
        3       198022430       500657651       60      61
        4       191154276       701980507       60      61
        5       180915260       896320740       60      61
        6       171115067       1080251307      60      61
    """

    fa = {}
    with open(infile) as I:
        for r in I:
            col = r.strip().split()
            fa[col[0]] = int(col[1])

    return fa


def get_nonN_region(infile, thd_n_size, fa):
    """
    Read N region to get non-n region
    """
    reg, pre_pos = {}, {}
    pre_id, pre_end = '', 0
    with open(infile) as I:
        # dealing with N region file

        for r in I:

            chromid, n_start, n_end = r.strip().split()[0:3]
            n_start, n_end = map(int, (n_start, n_end))
            n_region_size = n_end-n_start+1

            pre_pos[chromid] = n_end

            if n_start == 1 and n_end == fa[chromid]:
                # The whole chromsome is all N -- Probably impossible
                continue

            if chromid not in reg:

                if len(pre_id) and pre_pos[pre_id] < fa[pre_id]:
                    reg[pre_id][-1][1] = fa[pre_id]

                reg[chromid] = []
                if n_start > 1:
                    reg[chromid] = [[1, n_start-1]]
                elif n_end+1 < fa[chromid]:
                    reg[chromid] = [[n_end+1, 0]]

                pre_id = chromid

            if n_region_size < thd_n_size:
                # It's a small N-region keep it.

                if reg[chromid][-1][1] != 0:
                    reg[chromid][-1][1] = n_end

            else:

                if reg[chromid][-1][1] == 0:
                    reg[chromid][-1][1] = n_start-1

                if n_end+1 == reg[chromid][-1][0]:
                    continue

                if n_end + 1 < fa[chromid]:
                    reg[chromid].append([n_end+1, 0])

    if pre_id in reg and pre_pos[pre_id] < fa[pre_id]:
        reg[pre_id][-1][1] = fa[pre_id]

    for chrom, chromsize in fa.items():
        # loading all other non N regions
        if chrom not in reg:
            reg[chrom] = [[1, chromsize]]

    return reg


def overlap_region(start1, end1, start2, end2):
    """
    """
    s, e = 0, 0
    if end1 <= end2 and start1 <= start2:
        s, e = start2, end1

    elif end1 <= end2 and start1 > start2:
        if end1 < start1:
            print >> sys.stderr, '[ERROR] Start1 > end1:', start1, end1
            sys.exit(1)
        s, e = start1, end1

    elif start1 <= start2 and end1 > end2:
        if end2 < start2:
            print >> sys.stderr, '[ERROR] Start2 > end2:', start2, end2
            sys.exit(1)
        s, e = start2, end2

    elif start1 <= end2 and end1 > end2:
        s, e = start1, end2

    else:
        print >> sys.stderr, '[ERROR] Region not right.'

    return s, e


if __name__ == '__main__':

    """
    loading all the command options.
    """
    usage = '\nUsage: python %prog [options] -f <fasta_faifile> > Output'
    optp = optparse.OptionParser(usage=usage)
    optp.add_option('-f', '--fai', dest='ref_fai_file', metavar='STR',
                    help='The reference fasta fai file.', default='')
    optp.add_option('-n', '--nbed', dest='n_reg_file', metavar='STR',
                    help='The N-base region bed file.', default='')
    optp.add_option('-N', '--nbase', dest='nsize', metavar='int',
                    help='The biggest continue N size could be allowed '
                         'in one region.', default=1000)

    opt, _ = optp.parse_args()
    if not opt.ref_fai_file: 
        optp.error('Required [-f ref_fai_file]\n')

    main(opt)
    print >> sys.stderr, '>> For the flowers bloom in the desert <<'





