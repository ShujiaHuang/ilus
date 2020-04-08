"""
Find the Low Speed Region along the human genome.
Cut genome region from N.
/zfssz3/SOLEXA_DATA/BC_RD_P2/USER/huangshujia/ifs5_BC_ONLINE_RD/project/201605_ODPS4WGS/20160626_find_lowspeed_region/bin
Author: Shujia Huang
Date: 2016-06-26
"""
import os
import sys
import optparse

MAXSIZE = 4550000 # 4.6Mb
SLOWCHUNKSIZE = 500000  # 0.5 Mb

def main(opts):
    """
    """
    #accbed, accidx = load_accbed(opts.accbed_file)
    fa = get_fasta_size_from_faifile(opts.ref_fa_file)

    # get non N region
    call_region = get_nonN_region(opts.n_reg_file, opts.nsize, fa)
    slow_reg, s_idx = load_slow(opts.slow_region_file)

    for chrom in call_region:
        if chrom[:2] == 'GL': continue

        for start_pos, end_pos in call_region[chrom]:

            reg = []
            for pos in range(start_pos, end_pos+1, MAXSIZE):
                s = pos
                e = end_pos if pos + MAXSIZE >= end_pos else pos+MAXSIZE-1

                flag = True
                ovlp_regs = []
                if chrom in s_idx:
                    for i in range(s_idx[chrom], len(slow_reg[chrom])):
                        
                        if s > slow_reg[chrom][i][1]: 
                            continue

                        if e < slow_reg[chrom][i][0]: 
                            break

                        if flag:
                            flag = False
                            s_idx[chrom] = i

                        o_s, o_e = overlap_region(s, e, slow_reg[chrom][i][0],
                                                  slow_reg[chrom][i][1])
                        ovlp_regs.append([o_s, o_e])


                if len(ovlp_regs) > 0:

                    is_first = True
                    pre_pos = s

                    for start, end in ovlp_regs:

                        if pre_pos < start:
                            reg.append([chrom, pre_pos, start - 1])

                        #reg.append([chrom, start, end])
                        for j in range(start-1, end, SLOWCHUNKSIZE):
                            tmp_s = j
                            if j + SLOWCHUNKSIZE > end:
                                tmp_e = end
                            else:
                                tmp_e = j + SLOWCHUNKSIZE
                            reg.append([chrom, tmp_s+1, tmp_e])

                        pre_pos = end + 1

                    if pre_pos <= e:
                        reg.append([chrom, pre_pos, e])

                else:
                    reg.append([chrom, s, e])

            for c, s, e in reg:
                print '%s\t%d\t%d' % (c, s, e)


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


def load_accbed(infile):
    """
    """
    region, index = {}, {}
    with open(infile) as I:
        for r in I:
            col = r.strip('\n').split()
            if col[0] not in region:
                region[col[0]] = []
                index[col[0]] = 0

            region[col[0]].append([col[0], int(col[1]), int(col[2])])

    return region, index


def get_nonN_region(infile, thd_n_size, fa):
    """
    Read N region to get non-n region
    """
    reg, pos = {}, {}
    pre_id, pre_end, curr_start = '', 0, 0
    with open(infile) as I:
        # dealing with N region file

        for r in I:

            col = r.strip().split()
            n_start, n_end = map(int, (col[1], col[2]))
            n_region_size = n_end-n_start+1

            if col[0] not in reg:

                reg[col[0]] = []
                pos[col[0]] = 0

                if len(pre_id):
                    reg[pre_id].append([pos[pre_id]+1, fa[pre_id]])

                pre_id, pre_end = col[0], n_end

            elif n_region_size < thd_n_size:
                # ignore this n region
                pass

            else:
                reg[pre_id].append([pos[pre_id]+1, n_start-1])
                pre_id, pre_end = col[0], n_end

            if n_region_size >= thd_n_size:
                pos[col[0]] = n_end

    for chrom, chromsize in fa.items():
        # loading all other non N regions
        if chrom not in reg:
            reg[chrom] = [1, chromsize]

    return reg


def load_slow(infile):
    """
    """
    reg, idx = {}, {}
    with open(infile) as I:
        for line in I:
            # Y  11234914  14674573
            if line[0] == '#': continue
            col = line.strip('\n').split()
            if col[0] not in reg:
                reg[col[0]] = []
                idx[col[0]] = 0

            reg[col[0]].append(map(int, [col[1], col[2]]))

    return reg, idx

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
    usage = '\nUsage: python %prog [options] -f <fastafile> > Output'
    optp = optparse.OptionParser(usage=usage)
    optp.add_option('-f', '--fa', dest='ref_fa_file', metavar='STR',
                    help='The reference fasta file.', default='')
    #optp.add_option('-b', '--bed', dest='accbed_file', metavar='STR',
    #                help = 'The accessible region bed file.', default='')
    optp.add_option('-s', '--slow', dest='slow_region_file', metavar='STR',
                    help='The alow region bed file.', default='')
    optp.add_option('-n', '--nbed', dest='n_reg_file', metavar='STR',
                    help='The N-base region bed file.', default='')
    optp.add_option('-N', '--nbase', dest='nsize', metavar='int',
                    help='The biggest continue N size could be allow '
                         'in one region.', default=1000)
    #optp.add_option('-w', '--win', dest='win', metavar='int',
    #                help = 'The window size.', default=100000)

    opt, _ = optp.parse_args()
    if not opt.ref_fa_file: 
        optp.error('Required [-f ref_fa_file]\n')

    if not opt.slow_region_file:
        optp.error('Required [-s slow_region_file]\n')
    #if not opt.accbed_file: optp.error('Required [-b accessible bed file.\n]')

    main(opt)
    print >> sys.stderr, '>> For the flowers bloom in the desert <<'





