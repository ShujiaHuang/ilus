"""Work with sambamba commandline
"""

from _sambamba import *


def sambamba_mark_duplicates(sambamba, align_bams, outbamfile, remove_dups=False):
    """Run GATK MarkDuplicates.

    Parameters:
        ``gatk``: GATKRunner
        ``align_bams``: A list like
            contain at least one BAM files
    """
    options = ["-t %d" % sambamba.num_cores]
    if remove_dups:
        options.append("-r")

    options += align_bams + [outbamfile]

    sambamba.run("markdup", options)
    return outbamfile


def sambamba_merge(sambamba, align_bams, outbamfile):
    options = ["-t %d" % sambamba.num_cores]
    options.append(outbamfile)
    options += align_bams

    sambamba.run("merge", options)
    return outbamfile


def sambamba_index(sambamba, align_bam):
    options = ["-t %d" % sambamba.num_cores, align_bam]
    sambamba.run("index", options)

    return