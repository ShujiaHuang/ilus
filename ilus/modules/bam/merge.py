"""Run functions by provided a name with arguments.
"""
import sys

from ilus.tools.sambamba import SambambaRunner, sambamba_merge


def mergebamfiles(kwargs, aione):
    if not kwargs["args"].outfile:
        sys.stderr.write("Error: missing the output BAM files when running "
                         "mergebamfiles.\n")
        sys.exit(1)

    sambamba = SambambaRunner(aione["config"])
    aione["merge_bam"] = sambamba_merge(sambamba,
                                        kwargs["args"].inbam,
                                        kwargs["args"].outfile)
    return aione