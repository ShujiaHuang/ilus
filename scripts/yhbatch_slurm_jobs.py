"""Create yhbatch jobs in GuangZhou TianHe2 computer cluster.
"""
import argparse
from ilus.utils import split_jobs

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="yhbatch jobs")
    cmdparser.add_argument("-I", "--input", dest="input", required=True,
                           help="Input shell file.")
    cmdparser.add_argument("-p", "--prefix", dest="prefix", type=str, default="work",
                           help="The prefix name of output sub-shell file. [work]")
    cmdparser.add_argument("-n", "--number", dest="number", type=int, required=True,
                           help="The number of sub shells.")
    cmdparser.add_argument("-t", "--parallel", dest="t", type=int, required=True,
                           help="The number of parallel tasks.")

    args = cmdparser.parse_args()
    split_jobs(args.input, args.number, args.t, prefix=args.prefix)
