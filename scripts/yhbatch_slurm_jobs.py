"""Create yhbatch jobs in GuangZhou TianHe2 computer cluster.
"""
import argparse
from ilus.modules.utils import split_jobs

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="yhbatch jobs")

    cmdparser.add_argument(
        "-I", "--input",
        dest="input",
        required=True,
        help="Input shell file."
    )
    cmdparser.add_argument(
        "-p", "--prefix",
        dest="prefix",
        type=str,
        default="work",
        help="The prefix name for output sub-shell. (default: %(default)s)"
    )
    cmdparser.add_argument(
        "-n", "--number",
        dest="number",
        type=int,
        required=True,
        help="Number of sub job."
    )
    cmdparser.add_argument(
        "-t", "--parallel",
        dest="t",
        type=int,
        required=True,
        help="Parallel number for per sub job."
    )

    args = cmdparser.parse_args()
    split_jobs(args.input, args.number, args.t, prefix=args.prefix)
