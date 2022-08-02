"""
Author: Shujia Huang
Date: 2020-04-25
"""
import argparse
from ilus.utils import check_jobs_status

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="Check the jobs have finished or not.")
    cmdparser.add_argument("-I", "--input", dest="input", required=True,
                           help="Input the log file of task.")

    args = cmdparser.parse_args()
    check_jobs_status(args.input)
