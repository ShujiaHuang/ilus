"""
Author: Shujia Huang
Date: 2020-04-25
"""
import re
import os
import sys

pattern = re.compile(r'^\[\S+\]\s+(\S+)\s+done')


def check(input_fname):
    is_finish = False
    if os.path.exists(input_fname):
        with open(input_fname) as I:
            last_line = ""
            for line in I:
                last_line = line.strip()

            match = pattern.match(last_line)
            if match:
                is_finish = True

    return is_finish


if __name__ == "__main__":

    task_log_file = sys.argv[1]

    is_all_finish = True
    with open(task_log_file) as I:
        for line in I:
            # sample log_file task_shell_file"
            if line.startswith("#"):
                continue

            sample, log_fname, task_shell = line.strip().split()
            is_finish = check(log_fname)

            if not is_finish:
                is_all_finish = False
                print ("[unfinish]\t%s\t%s" % (sample, task_shell))

    if is_all_finish:
        print ("** All Jobs done **")
