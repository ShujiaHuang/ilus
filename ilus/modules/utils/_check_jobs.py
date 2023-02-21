"""Function for checking jobs of pipeline are done or not.
"""
import sys
import re

from pathlib import Path


def check_jobs_status(task_log_file):
    # Match anything looks like: "[xx] xxx done" or "[xx] xx xxx ... done"
    pattern = re.compile(r'^\[\S+\].*?\s+done$')

    def check(input_fname):
        if not Path(input_fname).exists():
            return False

        with open(input_fname) as I:
            for last_line in I:
                pass

        return pattern.match(last_line.strip()) is not None

    all_done = True
    with open(task_log_file) as I:
        for line in I:
            # sample log_file task_shell_file"
            if line.startswith("#"):
                continue

            sample, log_fname, task_shell = line.strip().split()
            if not check(log_fname):
                all_done = False
                print(f"[unfinish]\t{sample}\t{task_shell}")

    if all_done:
        print("** All Jobs done **\n", file=sys.stderr)
