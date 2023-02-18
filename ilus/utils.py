"""Useful utilities functions for building analysis pipeline
"""
import os
import sys
import re
import stat

from datetime import datetime


def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname

    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # Could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1

    return dname


def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False


def which(program):
    """Mimics behavior of UNIX which command.
    """
    # Add .exe program extension for windows support
    if os.name == "nt" and not program.endswith(".exe"):
        program += ".exe"

    envdir_list = [os.curdir] + os.environ["PATH"].split(os.pathsep)
    for envdir in envdir_list:
        program_path = os.path.join(envdir, program)
        if os.path.isfile(program_path) and os.access(program_path, os.X_OK):
            return program_path

    return


def split_jobs(input_shell_file, job_num, thread_num, prefix="work"):
    """Create sub jobs for input_shell_file."""
    if not os.path.isfile(input_shell_file):
        raise ValueError(f"{input_shell_file} is not a file.")

    # Read commands from input shell file
    with open(input_shell_file) as f:
        commands = [line.strip() for line in f if not line.startswith("#")]

    total_cmd_num = len(commands)
    # Ensure job_num is not greater than total_cmd_num
    job_num = min(job_num, total_cmd_num)

    # Calculate number of commands per job
    cmd_num_perjob, remainder = divmod(total_cmd_num, job_num)
    if remainder:
        cmd_num_perjob += 1

    # Create sub shell directory
    _, shell_fname = os.path.split(input_shell_file)
    now_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    sub_shell_dir = os.path.abspath(f"{shell_fname}.{now_time}")
    os.makedirs(sub_shell_dir, exist_ok=True)

    # Calculate number of parts per job
    part_num, remainder = divmod(cmd_num_perjob, thread_num)
    if remainder:
        part_num += 1

    # Calculate number of commands in the last job
    cmd_num_lastjob = total_cmd_num % cmd_num_perjob if total_cmd_num % cmd_num_perjob else cmd_num_perjob
    last_part_num = cmd_num_lastjob // thread_num + 1 \
        if cmd_num_lastjob % thread_num else cmd_num_lastjob // thread_num

    # Split commands into jobs and write to sub shell files
    for i, start_idx in enumerate(range(0, total_cmd_num, cmd_num_perjob)):
        end_idx = min(start_idx + cmd_num_perjob, total_cmd_num)
        sub_shell_fname = os.path.join(sub_shell_dir, f"{prefix}.{i + 1}.sh")

        # Determine number of parts in last job
        if end_idx >= total_cmd_num:
            part_num = last_part_num

        # Write commands to sub shell file
        with open(sub_shell_fname, "w") as f:
            f.write("#!/bin/bash\n")
            n = 0
            for j, cmd in enumerate(commands[start_idx:end_idx]):
                if thread_num > 1:
                    f.write(f"{cmd} &\n")
                else:
                    f.write(f"{cmd}\n")

                if (j + 1) % thread_num == 0:
                    n += 1
                    if thread_num > 1:
                        f.write("wait\n")
                    f.write(f'echo "----------- {n}/{part_num} ----------"\n')

            if (j + 1) % thread_num > 0:
                f.write("wait\n")
                f.write(f'echo "----------- {n + 1}/{part_num} ----------"\n')

        # Set file permissions to 0700
        os.chmod(sub_shell_fname, stat.S_IRWXU)

    return sub_shell_dir


def check_jobs_status(task_log_file):
    # Match anything looks like: "[xx] xxx done" or "[xx] xx xxx ... done"
    pattern = re.compile(r'^\[\S+\].*?\s+done$')

    def check(input_fname):
        if not os.path.exists(input_fname):
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
