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


def split_jobs(args):

    commands = []
    with open(args.input) as I:
        for line in I:
            if line.startswith("#"):
                continue

            commands.append(line.strip())

    total_cmd_num = len(commands)
    if args.number > total_cmd_num:
        args.number = total_cmd_num

    cmd_num_perjob = total_cmd_num // args.number
    if total_cmd_num % args.number > 0:
        cmd_num_perjob += 1

    _, shell_fname = os.path.split(args.input)
    now_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    sub_shell_dir = os.path.abspath("{shell_fname}.{now_time}".format(**locals()))
    if not os.path.exists(sub_shell_dir):
        os.makedirs(sub_shell_dir)

    part_num = cmd_num_perjob // args.t
    if cmd_num_perjob % args.t > 0:
        part_num += 1

    for i, k in enumerate(range(0, total_cmd_num, cmd_num_perjob)):
        sub_shell_fname = os.path.join(sub_shell_dir, "work.%d.sh" % (i + 1))
        with open(sub_shell_fname, "w") as OUT:

            n = 0
            OUT.write("#!/bin/bash\n")
            for j, cmd in enumerate(commands[k:k + cmd_num_perjob]):
                OUT.write("%s &\n" % cmd if args.t > 1 else "%s\n" % cmd)
                if (j + 1) % args.t == 0:
                    n += 1
                    if args.t > 1:
                        OUT.write("wait\n")
                    OUT.write("echo \"----------- %d/%d ----------\"\n" % (n, part_num))

            if (j + 1) % args.t > 0:
                OUT.write("wait\n")
                OUT.write("echo \"----------- %d/%d ----------\"\n" % (n + 1, part_num))

        os.chmod(sub_shell_fname, stat.S_IRWXU)  # 0700


def check_jobs_status(args):

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

    # Match anything looks like: "[xx] xxx done" or "[xx] xx xxx ... done"
    pattern = re.compile(r'^\[\S+\]\s+\S+.*?\s+done$')

    task_log_file = args.input
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
                print("[unfinish]\t%s\t%s" % (sample, task_shell))

    if is_all_finish:
        sys.stderr.write("** All Jobs done **\n")
