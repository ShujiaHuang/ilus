"""Useful utilities functions for building analysis pipeline
"""
import stat

from pathlib import Path
from datetime import datetime

from ._utils import safe_makedir


def split_jobs(input_shell_file, job_num, thread_num, prefix="work"):
    """Create sub jobs for input_shell_file."""
    if not Path(input_shell_file).is_file():
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
    shell_fname = Path(input_shell_file).name
    now_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    sub_shell_dir = Path(f"{shell_fname}.{now_time}").resolve()
    safe_makedir(sub_shell_dir, exist_ok=True)

    # Calculate number of parts per job
    part_num, remainder = divmod(cmd_num_perjob, thread_num)
    if remainder:
        part_num += 1

    # Calculate number of commands for the last job
    cmd_num_lastjob = total_cmd_num % cmd_num_perjob \
        if total_cmd_num % cmd_num_perjob else cmd_num_perjob
    last_part_num = cmd_num_lastjob // thread_num + 1 \
        if cmd_num_lastjob % thread_num else cmd_num_lastjob // thread_num

    # Split commands into jobs and write to sub shell files
    for i, start_idx in enumerate(range(0, total_cmd_num, cmd_num_perjob)):
        end_idx = min(start_idx + cmd_num_perjob, total_cmd_num)
        sub_shell_fname = sub_shell_dir.joinpath(f"{prefix}.{i + 1}.sh")

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
        sub_shell_fname.chmod(stat.S_IRWXU)

    return sub_shell_dir
