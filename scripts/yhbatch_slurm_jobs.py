"""Create yhbatch jobs in GuangZhou TianHe2 computer cluster.
"""
import os
import stat
import argparse

from datetime import datetime


if __name__ == "__main__":

    cmdparser = argparse.ArgumentParser(description="yhbatch jobs")
    cmdparser.add_argument("-I", "--input", dest="input", required=True,
                           help="Input shell file.")
    cmdparser.add_argument("-n", "--number", dest="number", type=int, required=True,
                           help="The number of sub shells.")
    cmdparser.add_argument("-t", "--parallel", dest="t", type=int, required=True,
                           help="The number of parallel tasks.")

    args = cmdparser.parse_args()

    commands = []
    with open(args.input) as I:
        for line in I:
            if line.startswith("#"):
                continue

            commands.append(line.strip())

    total_cmd = len(commands)
    if args.number > total_cmd:
        args.number = total_cmd

    inter_step_size = int(total_cmd/args.number)
    if inter_step_size * args.number < total_cmd:
        inter_step_size += 1

    _, shell_fname = os.path.split(args.input)
    now_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    sub_shell_dir = os.path.abspath("{shell_fname}.{now_time}".format(**locals()))
    if not os.path.exists(sub_shell_dir):
        os.makedirs(sub_shell_dir)

    for i, k in enumerate(range(0, total_cmd, inter_step_size)):
        sub_shell_fname = os.path.join(sub_shell_dir, "work.%d.sh" % (i+1))
        with open(sub_shell_fname, "w") as OUT:
            OUT.write("#!/bin/bash\n")

            n = 0
            for j, cmd in enumerate(commands[k:k+inter_step_size]):
                OUT.write("%s &\n" % cmd if args.t > 1 else "%s\n" % cmd)
                if (j + 1) % args.t == 0:
                    n += 1
                    if args.t > 1:
                        OUT.write("wait\n")
                    OUT.write("echo \"----------- %d ----------\"\n" %n)

            OUT.write("wait\n")
            OUT.write("echo \"----------- %d ----------\"\n" % (n+1))

        os.chmod(sub_shell_fname, stat.S_IRWXU)  # 0700











