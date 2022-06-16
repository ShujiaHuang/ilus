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
        sub_shell_fname = os.path.join(sub_shell_dir, "work.%d.sh" % (i+1))
        with open(sub_shell_fname, "w") as OUT:

            n = 0
            OUT.write("#!/bin/bash\n")
            for j, cmd in enumerate(commands[k:k+cmd_num_perjob]):
                OUT.write("%s &\n" % cmd if args.t > 1 else "%s\n" % cmd)
                if (j + 1) % args.t == 0:
                    n += 1
                    if args.t > 1:
                        OUT.write("wait\n")
                    OUT.write("echo \"----------- %d/%d ----------\"\n" % (n, part_num))

            if (j + 1) % args.t > 0:
                OUT.write("wait\n")
                OUT.write("echo \"----------- %d/%d ----------\"\n" % (n+1, part_num))

        os.chmod(sub_shell_fname, stat.S_IRWXU)  # 0700











