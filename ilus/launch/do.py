"""Centralize running of external commands
"""
import collections
import os
import subprocess
from datetime import datetime

from ilus.log import logger, logger_cmd, logger_stdout


def run(cmd, descr=None, checks=None, log_stdout=False, env=None):
    """Run the provided command, logging details and checking for error
    """

    start_time = datetime.now()
    if descr:
        logger.info(descr)

    try:
        logger_cmd.debug(" ".join(map(str, cmd))
                         if not isinstance(cmd, basestring) else cmd)
        _do_run(cmd, checks, log_stdout=log_stdout, env=env)
    except:
        logger.exception()
        raise

    # return the elsape time
    return datetime.now() - start_time


def find_bash():
    for test_bash in [find_cmd("bash"), "/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash

    raise IOError("Error: Could not find bash in any standard location. Needed for unix pipes")


def find_cmd(cmd):
    try:
        return subprocess.check_output(["which", cmd]).strip()
    except subprocess.CalledProcessError:
        return None


def _normalize_cmd_args(cmd):
    """Normalize subprocess arguments to handle list commands, string and pipes.
    Piped commands set pipefail and require use of bash to help with debugging
    intermediate errors.
    """
    if isinstance(cmd, basestring):
        # check for standard or anonymous named pipes
        if cmd.find(" | ") > 0 or cmd.find(">(") or cmd.find("<("):
            return "set -o pipefail; " + cmd, True, find_bash()
        else:
            return cmd, True, None
    else:
        return map(str, cmd), False, None


def _do_run(cmd, checks, log_stdout=False, env=None):
    """Perform running and check results, raising errors for issues.
    """
    cmd, is_shell_args, executable_args = _normalize_cmd_args(cmd)
    p = subprocess.Popen(
        cmd,
        shell=is_shell_args,
        executable=executable_args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        close_fds=True,
        env=env
    )

    debug_stdout = collections.deque(maxlen=100)
    while 1:
        line = p.stdout.readline()
        if line.rstrip():
            debug_stdout.append(line)
            if log_stdout:
                logger_stdout.debug(line.rstrip())
            else:
                logger.debug(line.rstrip())

        exitcode = p.poll()
        if exitcode is not None:
            for line in p.stdout:
                debug_stdout.append(line)

            if exitcode is not None and exitcode != 0:

                error_msg = " ".join(cmd) if not isinstance(cmd, basestring) else cmd
                error_msg += "\n"
                error_msg += "".join(debug_stdout)

                p.communicate()
                p.stdout.close()

                raise subprocess.CalledProcessError(exitcode, error_msg)

            else:

                # output logo information
                print ("".join(debug_stdout))
                break

    p.communicate()
    p.stdout.close()

    # Check for problems not identified by shell return codes
    if checks:
        for check in checks:
            if not check():
                raise IOError("External command failed")