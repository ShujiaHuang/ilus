"""Setting configurations from .yaml files and environment 
variables for pipeline.
"""
import sys
import os
import math


class CmdNotFound(Exception):
    pass


def expand_path(path):
    """ Combines os.path.expandvars with replacing ~ with $HOME.
    """
    try:
        return os.path.expandvars(path.replace("~", "$HOME"))
    except AttributeError:
        return path


def get_program(name, config, ptype="cmd", default=None):
    """Retrieve program information from the configuration.

    This handles back compatible location specification in input
    YAML. The preferred location for program information is in
    `resources` but the older `program` tag is also supported.
    """
    # support taking in the data dictionary
    config = config.get("config", config)
    try:
        pconfig = config.get("resources", {})[name]
    except KeyError:
        pconfig = {}

    old_config = config.get("program", {}).get(name, None)
    if old_config:
        for key in ["dir", "cmd"]:
            if not key in pconfig:
                pconfig[key] = old_config

    if ptype == "cmd":
        return _get_program_cmd(name, pconfig, config, default)
    elif ptype == "dir":
        return _get_program_dir(name, pconfig)
    else:
        raise ValueError("Don't understand program type: %s" % ptype)


def _get_check_program_cmd(fn):
    def wrap(name, pconfig, config, default):
        is_ok = lambda f: os.path.isfile(f) and os.access(f, os.X_OK)
        bcbio_system = config.get("bcbio_system", None)
        if bcbio_system:
            system_bcbio_path = os.path.join(os.path.dirname(bcbio_system),
                                             os.pardir, "anaconda", "bin", name)
            if is_ok(system_bcbio_path):
                return system_bcbio_path

        # support bioconda installed programs
        if is_ok(os.path.join(os.path.dirname(sys.executable), name)):
            return (os.path.join(os.path.dirname(sys.executable), name))
        # find system bioconda installed programs if using private code install
        program = expand_path(fn(name, pconfig, config, default))
        if is_ok(program):
            return program
        # search the PATH now
        for adir in os.environ['PATH'].split(":"):
            if is_ok(os.path.join(adir, program)):
                return os.path.join(adir, program)

        raise CmdNotFound(" ".join(map(repr, (fn.func_name, name, pconfig, default))))
    return wrap


@_get_check_program_cmd
def _get_program_cmd(name, pconfig, config, default):
    """Retrieve commandline of a program.
    """
    if pconfig is None:
        return name
    elif isinstance(pconfig, basestring):
        return pconfig
    elif "cmd" in pconfig:
        return pconfig["cmd"]
    elif default is not None:
        return default
    else:
        return name


def _get_program_dir(name, config):
    """Retrieve directory for a program (local installs/java jars).
    """
    if config is None:
        raise ValueError("Could not find directory in config for %s" % name)
    elif isinstance(config, basestring):
        return config
    elif "dir" in config:
        return expand_path(config["dir"])
    else:
        raise ValueError("Could not find directory in config for %s" % name)


def adjust_memory(val, magnitude, direction="increase", out_modifier="", maximum=None):
    """Adjust memory based on number of cores utilized.
    """
    modifier = val[-1:]
    amount = float(val[:-1])

    if direction == "decrease":
        new_amount = amount / float(magnitude)
        # dealing with a specifier like 1G, need to scale to Mb
        if (new_amount < 1 or (out_modifier.upper().startswith("M") and
                                   modifier.upper().startswith("G"))):

            if modifier.upper().startswith("G"):
                new_amount = (amount * 1024) / magnitude
                modifier = "M" + modifier[1:]

            else:
                raise ValueError("Unexpected decrease in memory: "
                                 "%s by %s" % (val, magnitude))

        amount = int(new_amount)

    elif direction == "increase" and magnitude > 1:
        # for increases with multiple cores, leave small percentage of
        # memory for system to maintain process running resource and
        # avoid OOM killers
        adjuster = 0.90
        amount = int(math.ceil(amount * (adjuster * magnitude)))

    else:
        pass

    if out_modifier.upper().startswith("G") and modifier.upper().startswith("M"):
        modifier = out_modifier
        amount = int(math.floor(amount / 1024.0))

    if out_modifier.upper().startswith("M") and modifier.upper().startswith("G"):
        modifier = out_modifier
        modifier = int(amount * 1024)

    if maximum:
        max_modifier = maximum[-1]
        max_amount = float(maximum[:-1])
        if modifier.upper() == "G" and max_modifier.upper() == "M":
            max_amount = max_amount / 1024.0

        elif modifier.upper() == "M" and max_modifier.upper() == "G":
            max_amount = max_amount * 1024.0

        amount = min([amount, max_amount])

    return "{amount}{modifier}".format(amount=int(math.floor(amount)), modifier=modifier)
