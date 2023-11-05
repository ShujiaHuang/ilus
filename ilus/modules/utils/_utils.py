"""Useful utilities functions for building analysis pipeline
"""
from pathlib import Path


def safe_makedir(dname, parents=True, exist_ok=False):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname

    num_tries = 0
    max_tries = 5
    while not Path(dname).exists():
        # Could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            Path(dname).mkdir(parents=parents, exist_ok=exist_ok)

        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1

    return dname


def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return Path(fname).is_file() and Path(fname).stat().st_size > 0
    except OSError:
        return False


def get_variant_calling_intervals(calling_interval_parameter):
    """get variant calling intervals into a bed format"""

    intervals = calling_interval_parameter
    if (type(calling_interval_parameter) is str) \
            and (Path(calling_interval_parameter).is_file()):
        # A file for recording interval
        interval_file = calling_interval_parameter
        with open(interval_file) as f:
            """Bed format:
            chr1	10001	207666
            chr1	257667	297968
            """
            # return the value to be a list of interval regions
            intervals = [line.strip().split()[:3] for line in f if not line.startswith("#")]

    elif type(calling_interval_parameter) is not list:
        raise ValueError(f"'variant_calling_interval' parameter could only be a file path or "
                         f"a list of chromosome id in the configure file (.yaml).\n")

    return intervals
