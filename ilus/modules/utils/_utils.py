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
