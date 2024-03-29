"""Setup file and install scripts.

Version 1.0.0 (May 23, 2020)
Copyright (c) 2020 Shujia Huang
"""
from argparse import Namespace
from pathlib import Path

try:
    from setuptools import setup, find_packages

    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages

DESCRIPTION = "ilus: A handy tools for generating WGS/WES analysis pipeline."
meta = Namespace(
    __DISTNAME__="ilus",
    __AUTHOR__="Shujia Huang",
    __AUTHOR_EMAIL__="huangshujia9@gmail.com",
    __URL__="https://github.com/ShujiaHuang/ilus",
    __LICENSE__="BSD (3-clause)",
    __DOWNLOAD_URL__="https://github.com/ShujiaHuang/ilus",
    __VERSION__="2.0.0",
)

if __name__ == "__main__":
    this_path = Path(__file__).resolve().parent
    long_description = this_path.joinpath("README.md")

    setup(name=meta.__DISTNAME__,
          author=meta.__AUTHOR__,
          author_email=meta.__AUTHOR_EMAIL__,
          maintainer=meta.__AUTHOR__,
          maintainer_email=meta.__AUTHOR_EMAIL__,
          description=DESCRIPTION,
          long_description=(open(long_description).read()),
          long_description_content_type="text/markdown",
          license=meta.__LICENSE__,
          url=meta.__URL__,
          version=meta.__VERSION__,
          download_url=meta.__DOWNLOAD_URL__,
          packages=find_packages(),
          python_requires=">=3.5",
          install_requires=[
              'PyYAML>=5.1.2',
          ],
          include_package_data=True,
          # scripts = [],
          entry_points={
              "console_scripts": [
                  'ilus = ilus.main:main'
              ]
          },
          classifiers=[
              'License :: OSI Approved :: BSD License',
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',
              'Programming Language :: Python :: 3.9',
              'Programming Language :: Python :: 3.10',
              'Programming Language :: Python :: 3.11',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Operating System :: POSIX',
              'Operating System :: POSIX :: Linux',
              'Operating System :: MacOS'],
          )
