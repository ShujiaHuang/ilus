"""Setup file and install scripts.

Version 1.0.0 (May 23, 2020)
Copyright (c) 2020 Shujia Huang
"""
import os
from argparse import Namespace

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


DESCRIPTION = "ilus: A handy tools for generating WGS/WES analysis pipeline."
meta = Namespace(
    __DISTNAME__     = "ilus",
    __AUTHOR__       = "Shujia Huang",
    __AUTHOR_EMAIL__ = "huangshujia9@gmail.com",
    __URL__          = "https://github.com/ShujiaHuang/ilus",
    __LICENSE__      = "BSD (3-clause)",
    __DOWNLOAD_URL__ = "https://github.com/ShujiaHuang/ilus",
    __VERSION__      = "1.1.0",
)


if __name__ == "__main__":

    THIS_PATH = os.path.abspath(os.path.dirname(__file__))
    long_description = os.path.join(THIS_PATH, "README.md")

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
             'Intended Audience :: Science/Research',
             'Programming Language :: Python :: 3.7',
             'Programming Language :: Python :: 3.8',
             'Programming Language :: Python :: 3.9',
             'License :: OSI Approved :: BSD License',
             'Topic :: Scientific/Engineering :: Bio-Informatics',
             'Operating System :: POSIX',
             'Operating System :: POSIX :: Linux',
             'Operating System :: MacOS'],
          )
