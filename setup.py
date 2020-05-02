"""Setup file and install script for NGS data analysis.

Version 0.1.0 (April 23, 2020)
Copyright (c) 2020 Shujia Huang
"""
import os

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


DESCRIPTION = "ilus: A python package for NGS data analysis."
DISTNAME = 'ilus'
MAINTAINER = 'Shujia Huang'
MAINTAINER_EMAIL = 'hshujia@qq.com'
URL = 'https://github.com/ShujiaHuang/ilus'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ShujiaHuang/ilus'
VERSION = "0.1.0"

ROOT_DIR = os.path.split(os.path.realpath(__file__))[0]

if __name__ == "__main__":
    setup(name=DISTNAME,
          author=MAINTAINER,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=(open(ROOT_DIR + "/README.rst").read()),
          license=LICENSE,
          url=URL,
          download_url=DOWNLOAD_URL,
          packages=find_packages(),
          install_requires=[
              'PyYAML>=5.1.2',
              'Logbook>=1.5.0',
          ],
          version=VERSION,
          include_package_data=True,
          # scripts = [],
          # entry_points={
          #     "console_scripts": []
          # },
          classifiers=[
             'Intended Audience :: Science/Research',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Python :: 3.7',
             'License :: OSI Approved :: BSD License',
             'Topic :: Scientific/Engineering :: Bio-Informatics',
             'Topic :: Scientific/Engineering :: Tools',
             'Topic :: Multimedia :: WGS',
             'Operating System :: POSIX',
             'Operating System :: Linux/Unix',
             'Operating System :: MacOS'],
          )
