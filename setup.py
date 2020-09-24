"""Setup file and install scripts.

Version 1.0.0 (May 23, 2020)
Copyright (c) 2020 Shujia Huang
"""
import os

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


DESCRIPTION = "ilus: A handy tools for generating WGS/WES analysis pipeline."
DISTNAME = 'ilus'
MAINTAINER = 'Shujia Huang'
MAINTAINER_EMAIL = 'huangshujia9@gmail.com'
URL = 'https://github.com/ShujiaHuang/ilus'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ShujiaHuang/ilus'
VERSION = "1.0.1"

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
          ],
          version=VERSION,
          include_package_data=True,
          # scripts = [],
          entry_points={
              "console_scripts": [
                  'ilus = ilus.main:main'
              ]
          },
          classifiers=[
             'Intended Audience :: Science/Research',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Python :: 3.7',
             'License :: OSI Approved :: BSD License',
             'Topic :: Scientific/Engineering :: Bio-Informatics',
             'Operating System :: POSIX',
             'Operating System :: POSIX :: Linux',
             'Operating System :: MacOS'],
          )
