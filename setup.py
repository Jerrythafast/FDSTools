#!/usr/bin/env python

#
# Copyright (C) 2016 Jerry Hoogenboom
#
# This file is part of FDSTools, data analysis tools for Next
# Generation Sequencing of forensic DNA markers.
#
# FDSTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FDSTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FDSTools.  If not, see <http://www.gnu.org/licenses/>.
#

from setuptools import setup, find_packages

requires = ["numpy"]

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append("argparse")


# Disabling hard linking as a workaround for this bug:
# http://bugs.python.org/issue8876
from sys import hexversion as sys_hexversion
if sys_hexversion < 0x020709C1:  # The bug is fixed in 2.7.9rc1.
    import os
    del os.link


import fdstools as distmeta
x = setup(
    name="fdstools",
    packages=find_packages(),
    package_data={
        "fdstools": ["vis/*.*", "vis/*/*"]
    },
    version=distmeta.__version__,
    install_requires=requires,
    description="Forensic DNA Sequencing Tools",
    long_description=distmeta.__doc__,
    author="Jerry Hoogenboom",
    author_email="jerryhoogenboom@outlook.com",
    url="https://git.lumc.nl/jerryhoogenboom/fdstools/blob/master/README.rst",
    license="GPLv3+",
    platforms=["any"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "License :: OSI Approved :: GNU General Public License v3 or "
            "later (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords='bioinformatics forensics stutter NGS MPS DNA sequencing STR',
    entry_points={
        'console_scripts': [
            "fdstools=fdstools.fdstools:main"
        ]
    }
)