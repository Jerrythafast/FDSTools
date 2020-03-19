#!/usr/bin/env python

#
# Copyright (C) 2020 Jerry Hoogenboom
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
from setuptools.extension import Extension
import fdstools as distmeta

import sys
if sys.hexversion >= 0x03000000:
    sys.stderr.write("error: This is FDSTools v%s, which is only compatible with Python2. "
                     "Please check FDSTools.nl to find out how to install FDSTools v2 on Python3, "
                     "or just try running 'pip3 install -U fdstools' to get the latest version." % distmeta.__version__)
    sys.exit(1)

x = setup(
    name="fdstools",
    packages=find_packages(),
    ext_modules=[
        Extension('fdstools.sg_align',
            sources=['fdstools/sg_align.c'],
            extra_compile_args=['-O3'])],
    package_data={
        "fdstools": ["vis/*.*", "vis/*/*"]
    },
    version=distmeta.__version__,
    install_requires=["numpy<1.17"],
    python_requires=">=2.7.9, <3",
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