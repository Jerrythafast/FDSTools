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

import setuptools, sys

with open("README.rst", "r") as fh:
    long_description = fh.read()

version = {}
with open("fdstools/__init__.py", "r") as fh:
    exec(fh.read(), version)

if sys.hexversion >= 0x03000000:
    sys.stderr.write("error: This is FDSTools v%s, which is only compatible with Python2. "
                     "Please check FDSTools.nl to find out how to install FDSTools v2 on Python3, "
                     "or just try running 'pip3 install -U fdstools' to get the latest version.\n"
                     % version["__version__"])
    sys.exit(1)

setuptools.setup(
    name="fdstools",
    version=version["__version__"],
    description="Forensic DNA Sequencing Tools",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://fdstools.nl",
    project_urls={
        "Bug Tracker": "https://github.com/Jerrythafast/FDSTools/issues",
        "Source Code": "https://github.com/Jerrythafast/FDSTools",
    },
    author="Jerry Hoogenboom",
    author_email="jerryhoogenboom@outlook.com",
    license="GPLv3+",
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords='bioinformatics forensics stutter NGS MPS DNA sequencing STR',
    packages=setuptools.find_packages(),
    ext_modules=[
        setuptools.extension.Extension('fdstools.sg_align',
            sources=['fdstools/sg_align.c'],
            extra_compile_args=['-O3'])],
    package_data={
        "fdstools": ["vis/*.*", "vis/*/*"]
    },
    install_requires=["numpy<1.17"],
    python_requires=">=2.7.9, <3",
    entry_points={
        'console_scripts': [
            "fdstools=fdstools.fdstools:main"
        ]
    }
)
