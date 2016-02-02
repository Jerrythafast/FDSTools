#!/usr/bin/env python
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
    url="https://git.lumc.nl/jhoogenboom/fdstools/blob/master/README.rst",
    license="LGPLv3+",
    platforms=["any"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or "
            "later (LGPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords='bioinformatics forensics stutter NGS sequencing STR',
    entry_points={
        'console_scripts': [
            "fdstools=fdstools.fdstools:main"
        ]
    }
)