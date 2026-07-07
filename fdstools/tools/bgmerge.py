#!/usr/bin/env python3

#
# Copyright (C) 2026 Netherlands Forensic Institute
#
# This file is part of FDSTools, data analysis tools for Massively
# Parallel Sequencing of forensic DNA markers.
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

"""
Merge multiple files containing background noise profiles.

Background noise profiles are merged in the order in which they are
specified.  If multple files specify a different value for the same
allele and sequence, the value of the first file is used.

It is convenient to pipe the output of bgpredict and/or bgestimate
into bgmerge to merge that with an existing file containing background
profiles.  Specify '-' as one of the input files to read from stdin
(i.e., read input from a pipe).  If only one input file is specified,
'-' is implicitly used as the second input file.  Note that as a result,
in case of conflicting values, the value in the specified input file
will take precedence over the value in the data that was piped in.

Example: fdstools bgpredict ... | fdstools bgmerge old.txt > out.txt
"""
import sys

from itertools import chain

from ..lib.cli import add_sequence_format_args, glob_path
from ..lib.io import std_or_open
from ..lib.noise import load_profiles

__version__ = "1.1.0"


def merge_profiles(infiles, outfile, library):
    merged_profiles = {}
    for infile in infiles:
        with std_or_open(infile) as handle:
            profiles = load_profiles(handle, library)
        for marker, markerprofile in profiles.items():
            if marker not in merged_profiles:
                merged_profiles[marker] = markerprofile
                continue
            mergedmarkerprofile = merged_profiles[marker]
            for allele, alleleprofile in markerprofile.items():
                if allele not in mergedmarkerprofile:
                    mergedmarkerprofile[allele] = alleleprofile
                    continue
                mergedalleleprofile = mergedmarkerprofile[allele]
                for sequence, sequenceprofile in alleleprofile.items():
                    if sequence not in mergedalleleprofile:
                        mergedalleleprofile[sequence] = sequenceprofile
                        continue
                    for strand in ("forward", "reverse", "total"):
                        if not mergedalleleprofile[sequence][strand]:
                            mergedalleleprofile[sequence][strand] = sequenceprofile[strand]
                            mergedalleleprofile[sequence]["tools"].update(sequenceprofile["tools"])

    outfile.write("\t".join(
        ("marker", "allele", "sequence", "fmean", "rmean", "tmean", "tools")) + "\n")
    for marker, markerprofile in merged_profiles.items():
        for allele, alleleprofile in markerprofile.items():
            for sequence, sequenceprofile in alleleprofile.items():
                outfile.write("\t".join(map(str, (
                    marker, allele, sequence,
                    sequenceprofile["forward"],
                    sequenceprofile["reverse"],
                    sequenceprofile["total"],
                    ",".join(sorted(sequenceprofile["tools"]))))) + "\n")
#merge_profiles


def add_arguments(parser):
    parser.add_argument("infiles", nargs="+", metavar="FILE",
        help="files containing the background noise profiles to combine; "
             "if a single file is given, it is merged with input from stdin; "
             "use '-' to use stdin as an explicit input source")
    parser.add_argument("-o", "--output", dest="outfile", metavar="FILE", default="-",
        help="file to write output to (default: write to stdout)")
    add_sequence_format_args(parser, default_format="raw", force=True)
#add_arguments


def run(args):
    infiles = list(chain.from_iterable(map(glob_path, args.infiles)))
    if len(infiles) < 2:
        if sys.stdin.isatty() or "-" in infiles:
            raise ValueError("please specify at least two input files")
        infiles.append("-")

    with std_or_open(args.outfile, "wt") as outfile:
        merge_profiles(infiles, outfile, args.library)
#run
