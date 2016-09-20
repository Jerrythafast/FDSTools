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
import argparse
import sys

from ..lib import load_profiles_new, ensure_sequence_format, glob_path,\
                  add_sequence_format_args

__version__ = "1.0.2"


def merge_profiles(infiles, outfile, seqformat, library):
    outfile.write("\t".join(
        ["marker", "allele", "sequence", "fmean", "rmean", "tool"]) + "\n")
    seen = {}
    for infile in infiles:
        if infile == "-":
            profiles = load_profiles_new(sys.stdin, library)
        else:
            with open(infile, "r") as handle:
                profiles = load_profiles_new(handle, library)
        for marker in profiles:
            for allele in profiles[marker]:
                for sequence in profiles[marker][allele]:
                    if marker not in seen:
                        seen[marker] = {}
                    if allele not in seen[marker]:
                        seen[marker][allele] = set()
                    elif sequence in seen[marker][allele]:
                        continue
                    outfile.write("\t".join([marker] + [
                        ensure_sequence_format(seq, seqformat, library=library,
                            marker=marker) for seq in (allele, sequence)] +
                        map(str, (
                            profiles[marker][allele][sequence]["forward"],
                            profiles[marker][allele][sequence]["reverse"])) +
                        [profiles[marker][allele][sequence]["tool"]]) + "\n")
                    seen[marker][allele].add(sequence)
#merge_profiles


def add_arguments(parser):
    parser.add_argument('infiles', nargs='+', metavar="FILE",
        help="files containing the background noise profiles to combine; "
             "if a single file is given, it is merged with input from stdin; "
             "use '-' to use stdin as an explicit input source")
    outgroup = parser.add_argument_group("output file options")
    outgroup.add_argument('-o', '--output', dest="outfile", metavar="FILE",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="file to write output to (default: write to stdout)")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
#add_arguments


def run(args):
    infiles = [x for x in args.infiles for x in glob_path(x)]
    if len(infiles) < 2:
        if sys.stdin.isatty() or "-" in infiles:
            raise ValueError("please specify at least two input files")
        infiles.append("-")

    merge_profiles(infiles, args.outfile, args.sequence_format, args.library)
#run
