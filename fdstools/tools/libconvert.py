#!/usr/bin/env python3

#
# Copyright (C) 2021 Jerry Hoogenboom
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
Convert a legacy TSSV library file (tab-separated) to the FDSTools library
format (ini-style).

This is a convenience tool for users migrating from the standalone
'TSSV' programme.  Use the 'library' tool if you wish to create a new,
empty FDSTools library file to start with.

Both FDSTools and the standalone 'TSSV' programme use a library file to
store the names, flanking (primer) sequences, and STR repeat structure
of forensic STR markers.  However, the TSSV library file format is not
well suited for non-STR markers and automatic generation of allele names.
FDSTools therefore employs a different (ini-style) library file format
that can store more details about the markers used.  The libconvert tool
can be used to convert old library files to the new format.

Please refer to the help of the 'library' tool for more information
about FDSTools library files.
"""
import argparse
import io
import os.path
import sys

from errno import EPIPE

from ..lib.library import PAT_STR_DEF, INI_COMMENT
from ..lib.seq import PAT_SEQ_RAW
from .library import make_empty_library_ini

__version__ = "1.2.0"


def convert_library(infile, outfile):
    markers = {}
    for line in infile:
        line = [x.strip() for x in line.rstrip("\r\n").split("\t")]
        if line == [""]:
            continue
        if len(line) < 4:
            raise ValueError(
                "Invalid library file: encountered line with %i columns, "
                "need at least 4" % len(line))
        marker = line[0]
        if PAT_SEQ_RAW.match(line[1]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" % (line[1], marker))
        if PAT_SEQ_RAW.match(line[2]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" % (line[2], marker))
        if PAT_STR_DEF.match(line[3]) is None:
            raise ValueError("STR definition '%s' of marker %s is invalid" % (line[3], marker))
        markers[marker] = line[1:4]

    outfile.write(INI_COMMENT.fill(
        "Lines beginning with a semicolon (;) are ignored by FDSTools.") + "\n\n")
    ini = make_empty_library_ini("str")

    # Enter flanking sequences and STR definitions.
    fmt = "%%-%is" % max(map(len, markers.keys() or [""]))
    for marker in sorted(markers):
        ini.set("flanks", fmt % marker, ", ".join(markers[marker][:2]))
        ini.set("repeat", fmt % marker, markers[marker][2])

    # Write INI file.
    ini.write(outfile)
#convert_library


def add_arguments(parser):
    parser.add_argument("infile", nargs="?", metavar="IN", default=sys.stdin,
        help="input library in the legacy TSSV format (default: read from stdin)")
    parser.add_argument("outfile", nargs="?", metavar="OUT",
        default=sys.stdout, type=argparse.FileType("tw", encoding="UTF-8"),
        help="the file to write the FDSTools library to (default: write to stdout)")
#add_arguments


def run(args):
    if args.infile == "-":
        args.infile = sys.stdin
    if args.infile != sys.stdin and args.outfile == sys.stdout and not os.path.exists(args.infile):
        # One filename given, and it does not exist.  Assume outfile.
        args.outfile = open(args.infile, "wt", encoding="UTF-8")
        args.infile = sys.stdin

    if args.infile != sys.stdin:
        # Open the specified input file.
        args.infile = open(args.infile, "rt", encoding="UTF-8")
    elif args.infile.isatty():
        # No input given.  Produce a default FDSTools library.
        args.infile = io.StringIO("")

    try:
        convert_library(args.infile, args.outfile)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#run
