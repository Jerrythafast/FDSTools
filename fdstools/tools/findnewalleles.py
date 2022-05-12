#!/usr/bin/env python3

#
# Copyright (C) 2022 Jerry Hoogenboom
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
Mark all sequences that are not in another list of sequences.

If not present, a new column 'flags' is added to the output.  Any sequence that
does not occur in the provided list of known sequences is flagged 'novel'.
"""
import argparse
import sys

from errno import EPIPE

from ..lib.cli import add_sequence_format_args, add_input_output_args, get_input_output_files
from ..lib.io import get_column_ids, parse_flags
from ..lib.seq import SEQ_SPECIAL_VALUES, ensure_sequence_format

__version__ = "1.1.1"


def read_known(infile, library):
    """Read sequence list from infile, return {marker: [sequences]}."""
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return {}  # Empty file.
    colid_marker, colid_sequence = get_column_ids(column_names, "marker", "sequence")

    data = {}
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        marker = line[colid_marker]
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
                                          library=library, marker=marker)
        if marker in data:
            data[marker].append(sequence)
        else:
            data[marker] = [sequence]

    return data
#read_known


def find_new(infile, outfile, known, library, remove_allele_flags):
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return  # Empty file.
    colid_marker, colid_sequence = get_column_ids(column_names, "marker", "sequence")
    colid_flags = get_column_ids(column_names, "flags", optional=True)
    if colid_flags is None:
        column_names.append("flags")
        colid_flags = -1

    outfile.write("\t".join(column_names) + "\n")
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if colid_flags == -1:
            line.append([])
        else:
            line[colid_flags] = parse_flags(line[colid_flags])
        if line[colid_sequence] not in SEQ_SPECIAL_VALUES:
            marker = line[colid_marker]
            if marker in known and ensure_sequence_format(line[colid_sequence], "raw",
                    library=library, marker=marker) not in known[marker]:
                line[colid_flags].append("novel")
                if remove_allele_flags and "allele" in line[colid_flags]:
                    line[colid_flags].remove("allele")
        line[colid_flags] = ",".join(line[colid_flags])
        outfile.write("\t".join(line) + "\n")
#find_new


def add_arguments(parser):
    parser.add_argument("known", metavar="KNOWN", type=argparse.FileType("rt", encoding="UTF-8"),
        help="file containing a list of known allelic sequences")
    parser.add_argument("-r", "--remove-allele-flags", action="store_true",
        help="remove the 'allele' flag from the alleles that are marked 'novel'")
    add_input_output_args(parser, single_in=True, batch_support=True, report_out=False)
    add_sequence_format_args(parser, default_format="raw", force=True)
#add_arguments


def run(args):
    gen = get_input_output_files(args, single_in=True, batch_support=True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output of another program")

    # Read list of known sequences once.
    known = read_known(args.known, args.library)

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError("multiple input files for sample '%s' specified " % tag)
        try:
            infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "rt", encoding="UTF-8")
            find_new(infile, outfile, known, args.library, args.remove_allele_flags)
            if infile != sys.stdin:
                infile.close()
        except IOError as e:
            if e.errno == EPIPE:
                continue
            raise
#run
