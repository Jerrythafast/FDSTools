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
Mark all sequences that are not in another list of sequences.

Adds a new column 'new_allele' to the input data.  An asterisk (*) will
be placed in this column for any sequence that does not occur in the
provided list of known sequences.
"""
import argparse, sys

from ..lib import ensure_sequence_format, add_sequence_format_args, \
                  get_column_ids, add_input_output_args, SEQ_SPECIAL_VALUES, \
                  get_input_output_files

__version__ = "1.0.0"


def read_known(infile, library, default_marker=None):
    """Read sequence list from infile, return {marker: [sequences]}"""
    data = {}

    # Get column numbers.  The marker name column is optional.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_sequence = get_column_ids(column_names, "sequence")
    colid_marker = get_column_ids(column_names, "marker", optional=True)

    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        marker = line[colid_marker] if colid_marker is not None \
            else default_marker
        if default_marker is not None and marker != default_marker:
            continue
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
                                          library=library, marker=marker)
        if marker in data:
            data[marker].append(sequence)
        else:
            data[marker] = [sequence]

    return data
#read_known


def find_new(infile, outfile, known, library):
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_marker, colid_sequence = get_column_ids(column_names, "marker",
        "sequence")
    column_names.insert(colid_sequence + 1, "new_allele")

    outfile.write("\t".join(column_names) + "\n")
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        marker = cols[colid_marker]
        cols.insert(colid_sequence + 1, "" if marker in known and
            ensure_sequence_format(cols[colid_sequence], "raw",
            library=library, marker=marker) in known[marker] else "*")
        outfile.write("\t".join(cols) + "\n")
#find_new


def add_arguments(parser):
    parser.add_argument('known', metavar="KNOWN",
        type=argparse.FileType('r'),
        help="file containing a list of known allelic sequences")
    add_input_output_args(parser, True, True, False)
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    add_sequence_format_args(parser, "raw", True)
#add_arguments


def run(args):
    gen = get_input_output_files(args, True, True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")

    # Read list of known sequences once.
    known = read_known(args.known, args.library, args.marker)

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError(
                "multiple input files for sample '%s' specified " % tag)
        infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "r")
        find_new(infile, outfile, known, args.library)
        if infile != sys.stdin:
            infile.close()
#run
