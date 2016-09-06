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
Link raw reads in a FastA or FastQ file to markers and count the number
of reads for each unique sequence.

This tool is basically a wrapper around the 'tssvl' program, offering
direct support for using FDSTools library files and allele name
generation.
"""
from __future__ import absolute_import  # Needed to import tssv package.
import sys
import math

# TSSV is only imported when actually running this tool.
#from tssv.tssv_lite import process_file, make_sequence_tables, \
#                           make_statistics_table, prepare_output_dir

from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_sequence_format_args, reverse_complement, PAT_SEQ_RAW,\
                  get_column_ids, ensure_sequence_format

__version__ = "1.0.2"


# Default values for parameters are specified below.

# Default maximum number of mismatches per nucleotide in the flanking
# sequences to allow.
# This value can be overridden by the -m command line option.
_DEF_MISMATCHES = 0.08

# Default penalty multiplier for insertions and deletions in the
# flanking sequences.
# This value can be overridden by the -n command line option.
_DEF_INDEL_SCORE = 1

# Default minimum number of reads to consider.
# This value can be overridden by the -a command line option.
_DEF_MINIMUM = 1


def convert_library(library, threshold):
    return {data[0]: (
        (data[1][0], reverse_complement(data[1][1])), (
            int(math.ceil(len(data[1][0]) * threshold)),
            int(math.ceil(len(data[1][1]) * threshold))))
                for data in (
                    (marker, library["flanks"][marker])
                        for marker in library["flanks"])}
#convert_library


def seq_pass_filt(sequence, reads, threshold, explen=None):
    """Return False if the sequence does not meet the criteria."""
    return (reads >= threshold and PAT_SEQ_RAW.match(sequence) is not None and
        (explen is None or explen[0] <= len(sequence) <= explen[1]))
#seq_pass_filt


def run_tssv_lite(infile, outfile, reportfile, is_fastq, library, seqformat,
                  threshold, minimum, aggregate_filtered,
                  missing_marker_action, dirname, indel_score):
    file_format = "fastq" if is_fastq else "fasta"
    tssv_library = convert_library(library, threshold)

    # Open output directory if we have one.
    if dirname:
        outfiles = prepare_output_dir(dirname, library["flanks"], file_format)
    else:
        outfiles = None

    total_reads, unrecognised, counters, sequences = process_file(
        infile, file_format, tssv_library, outfiles, indel_score)

    # Filter out sequences with low read counts and invalid bases now.
    if aggregate_filtered:
        aggregates = {}
        for marker in sequences:
            for sequence in sequences[marker]:
                if not seq_pass_filt(sequence,
                        sum(sequences[marker][sequence]), minimum,
                        library.get("expected_length", {}).get(marker)):
                    if marker not in aggregates:
                        aggregates[marker] = [0, 0]
                    aggregates[marker][0] += sequences[marker][sequence][0]
                    aggregates[marker][1] += sequences[marker][sequence][1]
    sequences = {marker:
        {sequence: sequences[marker][sequence]
            for sequence in sequences[marker]
            if seq_pass_filt(sequence,
                sum(sequences[marker][sequence]), minimum,
                library.get("expected_length", {}).get(marker))}
        for marker in sequences}

    # Add aggregate rows if the user requested so.
    if aggregate_filtered:
        for marker in aggregates:
            sequences[marker]["Other sequences"] = aggregates[marker]

    # Check presence of all markers.
    if missing_marker_action != "exclude":
        for marker in library["flanks"]:
            if not sequences[marker]:
                if missing_marker_action == "include":
                    sequences[marker]["No data"] = [0, 0]
                else:
                    raise ValueError("Marker %s was not detected!" % marker)

    column_names, tables = make_sequence_tables(sequences, 0)

    # Convert sequences to the desired format.
    colid_sequence = get_column_ids(column_names, "sequence")
    if seqformat != "raw":
        for marker in tables:
            for line in tables[marker]:
                line[colid_sequence] = ensure_sequence_format(
                    line[colid_sequence], seqformat, library=library,
                    marker=marker)

    # Write sequence tables.
    column_names = "\t".join(column_names)
    for marker in sorted(tables):
        tables[marker] = "\n".join(
            "\t".join(map(str, line)) for line in tables[marker])
        if outfiles:
            outfiles["markers"][marker]["sequences"].write(
                "\n".join((column_names, tables[marker])))
    tables = "\n".join(
        [column_names] + [tables[marker] for marker in sorted(tables)])
    if outfiles:
        outfiles["sequences"].write(tables)
    outfile.write(tables)
    outfile.write("\n")

    # Write statistics table.
    statistics = "\n".join((
        make_statistics_table(counters),
        "",  # Empty line.
        "total reads\t%i" % total_reads,
        "unrecognised reads\t%i" % unrecognised))
    if outfiles:
        outfiles["statistics"].write(statistics)
    reportfile.write(statistics)
    reportfile.write("\n")
#run_tssv_lite


def add_arguments(parser):
    add_sequence_format_args(parser, "raw", False, True)
    add_input_output_args(parser, True, False, True)
    parser.add_argument("-q", "--is-fastq", action="store_true",
        help="if specified, treat the input as a FASTQ file instead of FASTA")
    parser.add_argument("-D", "--dir",
        help="output directory for verbose output; when given, a subdirectory "
             "will be created for each marker, each with a separate "
             "sequences.csv file and a number of FASTA/FASTQ files containing "
             "unrecognised reads (unknown.fa), recognised reads "
             "(Marker/paired.fa), and reads that lack one of the flanks of a "
             "marker (Marker/noend.fa and Marker/nostart.fa)")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--mismatches", type=float,
        default=_DEF_MISMATCHES,
        help="number of mismatches per nucleotide to allow in flanking "
             "sequences (default: %(default)s)")
    filtergroup.add_argument("-n", "--indel-score", metavar="N",
        type=pos_int_arg, default=_DEF_INDEL_SCORE,
        help="insertions and deletions in the flanking sequences are "
             "penalised this number of times more heavily than mismatches "
             "(default: %(default)s)")
    filtergroup.add_argument("-a", "--minimum", metavar="N", type=pos_int_arg,
        default=_DEF_MINIMUM,
        help="report only sequences with this minimum number of reads "
             "(default: %(default)s)")
    filtergroup.add_argument("-A", "--aggregate-filtered", action="store_true",
        help="if specified, sequences that have been filtered (as per the "
             "-a/--minimum option, the expected_allele_length section in the "
             "library file, as well as all sequences with ambiguous bases) "
             "will be aggregated per marker and reported as 'Other sequences'")
    filtergroup.add_argument("-M", "--missing-marker-action", metavar="ACTION",
        choices=("include", "exclude", "halt"),
        default="include",
        help="action to take when no sequences are linked to a marker: one of "
             "%(choices)s (default: %(default)s)")
#add_arguments


def run(args):
    # Import TSSV now.
    global process_file, make_sequence_tables, make_statistics_table
    global prepare_output_dir
    try:
        from tssv.tssv_lite import process_file, make_sequence_tables, \
            make_statistics_table, prepare_output_dir
    except ImportError:
        raise ValueError(
            "This tool requires version 0.4.0 or later of the 'tssvl' program "
            "(TSSV-Lite) to be installed. Please download and install the "
            "latest version of TSSV from https://pypi.python.org/pypi/tssv.")

    files = get_input_output_files(args, True, False)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    infile = sys.stdin if files[0] == "-" else open(files[0], "r")
    run_tssv_lite(infile, files[1], args.report, args.is_fastq, args.library,
                  args.sequence_format, args.mismatches, args.minimum,
                  args.aggregate_filtered, args.missing_marker_action,
                  args.dir, args.indel_score)
    if infile != sys.stdin:
        infile.close()
#run