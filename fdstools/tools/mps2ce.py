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
Convert sample data file to GeneMapper/GeneMarker format,
so that MPS data can be used with tools developed for CE data.
"""
from ..lib.cli import add_batch_input_output_args, add_sequence_format_args, iterate_samples
from ..lib.io import parse_flags, read_sample_data_file, std_or_open
from ..lib.seq import ensure_sequence_format

__version__ = "1.0.0"

# Default values for parameters are specified below.
_DEF_SEPARATOR = "\t"
_DEF_SEQ_FORMAT = "ce"


def add_missing_columns(sep, columns_wanted, columns_present):
    missing_columns = columns_wanted - columns_present

    if missing_columns < 0:
        raise ValueError("max_alleles (" + str(columns_wanted) + ") "
                         "should not be lower than number of current alleles"
                         " (" + str(columns_present) + ").")

    return sep * missing_columns
# add_missing_columns


def write_header(outfile, sep, max_alleles, pair_columns_by_allele):
    outfile.write("Sample Name" + sep + "Marker" + sep)

    if pair_columns_by_allele:
        outfile.write(sep.join(["Allele %i%sHeight %i" % (i, sep, i)
                      for i in range(1, max_alleles + 1)]))
    else:
        outfile.write(sep.join(
            ["Allele %i" % i for i in range(1, max_alleles + 1)] +
            ["Height %i" % i for i in range(1, max_alleles + 1)]))

    outfile.write("\n")
# write_header


def write_line(outfile, sep, sample, marker,
               pair_columns_by_allele, marker_data, max_alleles):
    outfile.write(sample + sep + marker + sep)

    if pair_columns_by_allele:
        outfile.write(sep.join(sep.join((k, str(round(v))))
                               for k, v in marker_data.items()))
        outfile.write(2 * add_missing_columns(sep, max_alleles, len(marker_data)))
    else:
        outfile.write(sep.join(marker_data.keys()))
        outfile.write(add_missing_columns(sep, max_alleles, len(marker_data)))
        outfile.write(sep)
        outfile.write(sep.join(str(round(v))
                               for v in marker_data.values()))
        outfile.write(add_missing_columns(sep, max_alleles, len(marker_data)))

    outfile.write("\n")
# write_line


def is_merge_output(samples):
    if len(samples) > 1:
        if samples[0][2] == samples[1][2]:
            return samples[0][2]
    return None
# is_merge_output


def is_ce_format(seqformat):
    if seqformat == "ce":
        return True, "allelename"
    else:
        return False, seqformat
# is_ce_format


def read_tssv(infile):
    with std_or_open(infile) as infile:
        return read_sample_data_file(
            infile,
            drop_special_seq=True,
            after_correction=True,
            combine_strands=True,
            extra_columns={"flags": True})
# read_tssv


def reformat_data(library, data, ce_format, seqformat, remove_non_str_markers, alleles_only):
    """
    reformat_data returns [{Marker: {Allele: Height}}, max_alleles]
    """
    reformatted_data = {}
    max_alleles = 0

    for (marker, allele), (height, extra_cols) in data.items():
        if remove_non_str_markers:
            if library is None:
                raise ValueError("The -r/--remove-non-str-markers setting is turned on, "
                                 "the detection of SNPs requires a library file")
            elif not library.get_range(marker).library:
                continue

        if alleles_only:
            if "flags" not in extra_cols:
                raise ValueError("Not all input files contain the 'flags' "
                                 "column, please turn off -a/--alleles-only "
                                 "or add the flags column to all input files.")
            if "allele" not in parse_flags(extra_cols["flags"]):
                continue

        allele = ensure_sequence_format(allele, seqformat, library=library, marker=marker)

        if ce_format and allele.startswith("CE"):
            allele = allele.split("_")[0][2:]

        if marker not in reformatted_data:
            reformatted_data[marker] = {allele: height}
        elif allele not in reformatted_data[marker]:
            reformatted_data[marker][allele] = height
        else:
            reformatted_data[marker][allele] += height

        if len(reformatted_data[marker]) > max_alleles:
            max_alleles = len(reformatted_data[marker])

    return reformatted_data, max_alleles
# reformat_data


def write_ce_file(data, outfile, sep, max_alleles,
                  pair_columns_by_allele):
    with std_or_open(outfile, "wt") as outfile:
        write_header(outfile, sep, max_alleles, pair_columns_by_allele)
        for tag in data:
            for marker in data[tag]:
                write_line(outfile, sep, tag, marker, pair_columns_by_allele,
                           data[tag][marker], max_alleles)
# write_ce_file


def mps2ce_tool(samples, library, sep, alleles_only,
                remove_non_str_markers, seqformat, pair_columns_by_allele):

    merge_output = is_merge_output(samples)
    ce_format, seqformat = is_ce_format(seqformat)

    outfile_data = {}  # {tag: {marker: data}}}
    max_alleles = 0
    for tag, infile, outfile in samples:
        tssv_data = read_tssv(infile)
        ce_data, sample_max_alleles = reformat_data(
            library, tssv_data, ce_format, seqformat, remove_non_str_markers, alleles_only)
        outfile_data[tag] = ce_data

        if merge_output is None:
            write_ce_file(outfile_data, outfile, sep,
                          sample_max_alleles, pair_columns_by_allele)
            outfile_data = {}
        elif sample_max_alleles > max_alleles:
            max_alleles = sample_max_alleles

    if merge_output is not None:
        write_ce_file(outfile_data, merge_output, sep,
                      max_alleles, pair_columns_by_allele)
# mps2ce_tool


def add_arguments(parser):
    add_batch_input_output_args(parser)

    output_group = parser.add_argument_group("output file format options")
    output_group.add_argument("-s", "--separator", metavar="SEPARATOR",
                              default=_DEF_SEPARATOR,
                              help="delimiter used to separate the columns in the "
                              "output file (default: %(default)s)")
    output_group.add_argument("-p", "--pair-columns-by-allele",
                              action="store_true",
                              help="by default, all Height columns come after "
                                   "all Allele columns; specify this option to "
                                   "place each Height column directly after "
                                   "the corresponding Allele column")

    filter_group = parser.add_argument_group("filtering options")
    filter_group.add_argument("-a", "--alleles-only", action="store_true",
                              help="if specified, only sequences flagged as "
                                   "'allele' are included in the output")
    filter_group.add_argument("-r", "--remove-non-str-markers", action="store_true",
                              help="if specified, non-STR markers are excluded from "
                                   "the output")

    add_sequence_format_args(parser, default_format=_DEF_SEQ_FORMAT, ce_format=True)
# add_arguments


def run(args):
    samples = list(iterate_samples(args, open_files=False))
    mps2ce_tool(samples, args.library, args.separator, args.alleles_only,
                args.remove_non_str_markers, args.sequence_format,
                args.pair_columns_by_allele)
# run
