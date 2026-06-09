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
Convert between raw sequences, TSSV-style sequences, and allele names.

FDSTools was built to be compatible with TSSV, which writes sequences of
known STR alleles in a shortened form referred to as 'TSSV-style
sequences'.  At the same time, FDSTools supports the creation of
human-readable allele names which are more suitable for display.

For example, the raw sequence
'AGCGTAAGATAGATAGATAGATAGATAGATACCTACCTACCTCTAGCT' might be rewritten as
the TSSV-style sequence 'AGCGTA(1)AGAT(6)ACCT(3)CTAGCT(1)', or as the
allele name 'CE9_AGAT[6]ACCT[3]'.

Seqconvert can be used to explicitly convert all sequences in a file to
the same output format.  Conversions are done using a library file, see
the help text of the library tool for details.

You can specify multiple input files using the -i/--input option.  This
is especially useful when generating allele names for many samples that
have many sequences in common.  To call the variants in the allele
names, FDSTools needs to do sequence alignments which can be rather
slow.  When generating allele names for many input files at once, the
results of the alignments are cached which may give a significant
speed-up compared to generating allele names for each sample separately.

Seqconvert can also be used with two different library files to rewrite
the allele names or TSSV-style sequences after a library update.
Currently, the only limitation to this is that the ending position of
the left flank and the starting position of the right flank must be the
same.

Note that FDSTools makes no assumptions about the sequence format in its
input files; instead it automatically performs any required conversions
while running any tool.  Explicitly running seqconvert is never a
necessity; use this tool for your own convenience.
"""
import logging

from functools import cache

from ..lib.cli import add_batch_input_output_args, iterate_samples, library_arg
from ..lib.io import get_column_ids
from ..lib.library import BUILTIN_NAMES, CUSTOM_REPORTEDRANGE_OPTION
from ..lib.seq import SEQ_SPECIAL_VALUES, reverse_complement, ensure_sequence_format, \
                      sequence_range_adjust

__version__ = "1.1.0"

LOG = logging.getLogger(__name__)


# Default values for parameters are specified below.

# Default name of the column that contains the marker name.
# This value can be overridden by the -m command line option.
_DEF_COLNAME_MARKER = "marker"

# Default name of the column that contains the sequence.
# This value can be overridden by the -a command line option.
_DEF_COLNAME_SEQUENCE = "sequence"


@cache
def should_adjust_range(marker, old, new):
    if old.get_option(CUSTOM_REPORTEDRANGE_OPTION) or new.get_option(CUSTOM_REPORTEDRANGE_OPTION):
        # The user specified custom configuration for this marker.
        LOG.debug("Not adjusting reported range for marker %s because it is based on customized "
                  "configuration", marker)
        return False
    if len(old.location) != 3 or len(new.location) != 3:
        # The marker must be and remain mapped to one consecutive piece
        # of reference sequence.
        LOG.debug("Not adjusting reported range for marker %s because it is not mapped to one "
                  "contiguous region of the GRCh38 reference sequence", marker)
        return False
    if old.location == new.location:
        # No adjustment is necessary.
        return False
    chr, orig_start, orig_end = old.location
    new_chr, new_start, new_end = new.location
    if new_chr != chr or min(orig_end, new_end) < max(orig_start, new_start):
        # The old and new range should overlap.
        LOG.debug("Not adjusting reported range for marker %s because it gets remapped to a "
                  "non-overlapping region of the GRCh38 reference sequence", marker)
        return False
    if orig_end < new_end or orig_start > new_start:
        LOG.warning("Sequences for marker %s will be extended with nucleotides from the GRCh38 "
                    "reference sequence", marker)
    return True
#should_adjust_range


def convert_sequences(infile, outfile, to_format, library=None, fixed_marker=None,
                      colname_marker=_DEF_COLNAME_MARKER, colname_sequence=_DEF_COLNAME_SEQUENCE,
                      colname_sequence_out=None, library2=None, revcomp_markers=[]):
    if colname_sequence_out is None:
        colname_sequence_out = colname_sequence
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_sequence = get_column_ids(column_names, colname_sequence)
    if library is None:
        fixed_marker = ""  # Don't need marker names without library.
    if fixed_marker is None:
        colid_marker = get_column_ids(column_names, colname_marker)
    try:
        colid_sequence_out = get_column_ids(column_names, colname_sequence_out)
    except:
        column_names.append(colname_sequence_out)
        colid_sequence_out = -1

    outfile.write("\t".join(column_names) + "\n")
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        if colid_sequence_out == -1:
            line.append("")
        marker = line[colid_marker] if fixed_marker is None else fixed_marker
        seq = line[colid_sequence]

        if seq not in SEQ_SPECIAL_VALUES:
            from_format = None

            if marker in revcomp_markers or library2 != library:
                # Need to convert input to raw sequences first.
                seq = ensure_sequence_format(seq, "raw", marker=marker, library=library)
                from_format = "raw"

            if marker in revcomp_markers:
                # Reverse-complement the sequence, no range adjustment.
                seq = reverse_complement(seq)
            elif library2 != library and should_adjust_range(
                    marker, library.get_range(marker), library2.get_range(marker)):
                # Modify sequence to adjust its reported range.
                # We can efficiently convert to allele name if requested.
                if to_format == "allelename":
                    output_name = True
                    from_format = "allelename"
                else:
                    output_name = False
                seq = sequence_range_adjust(seq, library.get_range(marker),
                    library2.get_range(marker), library2, output_name=output_name)

            seq = ensure_sequence_format(
                seq, to_format, from_format=from_format, marker=marker, library=library2)

        line[colid_sequence_out] = seq
        outfile.write("\t".join(line) + "\n")
#convert_sequences


def add_arguments(parser):
    parser.add_argument("sequence-format", metavar="FORMAT", choices=("raw", "tssv", "allelename"),
        help="the format to convert to: one of %(choices)s")
    add_batch_input_output_args(parser)
    parser.add_argument("-m", "--marker-column", metavar="COLNAME", default=_DEF_COLNAME_MARKER,
        help="name of the column that contains the marker name (default: '%(default)s')")
    parser.add_argument("-a", "--allele-column", metavar="COLNAME", default=_DEF_COLNAME_SEQUENCE,
        help="name of the column that contains the sequence (default: '%(default)s')")
    parser.add_argument("-c", "--output-column", metavar="COLNAME",
        help="name of the column to write the output to (default: same as -a/--allele-column)")
    parser.add_argument("-M", "--marker", metavar="MARKER",
        help="assume the specified marker for all sequences")
    parser.add_argument("-l", "--library", metavar="LIBRARY", type=library_arg,
        help="library file with marker definitions; custom file or built-in: '%s'" %
             "', '".join(BUILTIN_NAMES))
    parser.add_argument("-L", "--library2", metavar="LIBRARY", type=library_arg,
        help="second library file to use for output; if specified, allele "
             "names can be conveniently updated to fit this new library file")
    parser.add_argument("-r", "--reverse-complement", metavar="MARKER", nargs="+", default=[],
        help="to be used together with -L/--library2; specify the markers for "
             "which the sequences are reverse-complemented in the new library")
#add_arguments


def run(args):
    library = args.library if args.library is not None else args.library2
    library2 = args.library2 if args.library2 is not None else library

    for tag, infile, outfile in iterate_samples(args):
        convert_sequences(infile, outfile, getattr(args, "sequence-format"), library,
                          args.marker, args.marker_column, args.allele_column,
                          args.output_column, library2, args.reverse_complement)
#run
