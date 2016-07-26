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
Mark potential stutter products by assuming a fixed maximum percentage
of stutter product vs the parent sequence.

Stuttermark adds a new column (named 'annotation' by default) to the
output.  The new column contains 'STUTTER' for possible stutter
products, or 'ALLELE' otherwise.  Lines that were not evaluated are
annotated as 'UNKNOWN'.  A sequence is considered a possible stutter
product if its total read count is less than or equal to the maximum
number of expected stutter reads.  The maximum number of stutter reads
is computed by assuming a fixed percentage of stutter product compared
to the originating sequence.

Stuttermark requires TSSV-style sequences as input (automatically
converting sequences to this format if necessary) and detects possible
stutter products by comparing sequences that have the same repeat blocks
but different numbers of repeats for one or more of their blocks.

The STUTTER annotation contains additional information.  For example:
'STUTTER:146.6x1(2-1):10.4x2(2-1x9-1)'.  This is a stutter product for
which at most 146.6 reads have come from the first sequence in the
output file ('146.6x1') and at most 10.4 reads have come from the second
sequence in the output file ('10.4x2').  This sequence differs from the
first sequence in the output file by a loss of one repeat of the second
repeat block ('2-1') and it differs from the second sequence by the loss
of one repeat in the second block and one repeat in the ninth block
('2-1x9-1').  (If this sequence would have more than 157 reads, it would
be annotated as 'ALLELE' instead.)
"""
import argparse, sys

from ..lib import pos_int_arg, print_db, PAT_TSSV_BLOCK, get_column_ids, \
                  add_input_output_args, get_input_output_files, \
                  ensure_sequence_format, add_sequence_format_args, \
                  SEQ_SPECIAL_VALUES

__version__ = "1.5.0"


# Default values for parameters are specified below.

# Link a difference in repeat count to a minimum occurrance percentage
# for stutter peaks.  Any repeat count differences other than listed
# here are taken to have zero probability.
# This value can be overridden by the -s command line option.
_DEF_STUTTERDEF = "-1:15,+1:4"

# Any peaks with less than this number of reads are not evaluated.
# This value can be overridden by the -m command line option.
_DEF_MIN_READS = 2

# Stutters are not expected to appear for blocks with less than
# _DEF_MIN_REPEATS repeats.  If two sequences have a different repeat
# count for any block with less than _DEF_MIN_REPEATS repeats, it is not
# possible that one of them is a stutter of the other.  Therefore, no
# comparison is made between any two such sequences.
# This value can be overridden by the -n command line option.
_DEF_MIN_REPEATS = 3

# Peaks are only annotated as a stutter peak of some other peak if the
# expected number of stutter occurances of this other peak is above
# this value.
# This value can be overridden by the -r command line option.
_DEF_MIN_REPORT = 0.1

# Default name of the column we are adding.
# This value can be overridden by the -c command line option.
_DEF_COLNAME = "annotation"


def load_data(infile, colname_annotation=_DEF_COLNAME, library=None):
    """
    Read data from the file handle infile and return a tuple
    (colum_names, data_rows).

    The input file is expected to have tab-delimited columns with
    column names on the first line and one sequence per data row.
    There are three required columns:
        marker    The name of the marker this sequence belongs to.
        sequence  The sequence in TSSV output format.
                  E.g., AGAT(7)TGAT(3).
        total     The total number of reads with this sequence.

    An additional column is appended.  All data rows
    are given a value of "UNKNOWN" for this column.

    :arg infile: Open readable handle to data file.
    :type infile: stream

    :arg colname_annotation: Name of the newly added column
    :type colname_annotation: str

    :arg library: Library for sequence format conversion
    :type library: dict

    :returns: A 2-tuple (column_names, data_rows).
    :rtype: tuple(list, list)
    """
    # Get column numbers.  We are adding a column as well.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_total, colid_sequence, = get_column_ids(column_names,
        "total", "sequence")
    colid_marker = get_column_ids(column_names, "marker", optional=True)
    column_names.append(colname_annotation)

    # Step through the file line by line to build the sequence list.
    allelelist = []
    for line in infile:
        columns = line.rstrip("\r\n").split("\t")

        # Skip empty/malformed lines (NOTE: allowing additional columns
        # beyond the expected 4 columns).
        if len(columns) != len(column_names)-1:
            if len(line.strip()):
                print("WARNING: skipped line: %s" % line)
            continue

        # String to integer conversion...
        columns[colid_total] = int(columns[colid_total])

        if columns[colid_sequence]:
            # Convert to TSSV-style sequences.
            marker = None if colid_marker is None else columns[colid_marker]
            columns[colid_sequence] = ensure_sequence_format(
                columns[colid_sequence], 'tssv', marker=marker,library=library)

            if columns[colid_sequence] not in SEQ_SPECIAL_VALUES:
                # Split the sequence column into a list of tuples:
                # [('ACTG', 4), ('CCTC', 12), ...]
                columns[colid_sequence] = [[x[0], int(x[1])] for x in
                    PAT_TSSV_BLOCK.findall(columns[colid_sequence])]

                # Repeat unit deduplication (data may contain stuff like
                # "AGAT(7)AGAT(5)" instead of "AGAT(12)").
                dedup = [columns[colid_sequence][0]]
                for i in range(1, len(columns[colid_sequence])):
                    if columns[colid_sequence][i][0] == dedup[-1][0]:
                        dedup[-1][1] += columns[colid_sequence][i][1]
                    else:
                        dedup.append(columns[colid_sequence][i])
                columns[colid_sequence] = dedup

        # Add the sequence to the list, including our new column.
        columns.append("UNKNOWN")
        allelelist.append(columns)
    return column_names, allelelist
#load_data


def annotate_alleles(infile, outfile, stutter, min_reads=_DEF_MIN_READS,
                     min_repeats=_DEF_MIN_REPEATS, min_report=_DEF_MIN_REPORT,
                     colname_annotation=_DEF_COLNAME, library=None,
                     debug=False):
    """
    Read data from the file handle infile and write annotated data to
    file handle outfile.

    The input file is expected to have tab-delimited columns with
    column names on the first line and one sequence per data row.
    There are three required columns:
        marker    The name of the marker this sequence belongs to.
        sequence  The sequence in TSSV output format.
                  E.g., AGAT(7)TGAT(3).
        total     The total number of reads with this sequence.

    An additional annotation column is appended.  All data rows
    are given a value of "ALLELE", "STUTTER:...", or "UNKNOWN".
    Sequences that are annotated as "STUTTER" are given a
    colon-separated description of their expected origins.  For each
    originating sequence a description of the form AxB(C) is given, with
    A the maximum expected number of stutter reads from this origin, B
    the sequence number of the originating sequence (where the first
    sequence in the output file is sequence 1), and C a 'x'-separated
    list of repeat blocks that had stuttered.

    Example stutter description:
        STUTTER:123.8x1(3-1):20.2x2(4+1)
    This stutter sequence has two originating sequences: a maximum of
    123.8 reads from sequence 1 plus a maximum of 20.2 reads from 2.
    Compared to sequence 1, this is a -1 stutter in the third repeat
    block.  Compared to sequence 2, this is a +1 stutter of the fourth
    repeat block.  (If this sequence had more than 144 total reads, it
    would have been annotated as "ALLELE".)

    :arg infile: Open readable handle to data file.
    :type infile: stream

    :arg outfile: Open writable handle for output data.
    :type outfile: stream

    :arg stutter: A dict mapping stutter size (keys) to a maximum
                  expected percentage of stutter of this size.
    :type stutter: dict{int: float}

    :arg min_reads: Any sequences with less than this number of reads are
                    not evaluated (annotation=UNKNOWN).
    :type min_reads: int

    :arg min_repeats: Stutters are not expected to appear for blocks
                      with less than this number of repeats.
    :type min_repeats: int

    :arg min_report: Stutters with an expected number of reads below
                     this value are not reported.
    :type min_report: float

    :arg colname_annotation: Name of the newly added column.
    :type colname_annotation: str

    :arg library: A parsed library file for sequence format conversion
    :type library: dict

    :arg debug: If True, print debug output to stdout.
    :type debug: bool
    """
    column_names, allelelist = load_data(infile, colname_annotation, library)
    colid_total, colid_sequence, = get_column_ids(column_names,
        "total", "sequence")
    try:
        colid_marker = get_column_ids(column_names, "marker")
    except:
        colid_marker = None

    # Sort (descending total reads).
    if colid_marker is not None:
        allelelist.sort(key=lambda allele:
            [allele[colid_marker], -allele[colid_total]])
    else:
        allelelist.sort(key=lambda allele: -allele[colid_total])

    for iCurrent in range(len(allelelist)):

        # See if this allele appeared a reasonable number of times.
        if allelelist[iCurrent][colid_total] < min_reads:
            continue

        # Skip special sequence values.
        if allelelist[iCurrent][colid_sequence] in SEQ_SPECIAL_VALUES:
            continue

        isStutterPeakOf = {}
        maximumOccurrenceExpected = 0

        # Find all ways to reach sequence iCurrent by mutating the other
        # sequences that appeared more often.
        for iOther in range(iCurrent):

            # Must be same marker.
            if colid_marker is not None and (allelelist[iCurrent][colid_marker]
                    != allelelist[iOther][colid_marker]):
                continue

            # Skip special sequence values.
            if allelelist[iOther][colid_sequence] in SEQ_SPECIAL_VALUES:
                continue

            print_db('%i vs %i' % (iCurrent+1, iOther+1), debug)

            # Must be same number of blocks.
            if (len(allelelist[iCurrent][colid_sequence]) !=
                    len(allelelist[iOther][colid_sequence])):
                print_db('Different number of blocks', debug)
                continue

            maximumOccurrenceExpectedForThisOther = 1.0
            blocksThatDiffer = {}

            # What is needed to get from sequence iOther to iCurrent?
            # We will look at each block in turn.
            for block in range(len(allelelist[iCurrent][colid_sequence])):

                # If mutations are needed to get from iOther to
                # iCurrent, iCurrent can't be a stutter of iOther.
                if (allelelist[iCurrent][colid_sequence][block][0] !=
                        allelelist[iOther][colid_sequence][block][0]):
                    print_db('Block %i has different sequence' % (block+1),
                             debug)
                    break

                # See how many times the current block appears more or
                # less often in iCurrent as compared to iOther.
                deltaRepeats = (allelelist[iCurrent][colid_sequence][block][1] -
                                allelelist[iOther][colid_sequence][block][1])
                if deltaRepeats == 0:
                    continue

                # iCurrent is not a stutter of iOther if any one of its
                # blocks differs in count while either count is below
                # min_repeats.  It is not a stutter of iOther because
                # those are expected not to occur for such low repeat
                # counts.  At the same time any further analysis
                # between the two is invalid because this block was not
                # the same.  Note that if this truly was a stutter of
                # iOther, maximumOccurrenceExpected is going to be a
                # bit too low now.
                if min(allelelist[iCurrent][colid_sequence][block][1],
                     allelelist[iOther][colid_sequence][block][1]) < min_repeats:
                    print_db('Block %i has low-count difference: %i and %i' %
                        (block+1, allelelist[iCurrent][colid_sequence][block][1],
                         allelelist[iOther][colid_sequence][block][1]), debug)
                    break

                # Depending on the repeat count difference, iCurrent
                # may appear at most a certain amount of times.
                # This is expressed as a percentage of the occurrence
                # of iOther.

                # If the repeat count difference is highly improbable,
                # this can't be a stutter.
                elif deltaRepeats not in stutter:
                    print_db('Improbable difference in number of repeats of '
                             'block %i: %i' % (block+1, deltaRepeats), debug)
                    break

                # Based on this block only, sequence iCurrent could be
                # a stutter of sequence iOther.  It may appear only a
                # limited number of times, however.
                else:
                    blocksThatDiffer[block] = deltaRepeats
                    maximumOccurrenceExpectedForThisOther *= \
                        stutter[deltaRepeats] / 100.0
                    print_db('Occurence percentage now %d' %
                             maximumOccurrenceExpectedForThisOther, debug)

            else:
                # If we end up here, iCurrent could very well be a
                # stutter peak of iOther, but it might be below the
                # reporting threshold.
                thisStutterOccurrenceExpected = \
                    (maximumOccurrenceExpectedForThisOther *
                     allelelist[iOther][colid_total])
                if thisStutterOccurrenceExpected < min_report:
                    print_db('Expected occurrence of stutter from '
                             '%i to %i is only %d' %
                             (iOther+1, iCurrent+1,
                              thisStutterOccurrenceExpected), debug)
                else:
                    isStutterPeakOf[iOther] = (blocksThatDiffer,
                        thisStutterOccurrenceExpected)
                    maximumOccurrenceExpected += thisStutterOccurrenceExpected
                    print_db('Expected occurence now %d' %
                             maximumOccurrenceExpected, debug)

        # Now we'll just need to check whether it does not appear
        # unrealistically often.
        if allelelist[iCurrent][colid_total] > maximumOccurrenceExpected:
            print_db('Occurs too often (%i > %d)' %
                     (allelelist[iCurrent][colid_total],
                      maximumOccurrenceExpected), debug)
            allelelist[iCurrent][-1] = "ALLELE"
        else:
            # Stutter peak detected.
            allelelist[iCurrent][-1] = "STUTTER:%s" % ":".join(
                "%.1fx%i(%s)" % (
                    isStutterPeakOf[iOther][1],
                    iOther+1,
                    "x".join(
                        "%i%+i" % (
                            block+1,
                            isStutterPeakOf[iOther][0][block])
                        for block in isStutterPeakOf[iOther][0]))
                for iOther in isStutterPeakOf)
            print_db("%2i is a stutter peak of %s; occur=%i, maxOccur=%s" %
                     (iCurrent+1, isStutterPeakOf,
                      allelelist[iCurrent][colid_total],
                      maximumOccurrenceExpected), debug)

    # Reconstruct the sequence and write out the findings.
    for i in range(len(allelelist)):
        if allelelist[i][colid_sequence] not in SEQ_SPECIAL_VALUES:
            allelelist[i][colid_sequence] = "".join(
                "%s(%i)" % tuple(block)
                    for block in allelelist[i][colid_sequence])
    outfile.write("\t".join(column_names))
    outfile.write("\n")
    outfile.write("\n".join(
        "\t".join(map(str, allele)) for allele in allelelist))
    outfile.write("\n")
#annotate_alleles


def stutter_def_arg(value):
    """
    Convert a stutter definition string to dict, raise ArgumentTypeError
    on failure.

    The stutter definition string is a comma-separated list of int:float
    pairs, e.g., "-1:15,+1:3.5".  Whitespace around interpunction is
    permitted, e.g., "-1: 15, +1: 3.5".

    :arg min_reads: A valid stutter definition string.
    :type min_reads: str

    :returns: A dict mapping stutter size (keys) to a maximum expected
              percentage of stutter of this size.
    :rtype: dict{int: float}
    """
    if ":" not in value:
        return {}
    try:
        return {int(y[0]): float(y[1])
                for y in (x.split(":") for x in value.split(","))}
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            "invalid stutter definition value: '%s' (%s)" % (value, str(e)))
#stutter_def_arg


def add_arguments(parser):
    add_input_output_args(parser, True, True, False)
    parser.add_argument('-s', '--stutter', metavar="DEF", type=stutter_def_arg,
        default=_DEF_STUTTERDEF,
        help="Define maximum expected stutter percentages.  The default value "
             "of '%(default)s' sets -1 stutter (loss of one repeat) to 15%%, "
             "+1 stutter (gain of one repeat) to 4%%.  Any unspecified "
             "stutter amount is assumed not to occur directly but e.g., a -2 "
             "stutter may still be recognised as two -1 stutters stacked "
             "together.  NOTE: It may be necessary to specify this option as "
             "'-s=%(default)s' (note the equals sign instead of a space).")
    parser.add_argument('-c', '--column-name', metavar="COLNAME",
        default=_DEF_COLNAME,
        help="name of the newly added column (default: '%(default)s')")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-m', '--min-reads', metavar="N",
        type=pos_int_arg, default=_DEF_MIN_READS,
        help="set minimum number of reads to evaluate (default: %(default)s)")
    filtergroup.add_argument('-n', '--min-repeats', metavar="N",
        type=pos_int_arg, default=_DEF_MIN_REPEATS,
        help="set minimum number of repeats of a block that can possibly "
             "stutter (default: %(default)s)")
    filtergroup.add_argument('-r', '--min-report', metavar="N", type=float,
        default=_DEF_MIN_REPORT,
        help="a sequence is only annotated as a stutter of some other "
             "sequence if the expected number of stutter occurances of this "
             "other sequence is above this value (default: %(default)s)")
    add_sequence_format_args(parser, "tssv", True)  # Force tssv seqs.
#add_arguments


def run(args):
    gen = get_input_output_files(args, True, True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        # This tool now only works properly with one infile per sample!
        if len(infiles) > 1:
            raise ValueError(
                "multiple input files for sample '%s' specified" % tag)
        infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "r")
        annotate_alleles(infile, outfile, args.stutter, args.min_reads,
                         args.min_repeats, args.min_report, args.column_name,
                         args.library, args.debug)
        if infile != sys.stdin:
            infile.close()
#run
