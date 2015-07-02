#!/usr/bin/env python
"""
Mark potential stutter products by assuming a fixed maximum percentage
of stutter product vs the parent allele.
"""
import argparse
import sys
import re

from ..lib import pos_int_arg, print_db, PAT_TSSV_BLOCK, get_column_ids

__version__ = "1.3"


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
# _DEF_MIN_REPEATS repeats.  If two alleles have a different repeat count
# for any block with less than _DEF_MIN_REPEATS repeats, it is not possible
# that one of them is a stutter of the other.  Therefore, no comparison
# is made between any two such sequences.
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


def load_data(infile, colname_annotation=_DEF_COLNAME):
    """
    Read data from the file handle infile and return a tuple
    (colum_names, data_rows).

    The input file is expected to have tab-delimited columns with
    column names on the first line and one allele per data row.
    There are three required columns:
        name    The name of the marker this allele belongs to.
        allele  The allele sequence in TSSV output format.
                E.g., AGAT(7)TGAT(3).
        total   The total number of reads with this sequence.

    An additional column "annotation" is appended.  All data rows
    are given a value of "UNKNOWN" for this column.

    :arg infile: Open readable handle to data file.
    :type infile: stream

    :returns: A 2-tuple (column_names, data_rows).
    :rtype: tuple(list, list)
    """
    # Get column numbers.  We are adding a column as well.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_total, colid_allele, colid_name = get_column_ids(column_names,
        "total", "allele", "name")
    column_names.append(colname_annotation)

    # Step through the file line by line to build the allele list.
    allelelist = []
    for line in infile:
        columns = line.rstrip("\r\n").split("\t")

        # Skip empty/malformed lines (NOTE: allowing additional columns
        # beyond the expected 4 columns).
        if len(columns) != len(column_names)-1:
            if len(line.strip()):
                print("WARNING: skipped line: %s" % line)
            continue

        # Split the allele column into a list of tuples:
        # [('ACTG','4'),('CCTC','12'),...]
        columns[colid_allele] = PAT_TSSV_BLOCK.findall(columns[colid_allele])
        if columns[colid_allele] == None:
            print("WARNING: skipped line: %s" % line)
            continue

        # String to integer conversion...
        columns[colid_allele] = map(
            lambda x: [x[0], int(x[1])], columns[colid_allele])
        columns[colid_total] = int(columns[colid_total])

        # Repeat unit deduplication (data may contain stuff like
        # "AGAT(7)AGAT(5)" instead of "AGAT(12)").
        if columns[colid_allele]:
            dedup = [columns[colid_allele][0]]
            for i in range(1, len(columns[colid_allele])):
                if columns[colid_allele][i][0] == dedup[-1][0]:
                    dedup[-1][1] += columns[colid_allele][i][1]
                else:
                    dedup.append(columns[colid_allele][i])
            columns[colid_allele] = dedup

        # Add the allele to the list, including our new column.
        columns.append("UNKNOWN")
        allelelist.append(columns)
    return column_names, allelelist
#load_data


def annotate_alleles(infile, outfile, stutter, min_reads=_DEF_MIN_READS,
                     min_repeats=_DEF_MIN_REPEATS, min_report=_DEF_MIN_REPORT,
                     colname_annotation=_DEF_COLNAME, debug=False):
    """
    Read data from the file handle infile and write annotated data to
    file handle outfile.

    The input file is expected to have tab-delimited columns with
    column names on the first line and one allele per data row.
    There are three required columns:
        name    The name of the marker this allele belongs to.
        allele  The allele sequence in TSSV output format.
                E.g., AGAT(7)TGAT(3).
        total   The total number of reads with this sequence.

    An additional annotation column is appended.  All data rows
    are given a value of "ALLELE", "STUTTER:...", or "UNKNOWN".
    Alleles that are annotated as "STUTTER" are given a colon-separated
    description of their expected origins.  For each originating allele
    a description of the form AxB(C) is given, with A the maximum
    expected number of stutter reads from this origin, B the allele
    number of the originating allele (where the first allele in the
    output file is allele 1), and C a 'x'-separated list of repeat
    blocks that had stuttered.

    Example stutter description:
        STUTTER:123.8x1(3-1):20.2x2(4+1)
    This stutter allele has two originating alleles: a maximum of 123.8
    reads from allele 1 plus a maximum of 20.2 reads from 2.  Compared
    to allele 1, this is a -1 stutter in the third repeat block.
    Compared to allele 2, this is a +1 stutter of the fourth repeat
    block.  (If this allele had more than 144 total reads, it would have
    been annotated as "ALLELE".)

    :arg infile: Open readable handle to data file.
    :type infile: stream

    :arg outfile: Open writable handle for output data.
    :type outfile: stream

    :arg stutter: A dict mapping stutter size (keys) to a maximum
                  expected percentage of stutter of this size.
    :type stutter: dict{int: float}

    :arg min_reads: Any alleles with less than this number of reads are
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

    :arg debug: If True, print debug output to stdout.
    :type debug: bool
    """
    column_names, allelelist = load_data(infile, colname_annotation)
    colid_total, colid_allele, colid_name = get_column_ids(column_names,
        "total", "allele", "name")

    # Sort (descending total reads).
    allelelist.sort(key=lambda allele:
        [allele[colid_name], -allele[colid_total]])

    for iCurrent in range(len(allelelist)):

        # See if this allele appeared a reasonable number of times.
        if allelelist[iCurrent][colid_total] < min_reads:
            continue

        isStutterPeakOf = {}
        maximumOccurrenceExpected = 0

        # Find all ways to reach allele iCurrent by mutating the other
        # alleles that appeared more often.
        for iOther in range(iCurrent):

            # Must be same marker.
            if (allelelist[iCurrent][colid_name] !=
                    allelelist[iOther][colid_name]):
                continue

            print_db('%i vs %i' % (iCurrent+1, iOther+1), debug)

            # Must be same number of blocks.
            if (len(allelelist[iCurrent][colid_allele]) !=
                    len(allelelist[iOther][colid_allele])):
                print_db('Different number of blocks', debug)
                continue

            maximumOccurrenceExpectedForThisOther = 1.0
            blocksThatDiffer = {}

            # What is needed to get from sequence iOther to iCurrent?
            # We will look at each block in turn.
            for block in range(len(allelelist[iCurrent][colid_allele])):

                # If mutations are needed to get from iOther to
                # iCurrent, iCurrent can't be a stutter of iOther.
                if (allelelist[iCurrent][colid_allele][block][0] !=
                        allelelist[iOther][colid_allele][block][0]):
                    print_db('Block %i has different sequence' % (block+1),
                             debug)
                    break

                # See how many times the current block appears more or
                # less often in iCurrent as compared to iOther.
                deltaRepeats = (allelelist[iCurrent][colid_allele][block][1] -
                                allelelist[iOther][colid_allele][block][1])
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
                if min(allelelist[iCurrent][colid_allele][block][1],
                     allelelist[iOther][colid_allele][block][1]) < min_repeats:
                    print_db('Block %i has low-count difference: %i and %i' %
                        (block+1, allelelist[iCurrent][colid_allele][block][1],
                         allelelist[iOther][colid_allele][block][1]), debug)
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
            allelelist[iCurrent][-1] = "STUTTER:%s" % ":".join(map(
                lambda iOther: "%.1fx%i(%s)" % (
                    isStutterPeakOf[iOther][1],
                    iOther+1,
                    "x".join(map(
                        lambda block: "%i%+i" % (
                            block+1,
                            isStutterPeakOf[iOther][0][block]),
                        isStutterPeakOf[iOther][0]))),
                isStutterPeakOf))
            print_db("%2i is a stutter peak of %s; occur=%i, maxOccur=%s" %
                     (iCurrent+1, isStutterPeakOf,
                      allelelist[iCurrent][colid_total],
                      maximumOccurrenceExpected), debug)

    # Reconstruct the allele sequence and write out the findings.
    for i in range(len(allelelist)):
        allelelist[i][colid_allele] = "".join(map(
            lambda block: "%s(%i)" % tuple(block),
            allelelist[i][colid_allele]))
    outfile.write("\t".join(column_names))
    outfile.write("\n")
    outfile.write("\n".join(map(
        lambda allele: "\t".join(map(str, allele)), allelelist)))
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
        return {int(x[0]): float(x[1])
                for x in map(lambda x: x.split(":"), value.split(","))}
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            "invalid stutter definition value: '%s' (%s)" % (value, str(e)))
#stutter_def_arg


def add_arguments(parser):
    parser.add_argument('infile', nargs='?', metavar="IN", default=sys.stdin,
        type=argparse.FileType('r'),
        help="the CSV data file to process (default: read from stdin)")
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-s', '--stutter', metavar="DEF", type=stutter_def_arg,
        default=_DEF_STUTTERDEF,
        help="Define maximum expected stutter percentages.  The default value "
             "of '%(default)s' sets -1 stutter (loss of one repeat) to 15%%, "
             "+1 stutter (gain of one repeat) to 4%%.  Any unspecified "
             "stutter amount is assumed not to occur directly but e.g., a -2 "
             "stutter may still be recognised as two -1 stutters stacked"
             "together.  NOTE: It may be necessary to specify this option as "
             "'-s=-1:15,+1:2' (note the equals sign instead of a space).")
    parser.add_argument('-m', '--min-reads', metavar="N", type=pos_int_arg,
        default=_DEF_MIN_READS,
        help="set minimum number of reads to evaluate (default: %(default)s)")
    parser.add_argument('-n', '--min-repeats', metavar="N", type=pos_int_arg,
        default=_DEF_MIN_REPEATS,
        help="set minimum number of repeats of a block that can possibly "
             "stutter (default: %(default)s)")
    parser.add_argument('-r', '--min-report', metavar="N", type=float,
        default=_DEF_MIN_REPORT,
        help="alleles are only annotated as a stutter of some other allele if "
             "the expected number of stutter occurances of this other allele "
             "is above this value (default: %(default)s)")
    parser.add_argument('-c', '--column-name', metavar="COLNAME",
        default=_DEF_COLNAME,
        help="name of the newly added column (default: %(default)s)")
    parser.add_argument('-d', '--debug', action="store_true",
        help="if specified, debug output is printed to stdout")
#add_arguments


def run(args):
    annotate_alleles(args.infile, args.outfile, args.stutter,
                     args.min_reads, args.min_repeats, args.min_report,
                     args.column_name, args.debug)
#run


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    try:
        add_arguments(parser)
        run(parser.parse_args())
    except OSError as error:
        parser.error(error)
#main


if __name__ == "__main__":
    main()
