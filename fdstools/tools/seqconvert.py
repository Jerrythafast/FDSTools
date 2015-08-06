#!/usr/bin/env python
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
the help text of the libconvert tool for details.

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
import argparse
import sys

from ..lib import get_column_ids, ensure_sequence_format, parse_library,\
                  reverse_complement

__version__ = "0.1dev"


# Default values for parameters are specified below.

# Default name of the column that contains the marker name.
# This value can be overridden by the -m command line option.
_DEF_COLNAME_MARKER = "name"

# Default name of the column that contains the allele.
# This value can be overridden by the -a command line option.
_DEF_COLNAME_ALLELE = "allele"

# Default name of the column to write the output to.
# This value can be overridden by the -o command line option.
_DEF_COLNAME_ALLELE_OUT = "allele"


def convert_sequences(infile, outfile, to_format, libfile=None,
                      fixed_marker=None, colname_marker=_DEF_COLNAME_MARKER,
                      colname_allele=_DEF_COLNAME_ALLELE,
                      colname_allele_out=_DEF_COLNAME_ALLELE_OUT,
                      libfile2=None, revcomp_markers=[]):
    libfile = libfile if libfile is not None else libfile2
    library = parse_library(libfile) if libfile is not None else None
    library2 = parse_library(libfile2) if libfile2 is not None else library
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_allele = get_column_ids(column_names, colname_allele)
    if library is None:
        fixed_marker = ""  # Don't need marker names without library.
    if fixed_marker is None:
        colid_marker = get_column_ids(column_names, colname_marker)
    try:
        colid_allele_out = get_column_ids(column_names, colname_allele_out)
    except:
        column_names.append(colname_allele_out)
        colid_allele_out = -1

    outfile.write("\t".join(column_names) + "\n")
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        if colid_allele_out == -1:
            line.append("")
        marker = line[colid_marker] if fixed_marker is None else fixed_marker

        seq = line[colid_allele]
        if library2 != library:
            seq = ensure_sequence_format(
                seq, "raw", marker=marker, library=library)
            if marker in revcomp_markers:
                seq = reverse_complement(seq)
            # TODO: The current implementation assumes identical
            # flanking sequences.  Introduce means to shift flanking
            # sequence in/out of prefix and/or suffix.

        line[colid_allele_out] = ensure_sequence_format(
            seq, to_format, marker=marker, library=library2)
        outfile.write("\t".join(line) + "\n")
#convert_sequences


def add_arguments(parser):
    parser.add_argument('format', metavar="FORMAT",
        choices=("raw", "tssv", "allelename"),
        help="the format to convert to: one of %(choices)s")
    parser.add_argument('infile', nargs='?', metavar="IN", default=sys.stdin,
        type=argparse.FileType('r'),
        help="the tab-separated data file to process (default: read from "
             "stdin)")
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-m', '--marker-column', metavar="COLNAME",
        default=_DEF_COLNAME_MARKER,
        help="name of the column that contains the marker name "
             "(default: '%(default)s')")
    parser.add_argument('-a', '--allele-column', metavar="COLNAME",
        default=_DEF_COLNAME_ALLELE,
        help="name of the column that contains the allele "
             "(default: '%(default)s')")
    parser.add_argument('-o', '--output-column', metavar="COLNAME",
        default=_DEF_COLNAME_ALLELE_OUT,
        help="name of the column to write the output to "
             "(default: '%(default)s')")
    parser.add_argument('-M', '--marker', metavar="MARKER",
        help="assume the specified marker for all sequences")
    parser.add_argument('-l', '--library', metavar="LIBRARY",
        type=argparse.FileType('r'),
        help="library file for sequence format conversion")
    parser.add_argument('-L', '--library2', metavar="LIBRARY",
        type=argparse.FileType('r'),
        help="second library file to use for output; if specified, allele "
             "names can be conveniently updated to fit this new library file")
    parser.add_argument('-R', '--reverse-complement', metavar="MARKER",
        nargs="+", default=[],
        help="to be used togethwer with -L/--library2; specify the markers "
             "for which the sequences are reverse-complemented in the new "
             "library")
#add_arguments


def run(args):
    if args.infile.isatty() and args.outfile.isatty():
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    convert_sequences(args.infile, args.outfile, args.format, args.library,
                      args.marker, args.marker_column, args.allele_column,
                      args.output_column, args.library2,
                      args.reverse_complement)
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
