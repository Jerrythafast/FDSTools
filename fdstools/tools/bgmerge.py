#!/usr/bin/env python
"""
Merge multiple files containing background noise profiles.

Background noise profiles are merged in the order in which they are
specified.  If multple files specify a different value for the same
allele and sequence, the value of the first file is used.

It is convenient to pipe the output of bgpredict and/or bgestimate
into bgmerge to merge that with an existing file containing background
profiles.  Specify '-' as one of the input files to read from stdin
(i.e., read input from a pipe).  If only one input file is specified,
'-' is implicitly used as the second input file.

Example: fdstools bgpredict ... | fdstools bgmerge old.txt > out.txt
"""
import argparse
import sys

from ..lib import load_profiles, ensure_sequence_format, parse_library,\
                  add_sequence_format_args

__version__ = "0.1dev"


def merge_profiles(infiles, outfile, crosstab, seqformat, library):
    # Parse library file.
    library = parse_library(library) if library is not None else None

    amounts = {}
    for infile in infiles:
        profiles = load_profiles(infile, library)
        for marker in profiles:
            if marker not in amounts:
                amounts[marker] = {}
            for i in range(profiles[marker]["n"]):
                for j in range(profiles[marker]["m"]):
                    key = (profiles[marker]["seqs"][i],
                           profiles[marker]["seqs"][j])
                    if key not in amounts[marker]:
                        this_amounts = (profiles[marker]["forward"][i][j],
                                        profiles[marker]["reverse"][i][j])
                        if sum(this_amounts):
                            amounts[marker][key] = this_amounts

    if not crosstab:
        # Tab-separated columns format.
        outfile.write("\t".join(
            ["marker", "allele", "sequence", "fmean", "rmean"]) + "\n")
        for marker in amounts:
            for allele, sequence in amounts[marker]:
                outfile.write("\t".join([marker] +
                    [ensure_sequence_format(seq, seqformat, library=library,
                        marker=marker) for seq in (allele, sequence)] +
                    map(str, amounts[marker][allele, sequence])) + "\n")
        return

    # Cross-tabular format (profiles in rows).
    for marker in amounts:
        alleles = set(allele for allele, sequence in amounts[marker])
        seqs = list(alleles) + list(set(
            sequence for allele, sequence in amounts[marker]
            if sequence not in alleles))
        outfile.write("\t".join([marker, "0"] +
            [ensure_sequence_format(seq, seqformat, library=library,
                marker=marker) for seq in seqs]) + "\n")
        forward = [[0 if (allele, sequence) not in amounts[marker] else
                    amounts[marker][allele, sequence][0] for sequence in seqs]
                   for allele in seqs[:len(alleles)]]
        reverse = [[0 if (allele, sequence) not in amounts[marker] else
                    amounts[marker][allele, sequence][1] for sequence in seqs]
                   for allele in seqs[:len(alleles)]]
        for i in range(len(alleles)):
            outfile.write("\t".join(
                [marker, str(i+1)] + map(str, forward[i])) + "\n")
            outfile.write("\t".join(
                [marker, str(-i-1)] + map(str, reverse[i])) + "\n")
#merge_profiles


def add_arguments(parser):
    parser.add_argument('infiles', nargs='+', metavar="FILE",
        type=argparse.FileType('r'),
        help="files containing the background noise profiles to combine; "
             "if a single file is given, it is merged with input from stdin; "
             "use '-' to use stdin as an explicit input source")
    outgroup = parser.add_argument_group("output file options")
    outgroup.add_argument('-o', '--output', dest="outfile", metavar="FILE",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="file to write output to (default: write to stdout)")
    parser.add_argument('-C', '--cross-tabular', action="store_true",
        help="if specified, a space-efficient cross-tabular output format is "
             "used instead of the default tab-separated columns format")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
#add_arguments


def run(args):
    if len(args.infiles) < 2:
        if sys.stdin.isatty() or sys.stdin in args.infiles:
            raise ValueError("please specify at least two input files")
        args.infiles.append(sys.stdin)

    merge_profiles(args.infiles, args.outfile, args.cross_tabular,
                   args.sequence_format, args.library)
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