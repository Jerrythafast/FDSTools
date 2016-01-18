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
'-' is implicitly used as the second input file.  Note that as a result,
in case of conflicting values, the value in the specified input file
will take precedence over the value in the data that was piped in.

Example: fdstools bgpredict ... | fdstools bgmerge old.txt > out.txt
"""
import argparse
import sys

from ..lib import load_profiles, ensure_sequence_format,\
                  add_sequence_format_args

__version__ = "0.1dev"


def merge_profiles(infiles, outfile, seqformat, library):
    amounts = {}
    for infile in infiles:
        if infile == "-":
            profiles = load_profiles(sys.stdin, library)
        else:
            with open(infile, "r") as handle:
                profiles = load_profiles(handle, library)
        for marker in profiles:
            if marker not in amounts:
                amounts[marker] = {}
            for i in range(profiles[marker]["n"]):
                for j in range(profiles[marker]["m"]):
                    key = (profiles[marker]["seqs"][i],
                           profiles[marker]["seqs"][j])
                    if key not in amounts[marker]:
                        this_amounts = (profiles[marker]["forward"][i][j],
                                        profiles[marker]["reverse"][i][j],
                                        profiles[marker]["tool"][i][j])
                        if sum(this_amounts[:2]):
                            amounts[marker][key] = this_amounts

    outfile.write("\t".join(
        ["marker", "allele", "sequence", "fmean", "rmean", "tool"]) + "\n")
    for marker in amounts:
        for allele, sequence in amounts[marker]:
            outfile.write("\t".join([marker] +
                [ensure_sequence_format(seq, seqformat, library=library,
                    marker=marker) for seq in (allele, sequence)] +
                map(str, amounts[marker][allele, sequence])) + "\n")
#merge_profiles


def add_arguments(parser):
    parser.add_argument('infiles', nargs='+', metavar="FILE",
        help="files containing the background noise profiles to combine; "
             "if a single file is given, it is merged with input from stdin; "
             "use '-' to use stdin as an explicit input source")
    outgroup = parser.add_argument_group("output file options")
    outgroup.add_argument('-o', '--output', dest="outfile", metavar="FILE",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="file to write output to (default: write to stdout)")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
#add_arguments


def run(args):
    if len(args.infiles) < 2:
        if sys.stdin.isatty() or "-" in args.infiles:
            raise ValueError("please specify at least two input files")
        args.infiles.append("-")

    merge_profiles(args.infiles, args.outfile, args.sequence_format,
                   args.library)
#run
