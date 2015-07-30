#!/usr/bin/env python
"""
Match background noise profiles to samples.
"""
import argparse
import sys
#import numpy as np  # Only imported when actually running this tool.

from ..lib import parse_library, load_profiles, ensure_sequence_format, nnls, \
                  get_column_ids

__version__ = "0.1dev"


def get_sample_data(infile, convert_to_raw=False, library=None):
    """
    Read data from the given file handle, corresponding to a single
    sample, and fill a dict with all sequences in the sample.
    """
    column_names = infile.readline().rstrip("\r\n").split("\t")
    column_names.append("forward_noise")
    column_names.append("reverse_noise")
    column_names.append("total_noise")
    column_names.append("forward_add")
    column_names.append("reverse_add")
    column_names.append("total_add")
    colid_name, colid_allele, colid_forward, colid_reverse = get_column_ids(
        column_names, "name", "allele", "forward", "reverse")
    data = {}
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        marker = cols[colid_name]
        if convert_to_raw:
            cols[colid_allele] = ensure_sequence_format(
                cols[colid_allele], "raw", library=library, marker=marker)
        cols[colid_forward] = int(cols[colid_forward])
        cols[colid_reverse] = int(cols[colid_reverse])
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        if marker not in data:
            data[marker] = []
        data[marker].append(cols)
    return column_names, data
#get_sample_data


def match_profile(column_names, data, profile, convert_to_raw, library,
                  marker):
    import numpy as np
    (colid_name, colid_allele, colid_forward, colid_reverse, colid_total,
     colid_forward_noise, colid_reverse_noise, colid_total_noise,
     colid_forward_add, colid_reverse_add, colid_total_add) = get_column_ids(
        column_names, "name", "allele", "forward", "reverse", "total",
        "forward_noise", "reverse_noise", "total_noise", "forward_add",
        "reverse_add", "total_add")

    # Enter profiles into P.
    P1 = np.matrix(profile["forward"])
    P2 = np.matrix(profile["reverse"])

    # Enter sample into C.
    seqs = []
    C1 = np.matrix(np.zeros([1, profile["m"]]))
    C2 = np.matrix(np.zeros([1, profile["m"]]))
    for line in data:
        if convert_to_raw:
            allele = ensure_sequence_format(line[colid_allele], "raw",
                                            library=library, marker=marker)
        else:
            allele = line[colid_allele]
        seqs.append(allele)
        try:
            i = profile["seqs"].index(allele)
        except ValueError:
            # Note: Not adding any new sequences to the profile, since
            # they will just be zeroes and have no effect on the result.
            continue
        C1[0, i] = line[colid_forward]
        C2[0, i] = line[colid_reverse]

    # Compute corrected read counts.
    A = nnls(np.hstack([P1, P2]).T, np.hstack([C1, C2]).T)
    np.fill_diagonal(P1, 0)
    np.fill_diagonal(P2, 0)
    forward_noise = A.T * P1
    reverse_noise = A.T * P2
    forward_add = np.multiply(A, P1.sum(1)).T
    reverse_add = np.multiply(A, P2.sum(1)).T

    j = 0
    for line in data:
        j += 1
        try:
            i = profile["seqs"].index(seqs[j-1])
        except ValueError:
            continue
        line[colid_forward_noise] = forward_noise[0, i]
        line[colid_reverse_noise] = reverse_noise[0, i]
        line[colid_total_noise] = forward_noise[0, i] + reverse_noise[0, i]
        if i < profile["n"]:
            line[colid_forward_add] = forward_add[0, i]
            line[colid_reverse_add] = reverse_add[0, i]
            line[colid_total_add] = forward_add[0, i] + reverse_add[0, i]

    # Add sequences that are in the profile but not in the sample.
    for i in range(profile["m"]):
        if profile["seqs"][i] in seqs:
            continue
        amount = forward_noise[0, i] + reverse_noise[0, i]
        if i < profile["n"]:
            amount += forward_add[0, i] + reverse_add[0, i]
        if amount > 0:
            line = [""] * len(column_names)
            line[colid_name] = marker
            line[colid_allele] = profile["seqs"][i]
            line[colid_forward] = 0
            line[colid_reverse] = 0
            line[colid_total] = 0
            line[colid_forward_noise] = forward_noise[0, i]
            line[colid_reverse_noise] = reverse_noise[0, i]
            line[colid_total_noise] = forward_noise[0, i] + reverse_noise[0, i]
            if i < profile["n"]:
                line[colid_forward_add] = forward_add[0, i]
                line[colid_reverse_add] = reverse_add[0, i]
                line[colid_total_add] = forward_add[0, i] + reverse_add[0, i]
            else:
                line[colid_forward_add] = 0
                line[colid_reverse_add] = 0
                line[colid_total_add] = 0
            data.append(line)
#match_profile


def match_profiles(profilefile, infile, outfile, libfile, seqformat, marker):
    library = parse_library(libfile) if libfile else None
    profiles = load_profiles(profilefile, library)
    if marker:
        profiles = {marker: profiles[marker]} if marker in profiles else {}

    column_names, data = get_sample_data(
        infile, convert_to_raw=seqformat=="raw", library=library)
    colid_allele = get_column_ids(column_names, "allele")

    outfile.write("\t".join(column_names) + "\n")
    for marker in data:
        if marker in profiles:
            match_profile(column_names, data[marker], profiles[marker],
                          seqformat!="raw", library, marker)
        for line in data[marker]:
            if seqformat is not None and seqformat != "raw":
                line[colid_allele] = ensure_sequence_format(line[colid_allele],
                    seqformat, library=library, marker=marker)
            outfile.write("\t".join(map(str, line)) + "\n")
#match_profiles


def add_arguments(parser):
    parser.add_argument('profiles', metavar="PROFILES",
        type=argparse.FileType('r'),
        help="file containing background noise profiles to match")
    parser.add_argument('infile', nargs='?', metavar="IN", default=sys.stdin,
        type=argparse.FileType('r'),
        help="the tab-separated data file to process (default: read from "
             "stdin)")
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-F', '--sequence-format', metavar="FORMAT",
        choices=["raw", "tssv", "allelename"],
        help="convert sequences to the specified format: one of %(choices)s "
             "(default: no conversion)")
    parser.add_argument('-l', '--library', metavar="LIBRARY",
        type=argparse.FileType('r'),
        help="library file for sequence format conversion")
    parser.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
#add_arguments


def run(args):
    if args.infile.isatty() and args.outfile.isatty():
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    match_profiles(args.profiles, args.infile, args.outfile, args.library,
                   args.sequence_format, args.marker)
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
