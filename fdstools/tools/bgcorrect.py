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
Match background noise profiles (obtained from e.g., bgestimate) to
samples.

Ten new columns are added to the output giving, for each sequence, the
number of reads attributable to noise from other sequences (_noise
columns) and the number of noise reads caused by the prescense of this
sequence (_add columns), as well as the resulting number of reads after
correction (_corrected columns: original minus _noise plus _add).

The correction_flags column contains one of the following values:
'not_corrected', no background noise profile was available for this
marker; 'not_in_ref_db', the sequence was not present in the noise
profiles given; 'corrected_as_background_only', the sequence was present
in the noise profiles given, but only as noise and not as genuine
allele; 'corrected_bgpredict', the sequence was present in the noise
profiles as a genuine allele, but its noise profile consists entirely of
predictions as opposed to direct observations;
'corrected_bgestimate'/'corrected_bghomstats', the sequence was present
in the noise profiles as a genuine allele and at least part of its noise
profile was based on direct observations.
"""
import argparse, sys
#import numpy as np  # Only imported when actually running this tool.

from ..lib import load_profiles, ensure_sequence_format, get_column_ids, \
                  nnls, add_sequence_format_args, SEQ_SPECIAL_VALUES, \
                  add_input_output_args, get_input_output_files

__version__ = "1.0.1"


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
    column_names.append("forward_corrected")
    column_names.append("reverse_corrected")
    column_names.append("total_corrected")
    column_names.append("correction_flags")
    column_names.append("weight")
    colid_marker, colid_sequence, colid_forward, colid_reverse =get_column_ids(
        column_names, "marker", "sequence", "forward", "reverse")
    data = {}
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        marker = cols[colid_marker]
        if convert_to_raw:
            cols[colid_sequence] = ensure_sequence_format(
                cols[colid_sequence], "raw", library=library, marker=marker)
        cols[colid_forward] = int(cols[colid_forward])
        cols[colid_reverse] = int(cols[colid_reverse])
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(0)
        cols.append(cols[colid_forward])
        cols.append(cols[colid_reverse])
        cols.append(cols[colid_forward] + cols[colid_reverse])
        cols.append("not_corrected")
        cols.append((cols[colid_forward] + cols[colid_reverse]) / 100.)
        if marker not in data:
            data[marker] = []
        data[marker].append(cols)
    return column_names, data
#get_sample_data


def match_profile(column_names, data, profile, convert_to_raw, library,
                  marker):
    (colid_marker, colid_sequence, colid_forward, colid_reverse, colid_total,
     colid_forward_noise, colid_reverse_noise, colid_total_noise,
     colid_forward_add, colid_reverse_add, colid_total_add,
     colid_forward_corrected, colid_reverse_corrected,
     colid_total_corrected, colid_correction_flags, colid_weight) = \
        get_column_ids(column_names, "marker", "sequence", "forward",
        "reverse", "total", "forward_noise", "reverse_noise", "total_noise",
        "forward_add", "reverse_add", "total_add", "forward_corrected",
        "reverse_corrected", "total_corrected", "correction_flags", "weight")

    # Enter profiles into P.
    P1 = np.matrix(profile["forward"])
    P2 = np.matrix(profile["reverse"])

    # Enter sample into C.
    seqs = []
    C1 = np.matrix(np.zeros([1, profile["m"]]))
    C2 = np.matrix(np.zeros([1, profile["m"]]))
    for line in data:
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        if convert_to_raw:
            sequence = ensure_sequence_format(line[colid_sequence], "raw",
                                            library=library, marker=marker)
        else:
            sequence = line[colid_sequence]
        seqs.append(sequence)
        try:
            i = profile["seqs"].index(sequence)
        except ValueError:
            # Note: Not adding any new sequences to the profile, since
            # they will just be zeroes and have no effect on the result.
            continue
        C1[0, i] = line[colid_forward]
        C2[0, i] = line[colid_reverse]

    # Stop if this sample has no explicit data for this marker.
    if not seqs:
        return

    # Compute corrected read counts.
    A = nnls(np.hstack([P1, P2]).T, np.hstack([C1, C2]).T).T
    np.fill_diagonal(P1, 0)
    np.fill_diagonal(P2, 0)
    forward_noise = A * P1
    reverse_noise = A * P2
    forward_add = np.multiply(A, P1.sum(1).T)
    reverse_add = np.multiply(A, P2.sum(1).T)

    # Round values to 3 decimal positions.
    A.round(3, A)
    forward_noise.round(3, forward_noise)
    reverse_noise.round(3, reverse_noise)
    forward_add.round(3, forward_add)
    reverse_add.round(3, reverse_add)

    j = 0
    for line in data:
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        j += 1
        try:
            i = profile["seqs"].index(seqs[j-1])
        except ValueError:
            line[colid_correction_flags] = "not_in_ref_db"
            continue
        line[colid_forward_noise] = forward_noise[0, i]
        line[colid_reverse_noise] = reverse_noise[0, i]
        line[colid_total_noise] = forward_noise[0, i] + reverse_noise[0, i]
        line[colid_forward_corrected] -= line[colid_forward_noise]
        line[colid_reverse_corrected] -= line[colid_reverse_noise]
        line[colid_total_corrected] -= line[colid_total_noise]
        if i < profile["n"]:
            line[colid_forward_add] = forward_add[0, i]
            line[colid_reverse_add] = reverse_add[0, i]
            line[colid_total_add] = forward_add[0, i] + reverse_add[0, i]
            line[colid_forward_corrected] += line[colid_forward_add]
            line[colid_reverse_corrected] += line[colid_reverse_add]
            line[colid_total_corrected] += line[colid_total_add]
            if "bgestimate" in profile["tool"][i]:
                line[colid_correction_flags] = "corrected_bgestimate"
            elif "bghomstats" in profile["tool"][i]:
                line[colid_correction_flags] = "corrected_bghomstats"
            elif "bgpredict" in profile["tool"][i]:
                line[colid_correction_flags] = "corrected_bgpredict"
            else:
                line[colid_correction_flags] = "corrected"
            line[colid_weight] = A[0, i]
        else:
            line[colid_correction_flags] = "corrected_as_background_only"
            line[colid_weight] = line[colid_total_corrected] / 100.

    # Add sequences that are in the profile but not in the sample.
    for i in range(profile["m"]):
        if profile["seqs"][i] in seqs:
            continue
        amount = forward_noise[0, i] + reverse_noise[0, i]
        if i < profile["n"]:
            amount += forward_add[0, i] + reverse_add[0, i]
        if amount > 0:
            line = [""] * len(column_names)
            line[colid_marker] = marker
            line[colid_sequence] = profile["seqs"][i]
            line[colid_forward] = 0
            line[colid_reverse] = 0
            line[colid_total] = 0
            line[colid_forward_noise] = forward_noise[0, i]
            line[colid_reverse_noise] = reverse_noise[0, i]
            line[colid_total_noise] = forward_noise[0, i] + reverse_noise[0, i]
            line[colid_forward_corrected] = -line[colid_forward_noise]
            line[colid_reverse_corrected] = -line[colid_reverse_noise]
            line[colid_total_corrected] = -line[colid_total_noise]
            if i < profile["n"]:
                line[colid_forward_add] = forward_add[0, i]
                line[colid_reverse_add] = reverse_add[0, i]
                line[colid_total_add] = forward_add[0, i] + reverse_add[0, i]
                line[colid_forward_corrected] += line[colid_forward_add]
                line[colid_reverse_corrected] += line[colid_reverse_add]
                line[colid_total_corrected] += line[colid_total_add]
                if "bgestimate" in profile["tool"][i]:
                    line[colid_correction_flags] = "corrected_bgestimate"
                elif "bghomstats" in profile["tool"][i]:
                    line[colid_correction_flags] = "corrected_bghomstats"
                elif "bgpredict" in profile["tool"][i]:
                    line[colid_correction_flags] = "corrected_bgpredict"
                else:
                    line[colid_correction_flags] = "corrected"
                line[colid_weight] = A[0, i]
            else:
                line[colid_forward_add] = 0
                line[colid_reverse_add] = 0
                line[colid_total_add] = 0
                line[colid_correction_flags] = "corrected_as_background_only"
                line[colid_weight] = line[colid_total_corrected] / 100.
            data.append(line)
#match_profile


def match_profiles(infile, outfile, profiles, library, seqformat):
    column_names, data = get_sample_data(
        infile, convert_to_raw=seqformat=="raw", library=library)
    colid_sequence = get_column_ids(column_names, "sequence")

    outfile.write("\t".join(column_names) + "\n")
    for marker in data:
        if marker in profiles:
            match_profile(column_names, data[marker], profiles[marker],
                          seqformat!="raw", library, marker)
        for line in data[marker]:
            if seqformat is not None and seqformat != "raw":
                line[colid_sequence] = ensure_sequence_format(
                    line[colid_sequence], seqformat, library=library,
                    marker=marker)
            outfile.write("\t".join(map(str, line)) + "\n")
#match_profiles


def add_arguments(parser):
    parser.add_argument('profiles', metavar="PROFILES",
        type=argparse.FileType('r'),
        help="file containing background noise profiles to match")
    add_input_output_args(parser, True, True, False)
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    add_sequence_format_args(parser)
#add_arguments


def run(args):
    # Import numpy now.
    global np
    import numpy as np

    gen = get_input_output_files(args, True, True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")

    # Read profiles once.
    profiles = load_profiles(args.profiles, args.library)
    if args.marker:
        profiles = {args.marker: profiles[args.marker]} \
                   if args.marker in profiles else {}

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError(
                "multiple input files for sample '%s' specified " % tag)
        infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "r")
        match_profiles(infile, outfile, profiles, args.library,
                       args.sequence_format)
        if infile != sys.stdin:
            infile.close()
#run
