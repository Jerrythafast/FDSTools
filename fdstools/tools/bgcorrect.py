#!/usr/bin/env python3

#
# Copyright (C) 2021 Jerry Hoogenboom
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
Match background noise profiles (obtained from e.g., bgestimate) to
samples.

Eleven new columns are added to the output giving, for each sequence,
the number of reads attributable to noise from other sequences (_noise
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

Finally, the weight column gives the number of times that the noise
profile of that allele fitted in the sample.
"""
import argparse
import sys
#import numpy as np  # Only imported when actually running this tool.

from errno import EPIPE

from ..lib.cli import add_sequence_format_args, add_input_output_args, get_input_output_files
from ..lib.io import get_column_ids
from ..lib.noise import load_profiles
from ..lib.seq import SEQ_SPECIAL_VALUES, ensure_sequence_format
from ..lib.util import nnls

__version__ = "1.1.0"


def get_sample_data(infile, convert_to_raw, library, combine_strands):
    """
    Read data from the given file handle, corresponding to a single
    sample, and fill a dict with all sequences in the sample.
    """
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return None, None  # Empty file.
    if not combine_strands:
        column_names.append("forward_noise")
        column_names.append("reverse_noise")
    column_names.append("total_noise")
    if not combine_strands:
        column_names.append("forward_add")
        column_names.append("reverse_add")
    column_names.append("total_add")
    if not combine_strands:
        column_names.append("forward_corrected")
        column_names.append("reverse_corrected")
    column_names.append("total_corrected")
    column_names.append("correction_flags")
    column_names.append("weight")
    colid_marker, colid_sequence, colid_forward, colid_reverse, colid_total = get_column_ids(
        column_names, "marker", "sequence", "forward", "reverse", "total")
    data = {}
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        marker = cols[colid_marker]
        if convert_to_raw:
            cols[colid_sequence] = ensure_sequence_format(
                cols[colid_sequence], "raw", library=library, marker=marker)
        cols[colid_forward] = int(cols[colid_forward])
        cols[colid_reverse] = int(cols[colid_reverse])
        cols[colid_total] = int(cols[colid_total])
        if not combine_strands:
            cols.append(0)
            cols.append(0)
        cols.append(0)
        if not combine_strands:
            cols.append(0)
            cols.append(0)
        cols.append(0)
        if not combine_strands:
            cols.append(cols[colid_forward])
            cols.append(cols[colid_reverse])
        cols.append(cols[colid_total])
        cols.append("not_corrected")
        cols.append(cols[colid_total] / 100)
        if marker not in data:
            data[marker] = []
        data[marker].append(cols)
    return column_names, data
#get_sample_data


def get_correction_tool_flags(profile_of_allele):
    return ",".join(("corrected_" + tool for tool in sorted(set(
        t for x in profile_of_allele.values() for t in x["tools"]))))
#get_correction_tool_flags


def get_noise_matrix(profile, seq_index, strand):
    Px = np.zeros((len(profile), len(seq_index)))
    np.fill_diagonal(Px, 100.)
    for allele in profile:
        allele_index = seq_index[allele]
        for sequence in profile[allele]:
            Px[allele_index, seq_index[sequence]] = profile[allele][sequence][strand]
    return Px
#get_noise_matrix


def match_profile(column_names, data, profile, seqformat, library, marker, combine_strands):
    (colid_marker, colid_sequence, colid_forward, colid_reverse, colid_total, colid_total_noise,
     colid_total_add, colid_total_corrected, colid_correction_flags,
     colid_weight) = get_column_ids(column_names,
        "marker", "sequence", "forward", "reverse", "total", "total_noise", "total_add",
        "total_corrected", "correction_flags", "weight")
    if not combine_strands:
        (colid_forward_noise, colid_reverse_noise, colid_forward_add, colid_reverse_add,
         colid_forward_corrected, colid_reverse_corrected) =  get_column_ids(column_names,
            "forward_noise", "reverse_noise", "forward_add", "reverse_add", "forward_corrected",
            "reverse_corrected")

    alleles = set(profile)
    other_seqs = set(seq for allele in profile for seq in profile[allele] if seq not in alleles)
    seq_index = {seq: i for i, seq in enumerate(
        sequence for subset in (alleles, other_seqs) for sequence in subset)}
    num_alleles = len(alleles)
    num_seqs = len(seq_index)

    # Enter profiles into P.
    if combine_strands:
        P = (get_noise_matrix(profile, seq_index, "total"),)
        C = (np.zeros((1, num_seqs)),)
    else:
        P = (get_noise_matrix(profile, seq_index, "forward"),
            get_noise_matrix(profile, seq_index, "reverse"))
        C = (np.zeros((1, num_seqs)), np.zeros((1, num_seqs)))

    # Enter sample into C.
    used_sequences = set()  # Sequences found in the sample.
    name2seq = {}  # Mapping names to raw sequences.
    for line in data:
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        if seqformat == "raw":
            # Already converted when reading the profile.
            sequence = line[colid_sequence]
        else:
            sequence = ensure_sequence_format(line[colid_sequence], "raw",
                                              library=library, marker=marker)
            name2seq[line[colid_sequence]] = sequence
        used_sequences.add(sequence)
        try:
            sequence_index = seq_index[sequence]
        except KeyError:
            # Note: Not adding any new sequences to the profile, since
            # they will just be zeroes and have no effect on the result.
            continue
        if combine_strands:
            C[0][0, sequence_index] = line[colid_total]
        else:
            C[0][0, sequence_index] = line[colid_forward]
            C[1][0, sequence_index] = line[colid_reverse]

    # Stop if this sample has no explicit data for this marker.
    if not used_sequences:
        return

    # Compute corrected read counts.
    A = nnls(np.hstack(P).T, np.hstack(C).T).T
    noises = []
    adds = []
    correcteds = []
    for Px, Cx in zip(P, C):
        np.fill_diagonal(Px, 0)
        noise = A @ Px
        add = A * Px.sum(1)
        corrected = Cx - noise
        corrected[:, :num_alleles] += add
        noises.append(noise)
        adds.append(add)
        correcteds.append(corrected)
    if combine_strands:
        colids_noise = (colid_total_noise,)
        colids_add = (colid_total_add,)
        colids_corr = (colid_total_corrected,)
    else:
        colids_noise = (colid_forward_noise, colid_reverse_noise, colid_total_noise)
        colids_add = (colid_forward_add, colid_reverse_add, colid_total_add)
        colids_corr = (colid_forward_corrected, colid_reverse_corrected, colid_total_corrected)
        noises.append(sum(noises))
        adds.append(sum(adds))
        correcteds.append(sum(correcteds))

    # Round values to 3 decimal positions.
    A.round(3, A)
    for noise, add, corrected in zip(noises, adds, correcteds):
        noise.round(3, noise)
        add.round(3, add)
        corrected.round(3, corrected)

    for line in data:
        if line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        sequence = name2seq.get(line[colid_sequence], line[colid_sequence])
        try:
            i = seq_index[sequence]
        except KeyError:
            line[colid_correction_flags] = "not_in_ref_db"
            continue
        for colid_noise, colid_corrected, noise, corrected in zip(colids_noise, colids_corr, noises, correcteds):
            line[colid_noise] = noise[0, i]
            line[colid_corrected] = corrected[0, i]
        if i < num_alleles:
            for colid_add, colid_corrected, add in zip(colids_add, colids_corr, adds):
                line[colid_add] = add[0, i]
            line[colid_correction_flags] = get_correction_tool_flags(profile[sequence])
            line[colid_weight] = A[0, i]
        else:
            line[colid_correction_flags] = "corrected_as_background_only"
            line[colid_weight] = line[colid_total_corrected] / 100

    # Add sequences that are in the profile but not in the sample.
    for sequence, i in seq_index.items():
        if sequence in used_sequences:
            continue
        if correcteds[-1][0, i] > 0:
            line = [""] * len(column_names)
            line[colid_marker] = marker
            line[colid_sequence] = sequence if seqformat in (None, "raw") \
                else ensure_sequence_format(sequence, seqformat, library=library, marker=marker)
            line[colid_forward] = 0
            line[colid_reverse] = 0
            line[colid_total] = 0
            for colid_noise, colid_corrected, noise, corrected in zip(colids_noise, colids_corr, noises, correcteds):
                line[colid_noise] = noise[0, i]
                line[colid_corrected] = corrected[0, i]
            if i < num_alleles:
                for colid_add, colid_corrected, add in zip(colids_add, colids_corr, adds):
                    line[colid_add] = add[0, i]
                line[colid_correction_flags] = get_correction_tool_flags(profile[sequence])
                line[colid_weight] = A[0, i]
            else:
                line[colid_forward_add] = 0
                line[colid_reverse_add] = 0
                line[colid_total_add] = 0
                line[colid_correction_flags] = "corrected_as_background_only"
                line[colid_weight] = line[colid_total_corrected] / 100
            data.append(line)
#match_profile


def match_profiles(infile, outfile, profiles, library, seqformat, combine_strands):
    column_names, data = get_sample_data(infile, seqformat=="raw", library, combine_strands)
    if column_names is None:
        return  # Empty file.
    colid_sequence = get_column_ids(column_names, "sequence")

    outfile.write("\t".join(column_names) + "\n")
    for marker in data:
        if marker in profiles:
            match_profile(column_names, data[marker], profiles[marker],
                          seqformat, library, marker, combine_strands)
        for line in data[marker]:
            if seqformat is not None and seqformat != "raw":
                line[colid_sequence] = ensure_sequence_format(
                    line[colid_sequence], seqformat, library=library, marker=marker)
            outfile.write("\t".join(map(str, line)) + "\n")
#match_profiles


def add_arguments(parser):
    parser.add_argument("profiles", metavar="PROFILES",
        type=argparse.FileType("tr", encoding="UTF-8"),
        help="file containing background noise profiles to match")
    add_input_output_args(parser, single_in=True, batch_support=True, report_out=False)
    parser.add_argument("-C", "--combine-strands", action="store_true",
        help="if specified, stutter noise correction will be done on the total number of reads, "
             "instead of separately for either strand")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-M", "--marker", metavar="MARKER", help="work only on MARKER")
    add_sequence_format_args(parser)
#add_arguments


def run(args):
    # Import numpy now.
    global np
    import numpy as np

    gen = get_input_output_files(args, single_in=True, batch_support=True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output of another program")

    # Read profiles once.
    profiles = load_profiles(args.profiles, args.library)
    if args.marker:
        profiles = {args.marker: profiles[args.marker]} if args.marker in profiles else {}

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError("multiple input files for sample '%s' specified " % tag)
        try:
            infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "rt", encoding="UTF-8")
            match_profiles(infile, outfile, profiles, args.library, args.sequence_format,
                args.combine_strands)
            if infile != sys.stdin:
                infile.close()
        except IOError as e:
            if e.errno == EPIPE:
                continue
            raise
#run
