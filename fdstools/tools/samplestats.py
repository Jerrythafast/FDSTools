#!/usr/bin/env python
"""
Compute various statistics for each sequence in the given sample data
file.

Adds the following columns to the input data (where 'X' is a placeholder
for 'forward', 'reverse', and 'total').  Some columns may be omitted
from the output if the input does not contain the required columns.

forward_pct: The percentage of reads of this sequence on the forward
strand.
forward_corrected_pct: The percentage of reads of this sequence on the
forward strand after noise correction.
X_correction_pct: The difference between the values of X_corrected and X
as a percentage of X.
X_mp: The number of X reads of this sequence as a percentage of the
total number of X reads of the corresponding marker.
X_noise_mp: The number of X_noise reads of this sequence as a percentage
of the total X_noise of the marker.
X_add_mp: The number of X_add reads of this sequence as a percentage of
the total X_add of the marker.
X_corrected_mp: The number of X_corrected reads of this sequence as a
percentage of the total number of X_corrected reads of the marker.
"""
import sys

from ..lib import ensure_sequence_format, add_sequence_format_args, \
                  add_input_output_args, get_input_output_files, get_column_ids

__version__ = "0.1dev"


def compute_stats(infile, outfile, library, seqformat):
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_name, colid_allele, colid_forward, colid_reverse, colid_total = \
        get_column_ids(column_names, "name", "allele", "forward", "reverse",
                       "total")
    (colid_forward_noise, colid_reverse_noise, colid_total_noise,
     colid_forward_add, colid_reverse_add, colid_total_add,
     colid_forward_corrected, colid_reverse_corrected,
     colid_total_corrected) = get_column_ids(column_names, "forward_noise",
        "reverse_noise", "total_noise", "forward_add", "reverse_add",
        "total_add", "forward_corrected", "reverse_corrected",
        "total_corrected", optional=True)

    # Add columns for which we have the required data.
    column_names.append("forward_pct")
    if colid_forward_corrected is not None:
        if colid_total_corrected is not None:
            column_names.append("forward_corrected_pct")
        column_names.append("forward_correction_pct")
    if colid_reverse_corrected is not None:
        column_names.append("reverse_correction_pct")
    if colid_total_corrected is not None:
        column_names.append("total_correction_pct")
    column_names.append("forward_mp")
    column_names.append("reverse_mp")
    column_names.append("total_mp")
    if colid_forward_noise is not None:
        column_names.append("forward_noise_mp")
    if colid_reverse_noise is not None:
        column_names.append("reverse_noise_mp")
    if colid_total_noise is not None:
        column_names.append("total_noise_mp")
    if colid_forward_add is not None:
        column_names.append("forward_add_mp")
    if colid_reverse_add is not None:
        column_names.append("reverse_add_mp")
    if colid_total_add is not None:
        column_names.append("total_add_mp")
    if colid_forward_corrected is not None:
        column_names.append("forward_corrected_mp")
    if colid_reverse_corrected is not None:
        column_names.append("reverse_corrected_mp")
    if colid_total_corrected is not None:
        column_names.append("total_corrected_mp")

    # Read data.
    data = {}
    for line in infile:
        cols = line.rstrip("\r\n").split("\t")
        marker = cols[colid_name]
        if seqformat is not None:
            cols[colid_allele] = ensure_sequence_format(
                cols[colid_allele], seqformat, library=library, marker=marker)
        cols[colid_forward] = int(cols[colid_forward])
        cols[colid_reverse] = int(cols[colid_reverse])
        cols[colid_total] = int(cols[colid_total])
        if colid_forward_corrected is not None:
            cols[colid_forward_corrected]=float(cols[colid_forward_corrected])
        if colid_reverse_corrected is not None:
            cols[colid_reverse_corrected]=float(cols[colid_reverse_corrected])
        if colid_total_corrected is not None:
            cols[colid_total_corrected] = float(cols[colid_total_corrected])
        if colid_forward_noise is not None:
            cols[colid_forward_noise] = float(cols[colid_forward_noise])
        if colid_reverse_noise is not None:
            cols[colid_reverse_noise] = float(cols[colid_reverse_noise])
        if colid_total_noise is not None:
            cols[colid_total_noise] = float(cols[colid_total_noise])
        if colid_forward_add is not None:
            cols[colid_forward_add] = float(cols[colid_forward_add])
        if colid_reverse_add is not None:
            cols[colid_reverse_add] = float(cols[colid_reverse_add])
        if colid_total_add is not None:
            cols[colid_total_add] = float(cols[colid_total_add])
        if marker not in data:
            data[marker] = []
        data[marker].append(cols)

    # Compute statistics.
    for marker in data:
        marker_forward = sum(row[colid_forward] for row in data[marker])
        marker_reverse = sum(row[colid_reverse] for row in data[marker])
        marker_total = sum(row[colid_total] for row in data[marker])
        if colid_forward_noise is not None:
            marker_forward_noise = sum(
                row[colid_forward_noise] for row in data[marker])
        if colid_reverse_noise is not None:
            marker_reverse_noise = sum(
                row[colid_reverse_noise] for row in data[marker])
        if colid_total_noise is not None:
            marker_total_noise = sum(
                row[colid_total_noise] for row in data[marker])
        if colid_forward_add is not None:
            marker_forward_add = sum(
                row[colid_forward_add] for row in data[marker])
        if colid_reverse_add is not None:
            marker_reverse_add = sum(
                row[colid_reverse_add] for row in data[marker])
        if colid_total_add is not None:
            marker_total_add = sum(
                row[colid_total_add] for row in data[marker])
        if colid_forward_corrected is not None:
            marker_forward_corrected = sum(
                row[colid_forward_corrected] for row in data[marker])
        if colid_reverse_corrected is not None:
            marker_reverse_corrected = sum(
                row[colid_reverse_corrected] for row in data[marker])
        if colid_total_corrected is not None:
            marker_total_corrected = sum(
                row[colid_total_corrected] for row in data[marker])
        for row in data[marker]:
            row.append(100.*row[colid_forward]/row[colid_total])
            if colid_forward_corrected is not None:
                if colid_total_corrected is not None:
                    row.append(100.*(
                        row[colid_forward_corrected]/row[colid_total_corrected]
                        if row[colid_total_corrected]
                        else row[colid_forward_corrected] > 0))
                row.append(
                    100.*row[colid_forward_corrected]/row[colid_forward]-100
                    if row[colid_forward]
                    else (row[colid_forward_corrected]>0)*200-100)
            if colid_reverse_corrected is not None:
                row.append(
                    100.*row[colid_reverse_corrected]/row[colid_reverse]-100
                    if row[colid_reverse]
                    else (row[colid_reverse_corrected]>0)*200-100)
            if colid_total_corrected is not None:
                row.append(
                    100.*row[colid_total_corrected]/row[colid_total]-100)
            row.append(100.*row[colid_forward]/marker_forward
                if marker_forward else 0)
            row.append(100.*row[colid_reverse]/marker_reverse
                if marker_reverse else 0)
            row.append(100.*row[colid_total]/marker_total
                if marker_total else 0)
            if colid_forward_noise is not None:
                row.append(100.*row[colid_forward_noise]/marker_forward_noise
                    if marker_forward_noise else 0)
            if colid_reverse_noise is not None:
                row.append(100.*row[colid_reverse_noise]/marker_reverse_noise
                    if marker_reverse_noise else 0)
            if colid_total_noise is not None:
                row.append(100.*row[colid_total_noise]/marker_total_noise
                    if marker_total_noise else 0)
            if colid_forward_add is not None:
                row.append(100.*row[colid_forward_add]/marker_forward_add
                    if marker_forward_add else 0)
            if colid_reverse_add is not None:
                row.append(100.*row[colid_reverse_add]/marker_reverse_add
                    if marker_reverse_add else 0)
            if colid_total_add is not None:
                row.append(100.*row[colid_total_add]/marker_total_add
                    if marker_total_add else 0)
            if colid_forward_corrected is not None:
                row.append(
                    100.*row[colid_forward_corrected]/marker_forward_corrected
                    if marker_forward_corrected else 0)
            if colid_reverse_corrected is not None:
                row.append(
                    100.*row[colid_reverse_corrected]/marker_reverse_corrected
                    if marker_reverse_corrected else 0)
            if colid_total_corrected is not None:
                row.append(
                    100.*row[colid_total_corrected]/marker_total_corrected
                    if marker_total_corrected else 0)

    # Write results.
    outfile.write("\t".join(column_names) + "\n")
    for marker in data:
        for line in data[marker]:
            outfile.write("\t".join(map(str, line)) + "\n")
#compute_stats


def add_arguments(parser):
    add_input_output_args(parser, True, True, False)
    add_sequence_format_args(parser)
#add_arguments


def run(args):
    gen = get_input_output_files(args, True, True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError(
                "multiple input files for sample '%s' specified " % tag)
        infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "r")
        compute_stats(infile, outfile, args.library, args.sequence_format)
        if infile != sys.stdin:
            infile.close()
#run
