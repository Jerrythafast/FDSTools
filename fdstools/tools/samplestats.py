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
Compute various statistics for each sequence in the given sample data
file and perform threshold-based allele calling.

Updates the 'flags' column (or adds it, if it was not present in the
input data) to include 'allele' for all sequences that meet various
allele calling thresholds.

Adds the following columns to the input data.  Some columns may be
omitted from the output if the input does not contain the required
columns.  In the column names below, 'X' is a placeholder for 'forward',
'reverse', and 'total', which refers to the strand of DNA for which the
statistics are calculated.  'Y' is a placeholder for 'corrected'
(statistics calculated on data after noise correction by e.g.,
BGCorrect), 'noise' (statistics calculated on the number of reads
attributed to noise), and 'add' (statistics calculated on the number
of reads recovered through noise correction).  Wherever the 'Y' part of
the column name is omitted, the values in the column are computed on
data prior to noise correction.

X_Y: The number of Y reads of this sequence on the X strand (this column
is not added by Samplestats, but should be present in the input).
X_Y_mp_sum: The value of X_Y, as a percentage of the sum of the X_Y of
the marker.
X_Y_mp_max: The value of X_Y, as a percentage of the maximum X_Y of the
marker.
forward_Y_pct: The number of Y reads on the forward strand, as a
percentage of the total number of Y reads of this sequence.
X_correction_pct: The difference between the values of X_corrected and
X, as a percentage of the value of X.
X_removed_pct: The value of X_noise, as a percentage of the value of X.
X_added_pct: The value of X_add, as a percentage of the value of X.
X_recovery: The value of X_add, as a percentage of the value of
X_corrected.
"""
import sys

from ..lib import add_sequence_format_args, add_input_output_args, \
                  get_input_output_files, get_column_ids

__version__ = "1.1.0"


# Default values for parameters are specified below.

# Default minimum number of reads to mark as allele.
# This value can be overridden by the -n command line option.
_DEF_MIN_READS = 30

# Default minimum number of reads per strand to mark as allele.
# This value can be overridden by the -b command line option.
_DEF_MIN_PER_STRAND = 1

# Default minimum percentage of reads w.r.t. the highest allele of the
# marker to mark as allele.
# This value can be overridden by the -m command line option.
_DEF_MIN_PCT_OF_MAX = 2.

# Default minimum percentage of reads w.r.t. the marker's total number
# of reads to mark as allele.
# This value can be overridden by the -p command line option.
_DEF_MIN_PCT_OF_SUM = 1.5

# Default minimum percentage of correction to mark as allele.
# This value can be overridden by the -c command line option.
_DEF_MIN_CORRECTION = 0

# Default minimum number of recovered reads as a percentage of the
# number of reads after correction to mark as allele.
# This value can be overridden by the -r command line option.
_DEF_MIN_RECOVERY = 0

# Default minimum number of reads for filtering.
# This value can be overridden by the -N command line option.
_DEF_MIN_READS_FILT = 1

# Default minimum number of reads per strand for filtering.
# This value can be overridden by the -B command line option.
_DEF_MIN_PER_STRAND_FILT = 1

# Default minimum percentage of reads w.r.t. the highest allele of the
# marker for filtering.
# This value can be overridden by the -M command line option.
_DEF_MIN_PCT_OF_MAX_FILT = 0.

# Default minimum percentage of reads w.r.t. the marker's total number
# of reads for filtering.
# This value can be overridden by the -P command line option.
_DEF_MIN_PCT_OF_SUM_FILT = 0.

# Default minimum percentage of correction for filtering.
# This value can be overridden by the -C command line option.
_DEF_MIN_CORRECTION_FILT = 0

# Default minimum number of recovered reads as a percentage of the
# number of reads after correction for filtering.
# This value can be overridden by the -R command line option.
_DEF_MIN_RECOVERY_FILT = 0

COLUMN_ORDER = [
    "total_corrected",
    "total_corrected_mp_sum",
    "total_corrected_mp_max",
    "total_correction_pct",
    "total_recovery",
    "weight",
    "forward_corrected_pct",
    "forward_corrected",
    "forward_corrected_mp_sum",
    "forward_corrected_mp_max",
    "forward_correction_pct",
    "forward_recovery",
    "reverse_corrected",
    "reverse_corrected_mp_sum",
    "reverse_corrected_mp_max",
    "reverse_correction_pct",
    "reverse_recovery",

    "total",
    "total_mp_sum",
    "total_mp_max",
    "forward_pct",
    "forward",
    "forward_mp_sum",
    "forward_mp_max",
    "reverse",
    "reverse_mp_sum",
    "reverse_mp_max",

    "total_noise",
    "total_noise_mp_sum",
    "total_noise_mp_max",
    "total_removed_pct",
    "forward_noise_pct",
    "forward_noise",
    "forward_noise_mp_sum",
    "forward_noise_mp_max",
    "forward_removed_pct",
    "reverse_noise",
    "reverse_noise_mp_sum",
    "reverse_noise_mp_max",
    "reverse_removed_pct",

    "total_add",
    "total_add_mp_sum",
    "total_add_mp_max",
    "total_added_pct",
    "forward_add_pct",
    "forward_add",
    "forward_add_mp_sum",
    "forward_add_mp_max",
    "forward_added_pct",
    "reverse_add",
    "reverse_add_mp_sum",
    "reverse_add_mp_max",
    "reverse_added_pct"
]


def compute_stats(infile, outfile, min_reads,
                  min_per_strand, min_pct_of_max, min_pct_of_sum,
                  min_correction, min_recovery, filter_action, filter_absolute,
                  min_reads_filt, min_per_strand_filt, min_pct_of_max_filt,
                  min_pct_of_sum_filt, min_correction_filt, min_recovery_filt):
    # Check presence of required columns.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    get_column_ids(column_names, "marker", "sequence", "forward", "reverse",
        "total")
    if "flags" not in column_names:
        column_names.append("flags")

    # Add columns for which we have the required data.
    if "total_corrected" in column_names:
        column_names.append("total_corrected_mp_sum")
        column_names.append("total_corrected_mp_max")
        column_names.append("total_correction_pct")
        if "forward_corrected" in column_names:
            column_names.append("forward_corrected_pct")
        if "total_add" in column_names:
            column_names.append("total_recovery")
    if "forward_corrected" in column_names:
        column_names.append("forward_corrected_mp_sum")
        column_names.append("forward_corrected_mp_max")
        column_names.append("forward_correction_pct")
        if "forward_add" in column_names:
            column_names.append("forward_recovery")
    if "reverse_corrected" in column_names:
        column_names.append("reverse_corrected_mp_sum")
        column_names.append("reverse_corrected_mp_max")
        column_names.append("reverse_correction_pct")
        if "reverse_add" in column_names:
            column_names.append("reverse_recovery")
    column_names.append("total_mp_sum")
    column_names.append("total_mp_max")
    column_names.append("forward_pct")
    column_names.append("forward_mp_sum")
    column_names.append("forward_mp_max")
    column_names.append("reverse_mp_sum")
    column_names.append("reverse_mp_max")
    if "total_noise" in column_names:
        column_names.append("total_noise_mp_sum")
        column_names.append("total_noise_mp_max")
        column_names.append("total_removed_pct")
        if "forward_noise" in column_names:
            column_names.append("forward_noise_pct")
    if "forward_noise" in column_names:
        column_names.append("forward_noise_mp_sum")
        column_names.append("forward_noise_mp_max")
        column_names.append("forward_removed_pct")
    if "reverse_noise" in column_names:
        column_names.append("reverse_noise_mp_sum")
        column_names.append("reverse_noise_mp_max")
        column_names.append("reverse_removed_pct")
    if "total_add" in column_names:
        column_names.append("total_add_mp_sum")
        column_names.append("total_add_mp_max")
        column_names.append("total_added_pct")
        if "forward_add" in column_names:
            column_names.append("forward_add_pct")
    if "forward_add" in column_names:
        column_names.append("forward_add_mp_sum")
        column_names.append("forward_add_mp_max")
        column_names.append("forward_added_pct")
    if "reverse_add" in column_names:
        column_names.append("reverse_add_mp_sum")
        column_names.append("reverse_add_mp_max")
        column_names.append("reverse_added_pct")

    # Build a column number lookup dictionary.
    ci = {column_names[i]: i for i in range(len(column_names))}

    # Read data.
    data = {}
    for line in infile:
        row = line.rstrip("\r\n").split("\t")
        marker = row[ci["marker"]]
        row[ci["forward"]] = int(row[ci["forward"]])
        row[ci["reverse"]] = int(row[ci["reverse"]])
        row[ci["total"]] = int(row[ci["total"]])
        for i in (ci[column] for column in (
                "forward_corrected", "reverse_corrected", "total_corrected",
                "forward_noise", "reverse_noise", "total_noise", "weight",
                "forward_add", "reverse_add", "total_add") if column in ci):
            row[i] = float(row[i])
        if len(row) == ci["flags"]:
            row.append([])
        else:
            row[ci["flags"]] = map(str.strip, row[ci["flags"]].split(","))
        if marker not in data:
            data[marker] = []
        data[marker].append(row)

    # Compute statistics.
    if filter_action != "off":
        filtered = {marker: [] for marker in data}
    for marker in data:
        if "total_corrected" in ci:
            marker_total_corrected_sum = sum(
                row[ci["total_corrected"]] for row in data[marker])
            marker_total_corrected_max = max(
                row[ci["total_corrected"]] for row in data[marker])
        if "forward_corrected" in ci:
            marker_forward_corrected_sum = sum(
                row[ci["forward_corrected"]] for row in data[marker])
            marker_forward_corrected_max = max(
                row[ci["forward_corrected"]] for row in data[marker])
        if "reverse_corrected" in ci:
            marker_reverse_corrected_sum = sum(
                row[ci["reverse_corrected"]] for row in data[marker])
            marker_reverse_corrected_max = max(
                row[ci["reverse_corrected"]] for row in data[marker])
        marker_total_sum = sum(row[ci["total"]] for row in data[marker])
        marker_total_max = max(row[ci["total"]] for row in data[marker])
        marker_forward_sum = sum(row[ci["forward"]] for row in data[marker])
        marker_forward_max = max(row[ci["forward"]] for row in data[marker])
        marker_reverse_sum = sum(row[ci["reverse"]] for row in data[marker])
        marker_reverse_max = max(row[ci["reverse"]] for row in data[marker])
        if "total_noise" in ci:
            marker_total_noise_sum = sum(
                row[ci["total_noise"]] for row in data[marker])
            marker_total_noise_max = max(
                row[ci["total_noise"]] for row in data[marker])
        if "forward_noise" in ci:
            marker_forward_noise_sum = sum(
                row[ci["forward_noise"]] for row in data[marker])
            marker_forward_noise_max = max(
                row[ci["forward_noise"]] for row in data[marker])
        if "reverse_noise" in ci:
            marker_reverse_noise_sum = sum(
                row[ci["reverse_noise"]] for row in data[marker])
            marker_reverse_noise_max = max(
                row[ci["reverse_noise"]] for row in data[marker])
        if "total_add" in ci:
            marker_total_add_sum = sum(
                row[ci["total_add"]] for row in data[marker])
            marker_total_add_max = max(
                row[ci["total_add"]] for row in data[marker])
        if "forward_add" in ci:
            marker_forward_add_sum = sum(
                row[ci["forward_add"]] for row in data[marker])
            marker_forward_add_max = max(
                row[ci["forward_add"]] for row in data[marker])
        if "reverse_add" in ci:
            marker_reverse_add_sum = sum(
                row[ci["reverse_add"]] for row in data[marker])
            marker_reverse_add_max = max(
                row[ci["reverse_add"]] for row in data[marker])
        for row in data[marker]:
            if "total_corrected" in ci:
                row.append(100.*row[ci["total_corrected"]] /
                    marker_total_corrected_sum
                    if marker_total_corrected_sum else 0)
                row.append(100.*row[ci["total_corrected"]] /
                    marker_total_corrected_max
                    if marker_total_corrected_max else 0)
                row.append(
                    100.*row[ci["total_corrected"]]/row[ci["total"]]-100
                    if row[ci["total"]]
                    else ((row[ci["total_corrected"]]>0)*200-100
                        if row[ci["total_corrected"]] else 0))
                if "forward_corrected" in ci:
                    row.append(100.*(
                        row[ci["forward_corrected"]]/row[ci["total_corrected"]]
                        if row[ci["total_corrected"]]
                        else row[ci["forward_corrected"]] > 0))
                if "total_add" in ci:
                    row.append(100.*row[ci["total_add"]] /
                        row[ci["total_corrected"]]
                        if row[ci["total_corrected"]] else 0)
            if "forward_corrected" in ci:
                row.append(100.*row[ci["forward_corrected"]] /
                    marker_forward_corrected_sum
                    if marker_forward_corrected_sum else 0)
                row.append(100.*row[ci["forward_corrected"]] /
                    marker_forward_corrected_max
                    if marker_forward_corrected_max else 0)
                row.append(
                    100.*row[ci["forward_corrected"]]/row[ci["forward"]]-100
                    if row[ci["forward"]]
                    else ((row[ci["forward_corrected"]]>0)*200-100
                        if row[ci["forward_corrected"]] else 0))
                if "forward_add" in ci:
                    row.append(100.*row[ci["forward_add"]] /
                        row[ci["forward_corrected"]]
                        if row[ci["forward_corrected"]] else 0)
            if "reverse_corrected" in ci:
                row.append(100.*row[ci["reverse_corrected"]] /
                    marker_reverse_corrected_sum
                    if marker_reverse_corrected_sum else 0)
                row.append(100.*row[ci["reverse_corrected"]] /
                    marker_reverse_corrected_max
                    if marker_reverse_corrected_max else 0)
                row.append(
                    100.*row[ci["reverse_corrected"]]/row[ci["reverse"]]-100
                    if row[ci["reverse"]]
                    else ((row[ci["reverse_corrected"]]>0)*200-100
                        if row[ci["reverse_corrected"]] else 0))
                if "reverse_add" in ci:
                    row.append(100.*row[ci["reverse_add"]] /
                        row[ci["reverse_corrected"]]
                        if row[ci["reverse_corrected"]] else 0)
            row.append(100.*row[ci["total"]]/marker_total_sum
                if marker_total_sum else 0)
            row.append(100.*row[ci["total"]]/marker_total_max
                if marker_total_max else 0)
            row.append(100.*(1.*row[ci["forward"]]/row[ci["total"]]
                if row[ci["total"]] else row[ci["forward"]] > 0))
            row.append(100.*row[ci["forward"]]/marker_forward_sum
                if marker_forward_sum else 0)
            row.append(100.*row[ci["forward"]]/marker_forward_max
                if marker_forward_max else 0)
            row.append(100.*row[ci["reverse"]]/marker_reverse_sum
                if marker_reverse_sum else 0)
            row.append(100.*row[ci["reverse"]]/marker_reverse_max
                if marker_reverse_max else 0)
            if "total_noise" in ci:
                row.append(100.*row[ci["total_noise"]]/marker_total_noise_sum
                    if marker_total_noise_sum else 0)
                row.append(100.*row[ci["total_noise"]]/marker_total_noise_max
                    if marker_total_noise_max else 0)
                row.append(100.*row[ci["total_noise"]]/row[ci["total"]]
                    if row[ci["total"]] else 0)
                if "forward_noise" in ci:
                    row.append(100.*(
                        row[ci["forward_noise"]]/row[ci["total_noise"]]
                        if row[ci["total_noise"]]
                        else row[ci["forward_noise"]] > 0))
            if "forward_noise" in ci:
                row.append(100.*row[ci["forward_noise"]]/
                    marker_forward_noise_sum
                    if marker_forward_noise_sum else 0)
                row.append(100.*row[ci["forward_noise"]]/
                    marker_forward_noise_max
                    if marker_forward_noise_max else 0)
                row.append(100.*row[ci["forward_noise"]]/row[ci["forward"]]
                    if row[ci["forward"]] else 0)
            if "reverse_noise" in ci:
                row.append(100.*row[ci["reverse_noise"]]/
                    marker_reverse_noise_sum
                    if marker_reverse_noise_sum else 0)
                row.append(100.*row[ci["reverse_noise"]]/
                    marker_reverse_noise_max
                    if marker_reverse_noise_max else 0)
                row.append(100.*row[ci["reverse_noise"]]/row[ci["reverse"]]
                    if row[ci["reverse"]] else 0)
            if "total_add" in ci:
                row.append(100.*row[ci["total_add"]]/marker_total_add_sum
                    if marker_total_add_sum else 0)
                row.append(100.*row[ci["total_add"]]/marker_total_add_max
                    if marker_total_add_max else 0)
                row.append(100.*row[ci["total_add"]]/row[ci["total"]]
                    if row[ci["total"]] else 0)
                if "forward_add" in ci:
                    row.append(100.*(
                        row[ci["forward_add"]]/row[ci["total_add"]]
                        if row[ci["total_add"]]
                        else row[ci["forward_add"]] > 0))
            if "forward_add" in ci:
                row.append(100.*row[ci["forward_add"]]/marker_forward_add_sum
                    if marker_forward_add_sum else 0)
                row.append(100.*row[ci["forward_add"]]/marker_forward_add_max
                    if marker_forward_add_max else 0)
                row.append(100.*row[ci["forward_add"]]/row[ci["forward"]]
                    if row[ci["forward"]] else 0)
            if "reverse_add" in ci:
                row.append(100.*row[ci["reverse_add"]]/marker_reverse_add_sum
                    if marker_reverse_add_sum else 0)
                row.append(100.*row[ci["reverse_add"]]/marker_reverse_add_max
                    if marker_reverse_add_max else 0)
                row.append(100.*row[ci["reverse_add"]]/row[ci["reverse"]]
                    if row[ci["reverse"]] else 0)

            # The 'No data' lines are fine like this.
            if row[ci["sequence"]] == "No data":
                row[ci["flags"]] = ",".join(row[ci["flags"]])
                continue

            # Get the values we will filter on.
            total_added = row[ci["total"]] if "total_corrected" not in ci \
                else row[ci["total_corrected"]]
            pct_of_sum = row[ci["total_mp_sum"]] if "total_corrected_mp_sum" \
                not in ci else row[ci["total_corrected_mp_sum"]]
            pct_of_max = row[ci["total_mp_max"]] if "total_corrected_mp_max" \
                not in ci else row[ci["total_corrected_mp_max"]]
            correction = 0 if "total_correction_pct" not in ci \
                else row[ci["total_correction_pct"]]
            recovery = 0 if "total_recovery" not in ci \
                else row[ci["total_recovery"]]
            strands = [
                row[ci["forward"]] if "forward_corrected" not in ci
                    else row[ci["forward_corrected"]],
                row[ci["reverse"]] if "reverse_corrected" not in ci
                    else row[ci["reverse_corrected"]]]
            fn = abs if filter_absolute else lambda x: x

            # Check if this sequence should be filtered out.
            # Always filter/combine existing 'Other sequences'.
            if filter_action != "off" and (
                    row[ci["sequence"]] == "Other sequences" or (
                    fn(total_added) < min_reads_filt or
                    fn(pct_of_max) < min_pct_of_max_filt or
                    fn(pct_of_sum) < min_pct_of_sum_filt or
                    (correction < min_correction_filt and
                    recovery < min_recovery_filt) or
                    min(map(fn, strands)) < min_per_strand_filt)):
                filtered[marker].append(row)

            # Check if this sequence is an allele.
            elif (row[ci["sequence"]] != "Other sequences" and
                    total_added >= min_reads and
                    pct_of_max >= min_pct_of_max and
                    pct_of_sum >= min_pct_of_sum and
                    (correction >= min_correction or
                    recovery >= min_recovery) and
                    min(strands) >= min_per_strand):
                row[ci["flags"]].append("allele")
            row[ci["flags"]] = ",".join(row[ci["flags"]])

    # Reorder columns.
    new_order = {}
    for i in range(len(COLUMN_ORDER)-1, -1, -1):
        if COLUMN_ORDER[i] in ci:
            new_order[len(ci) - len(new_order) - 1] = ci[COLUMN_ORDER[i]]
    for i in range(len(column_names)-1, -1, -1):
        if column_names[i] not in COLUMN_ORDER:
            new_order[len(ci) - len(new_order) - 1] = i

    # Write results.
    outfile.write("\t".join(
        column_names[new_order[i]] for i in range(len(column_names))) + "\n")
    for marker in sorted(data):
        if filter_action == "combine":
            have_combined = False
            combined = [""] * len(column_names)
            combined[ci["marker"]] = marker
            combined[ci["sequence"]] = "Other sequences"
            for i in (ci[column] for column in COLUMN_ORDER if column in ci):
                # Set known numeric columns to 0.
                combined[i] = 0

        for row in data[marker]:
            if filter_action == "combine" and row in filtered[marker]:
                have_combined = True
                if "total_corrected" in ci:
                    combined[ci["total_corrected"]] += \
                        row[ci["total_corrected"]]
                    combined[ci["total_corrected_mp_sum"]] += \
                        row[ci["total_corrected_mp_sum"]]
                    combined[ci["total_corrected_mp_max"]] += \
                        row[ci["total_corrected_mp_max"]]
                if "forward_corrected" in ci:
                    combined[ci["forward_corrected"]] += \
                        row[ci["forward_corrected"]]
                    combined[ci["forward_corrected_mp_sum"]] += \
                        row[ci["forward_corrected_mp_sum"]]
                    combined[ci["forward_corrected_mp_max"]] += \
                        row[ci["forward_corrected_mp_max"]]
                if "reverse_corrected" in ci:
                    combined[ci["reverse_corrected"]] += \
                        row[ci["reverse_corrected"]]
                    combined[ci["reverse_corrected_mp_sum"]] += \
                        row[ci["reverse_corrected_mp_sum"]]
                    combined[ci["reverse_corrected_mp_max"]] += \
                        row[ci["reverse_corrected_mp_max"]]
                combined[ci["total"]] += row[ci["total"]]
                combined[ci["total_mp_sum"]] += row[ci["total_mp_sum"]]
                combined[ci["total_mp_max"]] += row[ci["total_mp_max"]]
                combined[ci["forward"]] += row[ci["forward"]]
                combined[ci["forward_mp_sum"]] += row[ci["forward_mp_sum"]]
                combined[ci["forward_mp_max"]] += row[ci["forward_mp_max"]]
                combined[ci["reverse"]] += row[ci["reverse"]]
                combined[ci["reverse_mp_sum"]] += row[ci["reverse_mp_sum"]]
                combined[ci["reverse_mp_max"]] += row[ci["reverse_mp_max"]]
                if "total_noise" in ci:
                    combined[ci["total_noise"]] += row[ci["total_noise"]]
                    combined[ci["total_noise_mp_sum"]] += \
                        row[ci["total_noise_mp_sum"]]
                    combined[ci["total_noise_mp_max"]] += \
                        row[ci["total_noise_mp_max"]]
                if "forward_noise" in ci:
                    combined[ci["forward_noise"]] += row[ci["forward_noise"]]
                    combined[ci["forward_noise_mp_sum"]] += \
                        row[ci["forward_noise_mp_sum"]]
                    combined[ci["forward_noise_mp_max"]] += \
                        row[ci["forward_noise_mp_max"]]
                if "reverse_noise" in ci:
                    combined[ci["reverse_noise"]] += row[ci["reverse_noise"]]
                    combined[ci["reverse_noise_mp_sum"]] += \
                        row[ci["reverse_noise_mp_sum"]]
                    combined[ci["reverse_noise_mp_max"]] += \
                        row[ci["reverse_noise_mp_max"]]
                if "total_add" in ci:
                    combined[ci["total_add"]] += row[ci["total_add"]]
                    combined[ci["total_add_mp_sum"]] += \
                        row[ci["total_add_mp_sum"]]
                    combined[ci["total_add_mp_max"]] += \
                        row[ci["total_add_mp_max"]]
                if "forward_add" in ci:
                    combined[ci["forward_add"]] += row[ci["forward_add"]]
                    combined[ci["forward_add_mp_sum"]] += \
                        row[ci["forward_add_mp_sum"]]
                    combined[ci["forward_add_mp_max"]] += \
                        row[ci["forward_add_mp_max"]]
                if "reverse_add" in ci:
                    combined[ci["reverse_add"]] += row[ci["reverse_add"]]
                    combined[ci["reverse_add_mp_sum"]] += \
                        row[ci["reverse_add_mp_sum"]]
                    combined[ci["reverse_add_mp_max"]] += \
                        row[ci["reverse_add_mp_max"]]
                if "weight" in ci:
                    combined[ci["weight"]] += row[ci["weight"]]
            elif filter_action == "off" or row not in filtered[marker]:
                for i in (ci[col] for col in COLUMN_ORDER if col in ci
                        and col not in ("total", "forward", "reverse")):
                    if row[i] >= 10000:
                        row[i] = "%#.5g" % row[i]
                    elif row[i] >= 1000:
                        row[i] = "%#.4g" % row[i]
                    else:
                        row[i] = "%#.3g" % row[i]
                outfile.write("\t".join(map(str,
                    (row[new_order[i]] for i in range(len(row))))) + "\n")

        # Add combined row for this marker.
        if filter_action == "combine" and have_combined:
            if "total_corrected" in ci:
                combined[ci["total_correction_pct"]] = (
                    100.*combined[ci["total_corrected"]]/
                        combined[ci["total"]]-100
                    if combined[ci["total"]]
                    else ((combined[ci["total_corrected"]]>0)*200-100
                        if combined[ci["total_corrected"]] else 0))
                if "forward_corrected" in ci:
                    combined[ci["forward_corrected_pct"]] = 100.*(
                        combined[ci["forward_corrected"]]/
                            combined[ci["total_corrected"]]
                        if combined[ci["total_corrected"]]
                        else combined[ci["forward_corrected"]] > 0)
                if "total_add" in ci:
                    combined[ci["total_recovery"]] = (
                        100.*combined[ci["total_add"]]/
                            combined[ci["total_corrected"]]
                            if combined[ci["total_corrected"]] else 0)
            if "forward_corrected" in ci:
                combined[ci["forward_correction_pct"]] = (
                    100.*combined[ci["forward_corrected"]]/
                        combined[ci["forward"]]-100
                    if combined[ci["forward"]]
                    else ((combined[ci["forward_corrected"]]>0)*200-100
                        if combined[ci["forward_corrected"]] else 0))
                if "forward_add" in ci:
                    combined[ci["forward_recovery"]] = (
                        100.*combined[ci["forward_add"]]/
                            combined[ci["forward_corrected"]]
                            if combined[ci["forward_corrected"]] else 0)
            if "reverse_corrected" in ci:
                combined[ci["reverse_correction_pct"]] = (
                    100.*combined[ci["reverse_corrected"]]/
                        combined[ci["reverse"]]-100
                    if combined[ci["reverse"]]
                    else ((combined[ci["reverse_corrected"]]>0)*200-100
                        if combined[ci["reverse_corrected"]] else 0))
                if "reverse_add" in ci:
                    combined[ci["reverse_recovery"]] = (
                        100.*combined[ci["reverse_add"]]/
                            combined[ci["reverse_corrected"]]
                            if combined[ci["reverse_corrected"]] else 0)
            combined[ci["forward_pct"]] = 100.*(
                1.*combined[ci["forward"]]/combined[ci["total"]]
                if combined[ci["total"]] else combined[ci["forward"]] > 0)
            if "total_noise" in ci:
                combined[ci["total_removed_pct"]] = (
                    100.*combined[ci["total_noise"]]/combined[ci["total"]]
                    if combined[ci["total"]] else 0)
                if "forward_noise" in ci:
                    combined[ci["forward_noise_pct"]] = 100.*(
                        combined[ci["forward_noise"]]/
                            combined[ci["total_noise"]]
                        if combined[ci["total_noise"]]
                        else combined[ci["forward_noise"]] > 0)
            if "forward_noise" in ci:
                combined[ci["forward_removed_pct"]] = (
                    100.*combined[ci["forward_noise"]]/combined[ci["forward"]]
                    if combined[ci["forward"]] else 0)
            if "reverse_noise" in ci:
                combined[ci["reverse_removed_pct"]] = (
                    100.*combined[ci["reverse_noise"]]/combined[ci["reverse"]]
                    if combined[ci["reverse"]] else 0)
            if "total_add" in ci:
                combined[ci["total_added_pct"]] = (
                    100.*combined[ci["total_add"]]/combined[ci["total"]]
                    if combined[ci["total"]] else 0)
                if "forward_add" in ci:
                    combined[ci["forward_add_pct"]] = 100.*(
                        combined[ci["forward_add"]]/combined[ci["total_add"]]
                        if combined[ci["total_add"]]
                        else combined[ci["forward_add"]] > 0)
            if "forward_add" in ci:
                combined[ci["forward_added_pct"]] = (
                    100.*combined[ci["forward_add"]]/combined[ci["forward"]]
                    if combined[ci["forward"]] else 0)
            if "reverse_add" in ci:
                combined[ci["reverse_added_pct"]] = (
                    100.*combined[ci["reverse_add"]]/combined[ci["reverse"]]
                    if combined[ci["reverse"]] else 0)

            for i in (ci[column] for column in COLUMN_ORDER if column in ci
                    and column not in ("total", "forward", "reverse")):
                if combined[i] >= 10000:
                    combined[i] = "%#.5g" % combined[i]
                elif combined[i] >= 1000:
                    combined[i] = "%#.4g" % combined[i]
                else:
                    combined[i] = "%#.3g" % combined[i]
            outfile.write("\t".join(map(str,
                (combined[new_order[i]] for i in range(len(combined))))) +"\n")
#compute_stats


def add_arguments(parser):
    add_input_output_args(parser, True, True, False)
    intergroup = parser.add_argument_group("interpretation options",
        "sequences that match the -c or -y option (or both) and all of the "
        "other settings are marked as 'allele'")
    intergroup.add_argument('-n', '--min-reads', metavar="N", type=float,
        default=_DEF_MIN_READS,
        help="the minimum number of reads (default: %(default)s)")
    intergroup.add_argument('-b', '--min-per-strand', metavar="N", type=float,
        default=_DEF_MIN_PER_STRAND,
        help="the minimum number of reads in both orientations (default: "
             "%(default)s)")
    intergroup.add_argument('-m', '--min-pct-of-max', metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_MAX,
        help="the minimum percentage of reads w.r.t. the highest allele of "
             "the marker (default: %(default)s)")
    intergroup.add_argument('-p', '--min-pct-of-sum', metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_SUM,
        help="the minimum percentage of reads w.r.t. the marker's total "
             "number of reads (default: %(default)s)")
    intergroup.add_argument('-c', '--min-correction', metavar="PCT",
        type=float, default=_DEF_MIN_CORRECTION,
        help="the minimum change in read count due to correction by e.g., "
             "bgcorrect (default: %(default)s)")
    intergroup.add_argument('-y', '--min-recovery', metavar="PCT",
        type=float, default=_DEF_MIN_RECOVERY,
        help="the minimum number of reads that was recovered thanks to "
             "noise correction (by e.g., bgcorrect), as a percentage of the "
             "total number of reads after correction (default: %(default)s)")
    filtergroup = parser.add_argument_group("filtering options",
        "sequences that match the -C or -Y option (or both) and all of the "
        "other settings are retained, all others are filtered")
    filtergroup.add_argument('-a', '--filter-action', metavar="ACTION",
        choices=("off", "combine", "delete"), default="off",
        help="filtering mode: 'off', disable filtering; 'combine', replace "
             "filtered sequences by a single line with aggregate values per "
             "marker; 'delete', remove filtered sequences without leaving a "
             "trace (default: %(default)s)")
    filtergroup.add_argument('-A', '--filter-absolute', action="store_true",
        help="if specified, apply filters to absolute read counts (i.e., with "
             "the sign removed), which may keep over-corrected sequences that "
             "would otherwise be filtered out")
    filtergroup.add_argument('-N', '--min-reads-filt', metavar="N", type=float,
        default=_DEF_MIN_READS_FILT,
        help="the minimum number of reads (default: %(default)s)")
    filtergroup.add_argument('-B', '--min-per-strand-filt', metavar="N",
        type=float, default=_DEF_MIN_PER_STRAND_FILT,
        help="the minimum number of reads in both orientations (default: "
             "%(default)s)")
    filtergroup.add_argument('-M', '--min-pct-of-max-filt', metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_MAX_FILT,
        help="the minimum percentage of reads w.r.t. the highest allele of "
             "the marker (default: %(default)s)")
    filtergroup.add_argument('-P', '--min-pct-of-sum-filt', metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_SUM_FILT,
        help="the minimum percentage of reads w.r.t. the marker's total "
             "number of reads (default: %(default)s)")
    filtergroup.add_argument('-C', '--min-correction-filt', metavar="PCT",
        type=float, default=_DEF_MIN_CORRECTION_FILT,
        help="the minimum change in read count due to correction by e.g., "
             "bgcorrect (default: %(default)s)")
    filtergroup.add_argument('-Y', '--min-recovery-filt', metavar="PCT",
        type=float, default=_DEF_MIN_RECOVERY_FILT,
        help="the minimum number of reads that was recovered thanks to "
             "noise correction (by e.g., bgcorrect), as a percentage of the "
             "total number of reads after correction (default: %(default)s)")
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
        compute_stats(infile, outfile,
                      args.min_reads, args.min_per_strand, args.min_pct_of_max,
                      args.min_pct_of_sum, args.min_correction,
                      args.min_recovery, args.filter_action,
                      args.filter_absolute, args.min_reads_filt,
                      args.min_per_strand_filt, args.min_pct_of_max_filt,
                      args.min_pct_of_sum_filt, args.min_correction_filt,
                      args.min_recovery_filt)
        if infile != sys.stdin:
            infile.close()
#run
