#!/usr/bin/env python3

#
# Copyright (C) 2022 Jerry Hoogenboom
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

from errno import EPIPE

from ..lib.cli import add_input_output_args, get_input_output_files, add_sequence_format_args, \
                      pos_int_arg
from ..lib.io import get_column_ids, parse_flags
from ..lib.library import get_max_expected_alleles
from ..lib.seq import SEQ_SPECIAL_VALUES, ensure_sequence_format

__version__ = "1.3.1"


# Default values for parameters are specified below.

# Default minimum number of reads to mark as allele.
# This value can be overridden by the -n command line option.
_DEF_MIN_READS = 30

# Default minimum number of reads per strand to mark as allele.
# This value can be overridden by the -b command line option.
_DEF_MIN_PER_STRAND = 0

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

# Default minimum total number of allelic reads.
# This value can be overridden by the -E command line option.
_DEF_MIN_ALLELE_READS = 0

# Default maximum amount of noise to allow, as a percentage of the
# total number of allelic reads of this marker.  If any noise
# (i.e., non-allelic sequences) above this threshold are detected, the
# sample is considered 'noisy' for this marker.
# This value can be overridden by the -D command line option.
_DEF_MAX_NONALLELE_PCT = 100.

# Default minimum number of reads for filtering.
# This value can be overridden by the -N command line option.
_DEF_MIN_READS_FILT = 1

# Default minimum number of reads per strand for filtering.
# This value can be overridden by the -B command line option.
_DEF_MIN_PER_STRAND_FILT = 0

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

COLUMN_ORDER = (
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
)


def max_in_sequence(data, colid_max, colid_sequence, value_if_empty=0.):
    """
    Return the maximum data[][colid_max], where data[][colid_sequence]
    is not a special sequence value, or value_if_empty if nothing found.
    """
    try:
        return max(row[colid_max] for row in data if row[colid_sequence] not in SEQ_SPECIAL_VALUES)
    except ValueError:
        return value_if_empty
#max_in_sequence


def get_sums_and_maxes(markerdata, ci, columns):
    """
    Return a dict with the sum and max of each of the given columns.
    """
    return {
        column: {
            "sum": sum(row[ci[column]] for row in markerdata),
            "max": max_in_sequence(markerdata, ci[column], ci["sequence"])}
        for column in columns if column in ci}
#get_sum_and_max


def pct_or_default(numerator, denominator, default=0):
    """Return 100 * numerator / denominator, or 100*default on error."""
    return 100 * (numerator / denominator if denominator else default)
#pct_or_default


def append_if_missing(l, item):
    """Append item to list l if not already added."""
    if item not in l:
        l.append(item)
#append_if_missing


def set_or_append(l, i, item):
    """Store item at index i of list l if empty; append in case of IndexError."""
    try:
        if not l[i]:
            l[i] = item
        else:
            l[i] = float(l[i])
    except IndexError:
        l.append(item)
#set_or_append


def add_pct_column_names(column_names, strand, state):
    """
    If strand + state is not in column_names, do nothing.
    Else, add all applicable related columns to column_names.
    """
    if strand + state not in column_names:
        return
    append_if_missing(column_names, strand + state + "_mp_sum")
    append_if_missing(column_names, strand + state + "_mp_max")
    if state == "_corrected":
        append_if_missing(column_names, strand + "_correction_pct")
    elif state == "_noise":
        append_if_missing(column_names, strand + "_removed_pct")
    elif state == "_add":
        append_if_missing(column_names, strand + "_added_pct")
    if strand == "total":
        fwdcol = "forward" + state
        if fwdcol in column_names:
            append_if_missing(column_names, fwdcol + "_pct")
    if state == "_corrected":
        if strand + "_add" in column_names:
            append_if_missing(column_names, strand + "_recovery")
#add_pct_column_names


def add_pct_columns(row, ci, strand, state, aggr):
    """
    If strand + state is not in ci, do nothing.  Else, calculate all
    related percentages for row (adding if not already present).
    """
    column = strand + state
    if column not in ci:
        return
    set_or_append(row, ci[column + "_mp_sum"],
        pct_or_default(row[ci[column]], aggr[column]["sum"]))
    set_or_append(row, ci[column + "_mp_max"],
        pct_or_default(row[ci[column]], aggr[column]["max"]))
    if state == "_corrected":
        set_or_append(row, ci[strand + "_correction_pct"], pct_or_default(
            row[ci[column]] - row[ci[strand]], row[ci[strand]],
            (row[ci[column]] > 0) * 2 - 1 if row[ci[column]] else 0))
    elif state == "_noise":
        set_or_append(row, ci[strand + "_removed_pct"],
            pct_or_default(row[ci[column]], row[ci[strand]]))
    elif state == "_add":
        set_or_append(row, ci[strand + "_added_pct"],
            pct_or_default(row[ci[column]], row[ci[strand]]))
    if strand == "total":
        fwdcol = "forward" + state
        if fwdcol in ci:
            set_or_append(row, ci[fwdcol + "_pct"],
                pct_or_default(row[ci[fwdcol]], row[ci[column]], int(row[ci[fwdcol]] > 0)))
    if state == "_corrected":
        addcol = strand + "_add"
        if addcol in ci:
            set_or_append(row, ci[strand + "_recovery"],
                pct_or_default(row[ci[addcol]], row[ci[column]]))
#add_pct_columns


def compute_stats(infile, outfile, min_reads, min_per_strand, min_pct_of_max, min_pct_of_sum,
                  min_correction, min_recovery, min_allele_reads, max_nonallele_pct, filter_action,
                  filter_absolute, max_alleles, min_reads_filt, min_per_strand_filt,
                  min_pct_of_max_filt, min_pct_of_sum_filt, min_correction_filt,
                  min_recovery_filt, library, sequence_format, uncall_alleles):
    # Check presence of required columns.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return  # Empty file.
    get_column_ids(column_names, "marker", "sequence", "forward", "reverse", "total")

    # Add columns for which we have the required data.
    append_if_missing(column_names, "flags")
    for state in ("_corrected", "", "_noise", "_add"):
        for strand in ("total", "forward", "reverse"):
            add_pct_column_names(column_names, strand, state)

    # Build a column number lookup dictionary.
    ci = {name: i for i, name in enumerate(column_names)}
    total_column = "total_corrected" if "total_corrected" in ci else "total"

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
            row[ci["flags"]] = parse_flags(row[ci["flags"]])
        if marker not in data:
            data[marker] = []
        data[marker].append(row)

    # Compute statistics.
    if filter_action != "off":
        filtered = {marker: [] for marker in data}
    for marker in data:
        # Calculate aggregate statistics (sums and maxes).
        aggr = get_sums_and_maxes(data[marker], ci, (strand + state
            for state in ("_corrected", "", "_noise", "_add")
            for strand in ("total", "forward", "reverse")))

        allele_reads = 0
        max_noise_reads = 0
        alleles = []
        for row in data[marker]:
            # Add various percentage columns to this row.
            for state in ("_corrected", "", "_noise", "_add"):
                for strand in ("total", "forward", "reverse"):
                    add_pct_columns(row, ci, strand, state, aggr)

            # The 'No data' lines are fine like this.
            if row[ci["sequence"]] == "No data":
                row[ci["flags"]] = ",".join(row[ci["flags"]])
                continue

            # Get the values we will filter on.
            is_aggregate = row[ci["sequence"]] == "Other sequences"
            total_reads = row[ci[total_column]]
            pct_of_sum = row[ci["total_mp_sum"]] if "total_corrected_mp_sum" \
                not in ci else row[ci["total_corrected_mp_sum"]]
            pct_of_max = row[ci["total_mp_max"]] if "total_corrected_mp_max" \
                not in ci else row[ci["total_corrected_mp_max"]]
            correction = row[ci["total_correction_pct"]] if "total_correction_pct" in ci else 0
            recovery = row[ci["total_recovery"]] if "total_recovery" in ci else 0
            strands = (
                row[ci["forward" if "forward_corrected" not in ci else "forward_corrected"]],
                row[ci["reverse" if "reverse_corrected" not in ci else "reverse_corrected"]])
            fn = abs if filter_absolute else lambda x: x

            # Check if this sequence should be filtered out.
            # Always filter/combine existing 'Other sequences'.
            if filter_action != "off" and (
                    is_aggregate or (
                    fn(total_reads) < min_reads_filt or
                    fn(pct_of_max) < min_pct_of_max_filt or
                    fn(pct_of_sum) < min_pct_of_sum_filt or
                    (correction < min_correction_filt and
                    recovery < min_recovery_filt) or
                    min(map(fn, strands)) < min_per_strand_filt)):
                filtered[marker].append(row)
                continue

            if is_aggregate:
                continue

            pass_thresholds = (
                    total_reads >= min_reads and
                    pct_of_max >= min_pct_of_max and
                    pct_of_sum >= min_pct_of_sum and
                    (correction >= min_correction or
                    recovery >= min_recovery) and
                    min(strands) >= min_per_strand)

            # Check if this sequence is an allele.
            if (pass_thresholds and ("allele" in row[ci["flags"]] or uncall_alleles != "only")) or (
                    not pass_thresholds and "allele" in row[ci["flags"]] and not uncall_alleles):
                append_if_missing(row[ci["flags"]], "allele")
                allele_reads += total_reads
                alleles.append(row)
            else:
                if uncall_alleles and "allele" in row[ci["flags"]]:
                    # Un-call allele.
                    row[ci["flags"]].remove("allele")

                # Check if this sequence is the highest noise.
                if total_reads > max_noise_reads:
                    max_noise_reads = total_reads

        if (not max_alleles == 0 and len(alleles) >
                get_max_expected_alleles(max_alleles, marker, library)) or (allele_reads and
                (allele_reads < min_allele_reads or
                100 * max_noise_reads / allele_reads > max_nonallele_pct)):
            # Remove allele flags if there are more alleles than expected,
            # total allele coverage is too low, or noise is too high.
            for allele in alleles:
                allele[ci["flags"]].remove("allele")

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
                for column in (strand + state + aggr
                        for state in ("_corrected", "", "_noise", "_add")
                        for strand in ("total", "forward", "reverse")
                        for aggr in ("", "_mp_sum", "_mp_max")):
                    if column in ci:
                        combined[ci[column]] += row[ci[column]]
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
                row[ci["flags"]] = ",".join(row[ci["flags"]])
                if sequence_format is not None:
                    row[ci["sequence"]] = ensure_sequence_format(
                        row[ci["sequence"]], sequence_format, library=library, marker=marker)
                outfile.write("\t".join(map(str,
                    (row[new_order[i]] for i in range(len(row))))) + "\n")

        # Add combined row for this marker.
        if filter_action == "combine" and have_combined:
            if "total_corrected" in ci:
                combined[ci["total_correction_pct"]] = pct_or_default(
                    combined[ci["total_corrected"]] - combined[ci["total"]],
                    combined[ci["total"]],
                    (combined[ci["total_corrected"]] > 0) * 2 - 1
                        if combined[ci["total_corrected"]] else 0)
                if "forward_corrected" in ci:
                    combined[ci["forward_corrected_pct"]] = pct_or_default(
                        combined[ci["forward_corrected"]], combined[ci["total_corrected"]],
                        int(combined[ci["forward_corrected"]] > 0))
                if "total_add" in ci:
                    combined[ci["total_recovery"]] = pct_or_default(
                        combined[ci["total_add"]], combined[ci["total_corrected"]])
            if "forward_corrected" in ci:
                combined[ci["forward_correction_pct"]] = pct_or_default(
                    combined[ci["forward_corrected"]] - combined[ci["forward"]],
                    combined[ci["forward"]],
                    (combined[ci["forward_corrected"]] > 0) * 2 - 1
                        if combined[ci["forward_corrected"]] else 0)
                if "forward_add" in ci:
                    combined[ci["forward_recovery"]] = pct_or_default(
                        combined[ci["forward_add"]], combined[ci["forward_corrected"]])
            if "reverse_corrected" in ci:
                combined[ci["reverse_correction_pct"]] = pct_or_default(
                    combined[ci["reverse_corrected"]] - combined[ci["reverse"]],
                    combined[ci["reverse"]],
                    (combined[ci["reverse_corrected"]] > 0) * 2 - 1
                        if combined[ci["reverse_corrected"]] else 0)
                if "reverse_add" in ci:
                    combined[ci["reverse_recovery"]] = pct_or_default(
                        combined[ci["reverse_add"]], combined[ci["reverse_corrected"]])
            combined[ci["forward_pct"]] = pct_or_default(
                combined[ci["forward"]], combined[ci["total"]], int(combined[ci["forward"]] > 0))
            if "total_noise" in ci:
                combined[ci["total_removed_pct"]] = pct_or_default(
                    combined[ci["total_noise"]], combined[ci["total"]])
                if "forward_noise" in ci:
                    combined[ci["forward_noise_pct"]] = pct_or_default(
                        combined[ci["forward_noise"]], combined[ci["total_noise"]],
                        int(combined[ci["forward_noise"]] > 0))
            if "forward_noise" in ci:
                combined[ci["forward_removed_pct"]] = pct_or_default(
                    combined[ci["forward_noise"]], combined[ci["forward"]])
            if "reverse_noise" in ci:
                combined[ci["reverse_removed_pct"]] = pct_or_default(
                    combined[ci["reverse_noise"]], combined[ci["reverse"]])
            if "total_add" in ci:
                combined[ci["total_added_pct"]] = pct_or_default(
                    combined[ci["total_add"]], combined[ci["total"]])
                if "forward_add" in ci:
                    combined[ci["forward_add_pct"]] = pct_or_default(
                        combined[ci["forward_add"]], combined[ci["total_add"]],
                        int(combined[ci["forward_add"]] > 0))
            if "forward_add" in ci:
                combined[ci["forward_added_pct"]] = pct_or_default(
                    combined[ci["forward_add"]], combined[ci["forward"]])
            if "reverse_add" in ci:
                combined[ci["reverse_added_pct"]] = pct_or_default(
                    combined[ci["reverse_add"]], combined[ci["reverse"]])

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
    add_input_output_args(parser, single_in=True, batch_support=True, report_out=False)
    intergroup = parser.add_argument_group("interpretation options",
        "sequences that match the -c or -y option (or both) and all of the "
        "other settings are marked as 'allele'")
    intergroup.add_argument("-U", "--uncall-alleles", nargs="?", metavar="only",
        choices=("only",), const=True, default=False,
        help="if specified and the input contains sequences with the 'allele' flag, the flag will "
             "be removed for sequences not meeting the requirements; with the optional keyword "
             "'only', no 'allele' flags will be added to any sequences that do meet the criteria")
    intergroup.add_argument("-n", "--min-reads", metavar="N", type=float, default=_DEF_MIN_READS,
        help="the minimum number of reads (default: %(default)s)")
    intergroup.add_argument("-b", "--min-per-strand", metavar="N", type=float,
        default=_DEF_MIN_PER_STRAND,
        help="the minimum number of reads in both orientations (default: %(default)s)")
    intergroup.add_argument("-m", "--min-pct-of-max", metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_MAX,
        help="the minimum percentage of reads w.r.t. the highest allele of "
             "the marker (default: %(default)s)")
    intergroup.add_argument("-p", "--min-pct-of-sum", metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_SUM,
        help="the minimum percentage of reads w.r.t. the marker's total "
             "number of reads (default: %(default)s)")
    intergroup.add_argument("-c", "--min-correction", metavar="PCT",
        type=float, default=_DEF_MIN_CORRECTION,
        help="the minimum percentage change in read count due to correction by e.g., "
             "bgcorrect (total_correction column; default: %(default)s)")
    intergroup.add_argument("-y", "--min-recovery", metavar="PCT",
        type=float, default=_DEF_MIN_RECOVERY,
        help="the minimum number of reads that was recovered thanks to "
             "noise correction (by e.g., bgcorrect), as a percentage of the total number "
             "of reads after correction (total_recovery column; default: %(default)s)")
    intergroup.add_argument('-E', '--min-allele-reads', metavar="N",
        type=float, default=_DEF_MIN_ALLELE_READS,
        help="force a minimum total number of reads for all alleles on a marker; don't call "
             "any alleles otherwise (default: %(default)s)")
    intergroup.add_argument('-D', '--max-nonallele-pct', metavar="PCT",
        type=float, default=_DEF_MAX_NONALLELE_PCT,
        help="drop all allele markings if the highest non-allelic sequence is at least this "
             "percentage of the total number of reads for all alleles on that marker "
             "(default: %(default)s)")
    intergroup.add_argument("-G", "--max-alleles", metavar="N", type=pos_int_arg,
        nargs="?", default=0,
        help="if specified, do not mark any alleles on a marker if more than N alleles meet the "
             "criteria; without N, the amounts given in the library file are used, which have a "
             "default value of 1 for markers on the mitochondrial genome and Y chromosome, or 2 "
             "otherwise (Note: don't forget to provide -l/--library!)")
    filtergroup = parser.add_argument_group("filtering options",
        "sequences that match the -C or -Y option (or both) and all of the "
        "other settings are retained, all others are filtered")
    filtergroup.add_argument("-a", "--filter-action", metavar="ACTION",
        choices=("off", "combine", "delete"), default="off",
        help="filtering mode: 'off', disable filtering; 'combine', replace "
             "filtered sequences by a single line with aggregate values per "
             "marker; 'delete', remove filtered sequences without leaving a "
             "trace (default: %(default)s)")
    filtergroup.add_argument("-A", "--filter-absolute", action="store_true",
        help="if specified, apply filters to absolute read counts (i.e., with "
             "the sign removed), which may keep over-corrected sequences that "
             "would otherwise be filtered out")
    filtergroup.add_argument("-N", "--min-reads-filt", metavar="N", type=float,
        default=_DEF_MIN_READS_FILT,
        help="the minimum number of reads (default: %(default)s)")
    filtergroup.add_argument("-B", "--min-per-strand-filt", metavar="N",
        type=float, default=_DEF_MIN_PER_STRAND_FILT,
        help="the minimum number of reads in both orientations (default: %(default)s)")
    filtergroup.add_argument("-M", "--min-pct-of-max-filt", metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_MAX_FILT,
        help="the minimum percentage of reads w.r.t. the highest allele of "
             "the marker (default: %(default)s)")
    filtergroup.add_argument("-P", "--min-pct-of-sum-filt", metavar="PCT",
        type=float, default=_DEF_MIN_PCT_OF_SUM_FILT,
        help="the minimum percentage of reads w.r.t. the marker's total "
             "number of reads (default: %(default)s)")
    filtergroup.add_argument("-C", "--min-correction-filt", metavar="PCT",
        type=float, default=_DEF_MIN_CORRECTION_FILT,
        help="the minimum percentage change in read count due to correction by e.g., "
             "bgcorrect (total_correction column; default: %(default)s)")
    filtergroup.add_argument("-Y", "--min-recovery-filt", metavar="PCT",
        type=float, default=_DEF_MIN_RECOVERY_FILT,
        help="the minimum number of reads that was recovered thanks to "
             "noise correction (by e.g., bgcorrect), as a percentage of the "
             "total number of reads after correction (total_recovery column; default: %(default)s)")
    add_sequence_format_args(parser)
#add_arguments


def run(args):
    gen = get_input_output_files(args, single_in=True, batch_support=True)
    if not gen:
        raise ValueError("please specify an input file, or pipe in the output of another program")

    for tag, infiles, outfile in gen:
        # TODO: Aggregate data from all infiles of each sample.
        if len(infiles) > 1:
            raise ValueError("multiple input files for sample '%s' specified " % tag)
        try:
            infile = sys.stdin if infiles[0] == "-" else open(infiles[0], "rt", encoding="UTF-8")
            compute_stats(infile, outfile, args.min_reads, args.min_per_strand,
                          args.min_pct_of_max, args.min_pct_of_sum, args.min_correction,
                          args.min_recovery, args.min_allele_reads, args.max_nonallele_pct,
                          args.filter_action, args.filter_absolute, args.max_alleles,
                          args.min_reads_filt, args.min_per_strand_filt, args.min_pct_of_max_filt,
                          args.min_pct_of_sum_filt, args.min_correction_filt,
                          args.min_recovery_filt, args.library, args.sequence_format,
                          args.uncall_alleles)
            if infile != sys.stdin:
                infile.close()
        except IOError as e:
            if e.errno == EPIPE:
                continue
            raise
#run
