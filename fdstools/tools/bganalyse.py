#!/usr/bin/env python3

#
# Copyright (C) 2020 Jerry Hoogenboom
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
Analyse the amount of noise in reference samples.

Use this tool after correcting the reference samples with BGCorrect to
analyse the amount of remaining noise after correction.  This way,
potentially contaminated or otherwise 'dirty' reference samples can be
detected.  The highest amount of remaining noise can be interpreted as a
lower bound to the reliable detection of a minor contributor's alleles
in mixed DNA samples.

In the default mode ('full'), the lowest, highest, and total number of
backgroud/noise reads as well as the respective percentages w.r.t. the
number of allelic reads of each marker in each sample is printed.  This
data can be visualised using 'fdstools vis bganalyse'.

In the alternative 'percentiles' mode, the highest and the total number
of background reads as a percentage of the number of allelic reads for
each marker is given at selected percentiles of the samples.  I.e., it
gives the highest and total remaining noise considering only the
cleanest x% of samples, for different values of x.
"""
import math

from errno import EPIPE

from ..lib.cli import add_sequence_format_args, add_input_output_args, get_input_output_files, \
                      add_allele_detection_args, comma_separated_arg
from ..lib.io import get_sample_data, parse_allelelist

__version__ = "1.1.0"


# Default values for parameters are specified below.

# Default percentiles.
# This value can be overridden by the -p command line option.
_DEF_PERCENTILES = "100,99,95,90"


def get_total_recovery(values):
    if "total_recovery" in values:
        return float(values["total_recovery"])
    elif ("total_add" in values and "total_corrected" in values and
            float(values["total_corrected"])):
        return 100 * float(values["total_add"]) / float(values["total_corrected"])
    return 0
#get_total_recovery


def process_sample(sample_data, sample_alleles, tag, library):
    # Get read counts of the true alleles.
    allele_reads = {}
    for marker in sample_alleles:
        allele_reads[marker] = 0
        for allele in sample_alleles[marker]:
            if (marker, allele) not in sample_data:
                raise ValueError(
                    "Missing allele %s of marker %s in sample %s!" % (allele, marker, tag))
            this_allele_reads = sample_data[marker, allele][0]
            if not this_allele_reads:
                raise ValueError("Allele %s of marker %s has 0 reads!" % (allele, marker))
            allele_reads[marker] += this_allele_reads

    # Find the highest, the lowest, and the total noise for each marker,
    # and the highest total recovery.
    noise = {}
    marker_reads = {}
    total_reads = 0
    for marker, sequence in sample_data:
        reads = sample_data[marker, sequence][0]
        total_reads += reads
        if marker not in sample_alleles or not sample_alleles[marker]:
            # Sample does not participate in this marker.
            continue
        if marker not in noise:
            noise[marker] = [0, 0, 0, 0]  # Max, min, sum, max(recovery)
            marker_reads[marker] = reads
        else:
            marker_reads[marker] += reads
        if sequence in sample_alleles[marker]:
            # Note: Some noise has recovery in the order of 1E17.
            recovery = get_total_recovery(sample_data[marker, sequence][1])
            if recovery > noise[marker][3]:
                noise[marker][3] = recovery
        else:
            # Update lowest/highest noise as needed.
            if reads > noise[marker][0]:
                noise[marker][0] = reads
            if reads < noise[marker][1]:
                noise[marker][1] = reads
            noise[marker][2] += reads

    return {
        marker: {
            "genotype": sorted(sample_alleles[marker]),
            "allele_reads": allele_reads[marker],
            "noise": noise[marker],
            "marker_reads": marker_reads[marker],
            "total_reads": total_reads
        }
        for marker in noise}
#process_sample


def write_full_table(outfile, data):
    outfile.write("\t".join((
        "sample", "marker", "genotype", "allelic_reads",
        "highest_remaining_bg_reads", "lowest_remaining_bg_reads",
        "total_remaining_bg_reads", "highest_as_pct_of_allelic",
        "lowest_as_pct_of_allelic", "total_as_pct_of_allelic",
        "highest_recovery", "total_reads_marker", "total_reads_sample")) + "\n")
    for tag in data:
        for marker in data[tag]:
            d = data[tag][marker]
            outfile.write("\t".join([
                tag, marker, ",".join(d["genotype"])] + [
                "%.5g" % x if abs(x) >= 10000 else
                "%.4g" % x if abs(x) >= 1000 else
                "%.3g" % x if abs(x) > 0.0000000001 else "0" for x in (
                    d["allele_reads"],
                    d["noise"][0],
                    d["noise"][1],
                    d["noise"][2],
                    100 * d["noise"][0] / d["allele_reads"],
                    100 * d["noise"][1] / d["allele_reads"],
                    100 * d["noise"][2] / d["allele_reads"],
                    d["noise"][3],
                    d["marker_reads"],
                    d["total_reads"])]) + "\n")
#write_full_table


def write_percentiles_table(outfile, data, percentiles):
    per_marker = {}
    for tag in data:
        for marker in data[tag]:
            d = data[tag][marker]
            if marker not in per_marker:
                per_marker[marker] = ([], [], [], [])
            per_marker[marker][0].append(100 * d["noise"][0] / d["allele_reads"])
            per_marker[marker][1].append(100 * d["noise"][1] / d["allele_reads"])
            per_marker[marker][2].append(100 * d["noise"][2] / d["allele_reads"])
            per_marker[marker][3].append(d["noise"][3])

    outfile.write("\t".join((
        "marker", "percentile", "highest_as_pct_of_allelic",
        "lowest_as_pct_of_allelic", "total_as_pct_of_allelic",
        "highest_recovery")) + "\n")
    for marker in per_marker:
        per_marker[marker] = tuple(map(sorted, per_marker[marker]))
        n = len(per_marker[marker][0])
        for percentile in percentiles:
            i = math.ceil(percentile * 0.01 * n) - 1
            outfile.write("\t".join([
                marker, "%g" % percentile] + [
                "%.5g" % x if abs(x) >= 10000 else
                "%.4g" % x if abs(x) >= 1000 else
                "%.3g" % x if abs(x) > 0.0000000001 else "0" for x in (
                    per_marker[marker][0][i],
                    per_marker[marker][1][n - i - 1],
                    per_marker[marker][2][i],
                    per_marker[marker][3][i])]) + "\n")
#write_percentiles_table


def analyse_background(samples_in, outfile, allelefile, annotation_column, seqformat, library,
                       mode, percentiles):

    # Parse allele list.
    allelelist = {} if allelefile is None else parse_allelelist(
        allelefile, convert=seqformat, library=library)

    # TODO: Should we consider forward and reverse reads separately?
    data = {}
    get_sample_data(
        samples_in,
        lambda tag, sample_data: data.__setitem__(tag, process_sample(
            sample_data, allelelist[tag], tag, library)),
        allelelist, annotation_column, seqformat, library,
        drop_special_seq=True, after_correction=True, combine_strands=True,
        extra_columns={
            "total_recovery": True,
            "total_add": True,
            "total_corrected": True
        })

    if mode == "full":
        write_full_table(outfile, data)
    else:
        write_percentiles_table(outfile, data, percentiles)
#analyse_background


def add_arguments(parser):
    parser.add_argument("-m", "--mode", metavar="MODE", default="full",
        choices=("full", "percentiles"),
        help="controls what kind of information is printed; 'full' (the "
             "default) prints the lowest, highest, and total number of "
             "backgroud reads as well as the respective percentages w.r.t. "
             "the number of allelic reads of each marker in each "
             "sample; 'percentiles' prints the highest and the total number "
             "of background reads as a percentage of the number of allelic "
             "reads for each marker at given percentiles")
    parser.add_argument("-p", "--percentiles", metavar="PCT", default=_DEF_PERCENTILES,
        type=comma_separated_arg(tuple, float),
        help="comma-separated list of percentiles to report when -m/--mode is "
             "set to 'percentiles' (default: %(default)s)")
    add_input_output_args(parser)
    add_allele_detection_args(parser)
    add_sequence_format_args(parser, default_format="raw")
#add_arguments


def run(args):
    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output of another program")
    try:
        analyse_background(files[0], files[1], args.allelelist, args.annotation_column,
                           args.sequence_format, args.library, args.mode, args.percentiles)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#run
