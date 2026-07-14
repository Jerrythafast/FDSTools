#!/usr/bin/env python3

#
# Copyright (C) 2026 Netherlands Forensic Institute
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
This tool calculates the pct_max of the lowest allele and the pct_allelic of the highest noise...

This script:
- reads the allelelist and saves the alleles for each sample
- reads the samplestats file (given version) for each sample and saves the total_corrected_pct_max of each allele
- determines for heterozygotes the pct_max of the lowest allele
- determines for both heterozygotes as homozygotes the highest noise pct_max
- prints the results in the terminal. To print to file, write '> filename.txt' after the command

CHANGES IN 1.1:
- checks if all alleles in allelelist are in samplestats file, otherwise uses 0 as value.
- if there is no noise, the marker would be skipped in previous version. Now fixed.
- gives a warning if a marker is not in the allelelist, because it will be skipped.
"""
import argparse, sys, re, os

from glob import glob

from ..lib.cli import add_sequence_format_args, pos_int_arg
from ..lib.io import std_or_open, parse_allelelist
from ..lib.seq import ensure_sequence_format
from ..lib.tssv import filter_sequences
from ..lib.util import batched, thread_safe_process_pool


__version__ = "1.0.0"


# Default minimum number of reads to consider.
# This value can be overridden by the -a command line option.
_DEF_MINIMUM = 2


# Sequences to be filtered out.
FILTERS = (
    (None, "Other sequences"),
    (None, "No data")
)

def read_sample_file(instream):
    """Reads data from the given input stream"""
    # Read first line, which contains the column names.
    # Remove the line ending and split it into a list on the tab chars.
    column_names = instream.readline().rstrip("\r\n").split("\t")

    # Find the column indices of the columns we need.
    ci = {colname: column_names.index(colname) for colname in (
        "marker", "sequence", "total_corrected_mp_max", "flags")}

    # Now read the data line by line.
    data = {}
    for line in instream:
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < len(column_names): # Otherwise it breaks on the stats on the bottom of the table
            continue
        marker = fields[ci["marker"]]
        sequence = fields[ci["sequence"]]
        if (marker, sequence) in FILTERS or (None, sequence) in FILTERS:
            continue

        total_corrected_mp_max = float(fields[ci["total_corrected_mp_max"]])

        flags = fields[ci["flags"]].split(',')
        if any(flag.startswith("STUTTER") for flag in flags):
            # Stutter noise
            category = "stutter"
        elif "novel" in flags:
            category = "novel" # will split out same_size later
        else:
            category = "other"

        data[marker, sequence] = [total_corrected_mp_max, category]
    return data

def determine_lowest_allele_highest_noise(sampledata, sample, results):
    markers = set(key[0] for key in sampledata.keys()) # Get unique markers

    for marker in markers:
        if marker not in allelelist[sample]:
            print(f"WARNING: Skipping marker {marker} in sample {sample}", file=sys.stderr)
            continue

        allele_data = [
            sampledata.get((marker, sequence), [0])[0]  # Get value from sampledata, or 0 if not detected
            for sequence in allelelist[sample][marker]]
        noise_data = {
            sequence: value
            for (m, sequence), value in sampledata.items()
            if m == marker and sequence not in allelelist[sample][marker]}

        # Change category of same size noise.
        allele_ce = set(allele.split("_")[0] for allele in allelelist[sample][marker])

        if len(allele_data) > 0:
            lowest_allele = min(allele_data) # For homozygotes, this number will be 100
        else:
            lowest_allele = 0

        if len(noise_data) > 0:
            for sequence, noise in noise_data.items():
                if noise[1] == "novel" and sequence.split("_")[0] in allele_ce:
                    noise[1] = "same_size"

            sorted_noise_data = dict(sorted(noise_data.items(), key=lambda item: item[1][0], reverse=True)) # sort noise_data dict based on pct_max value (first value in tuple)
            highest_stutter = 100. * next((mp_max for mp_max, category in sorted_noise_data.values() if category == "stutter"), 0.) / sum(allele_data)
            highest_novel = 100. * next((mp_max for mp_max, category in sorted_noise_data.values() if category == "novel"), 0.) / sum(allele_data)
            highest_samesize = 100. * next((mp_max for mp_max, category in sorted_noise_data.values() if category == "same_size"), 0.) / sum(allele_data)
            highest_other = 100. * next((mp_max for mp_max, category in sorted_noise_data.values() if category == "other"), 0.) / sum(allele_data)
        else:
            highest_stutter = 0
            highest_novel = 0
            highest_samesize = 0
            highest_other = 0


        if marker not in results:
            results[marker] = list()

        results[marker].append((sample, lowest_allele, highest_stutter, highest_novel, highest_samesize, highest_other))




def analyse_thresholds(samples_in, outfile, allelefile, library):

    allelelist = parse_allelelist(allelefile, convert="allelename", library=library)

    # Read data file of each sample and calculate for each sample the heterozygote balance per marker
    results_noise_alleles = {}
    for tag, infile in samples_in:
        sampledata = read_sample_file(infile)

        # Determine the pct_max for the lowest allele for heterozygotes for each marker of this sample
        determine_lowest_allele_highest_noise(sampledata, tag, results_noise_alleles)



    # Wrting lowest alleles per marker to output (either file or terminal)
    outfile.write("Pct_max of lowest alleles and pct_of_allelic of highest noise variants in the different OmniSTR markers (if lowest = 100, then marker is homozygous):\n")
    outfile.write("\t".join(["samplename", "marker", "lowest allele (pct_max)", "highest stutter (pct_allelic)", "highest novel (pct_allelic)", "highest samesize novel (pct_allelic)", "highest other noise (pct_allelic)"]))
    outfile.write("\n")

    for marker in results_noise_alleles:
        for (sample, lowest_allele, highest_stutter, highest_novel, highest_samesize, highest_other) in results_noise_alleles[marker]:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, marker, lowest_allele, highest_stutter, highest_novel, highest_samesize, highest_other))
#analyse_thresholds


def add_arguments(parser):
    add_sequence_format_args(parser, default_format="allelename", force=True)
    add_bulk_input_output_args(parser)
    parser.add_argument("-a", "--allelelist", metavar="ALLELEFILE",
        help="file containing a list of the true alleles of each sample (e.g., obtained from "
             "allelefinder); if not given, the 'allele' flag from the 'flags' column is used")
#add_arguments


def run(args):

    SAMPLESTAT_LOCATION = (args.d)
    VERSION = (args.v) # Should be in format 'v1' or 'v1a' - no underscores!
    TAG_PATTERN = re.compile(r"((?:NS|Mi).+)_(E\d{4})") # It can handle both NextSeq and MiSeq runs
    SAMPLESTATS_NAME_BEFORE_TAG = f"samplestats_{VERSION}_"
    SAMPLESTATS_NAME_AFTER_TAG = "-done.txt"

    samples_in = iterate_input_samples(args.infiles, args.tag_expr, args.tag_format)
    with std_or_open(args.outfile, "wt") as outfile:
        analyse_thresholds(samples_in, outfile, args.allelelist, args.library)

#run


