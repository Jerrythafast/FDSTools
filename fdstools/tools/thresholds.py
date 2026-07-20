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
Calculate the pct_of_highest_allele of the lowest allele and the pct_of_allelic of the highest noise
in different categories (stutter, novel, same size novel and other) from infile. 

Return a table with highest percentages for each category in each samplefile given.
For homozygotes, the lowest allele will be 100%. 

If allelelist is not given, the flags column is used to determine alleles.

Default sequence format is raw.
"""
import logging

from ..lib.cli import add_bulk_input_output_args, add_sequence_format_args, iterate_input_samples
from ..lib.io import std_or_open, parse_allelelist, parse_flags, read_sample_data_file
from ..lib.seq import SEQ_SPECIAL_VALUES

__version__ = "1.0.0"

LOG = logging.getLogger(__name__)

def process_sample_data(sampledata, sample_alleles):
    marker_allele_sizes = {
        marker: {len(allele) for allele in alleles}
        for marker, alleles in sample_alleles.items()}
    marker_max = {
        marker: max((
            reads
            for (marker2, sequence), (reads, extra_col) in sampledata.items()
            if marker2 == marker),
            default=0)
        for marker in sample_alleles}
    processed_data = {
        marker: {category: [] for category in ("allele", "stutter", "same_size", "novel", "other")}
        for marker in sample_alleles}
    for (marker, sequence), (reads, extra_col) in sampledata:
        if marker not in sample_alleles:
            continue
        mp_max = reads / marker_max[marker]
        flags = parse_flags(extra_col["flags"])
        if sequence in sample_alleles[marker]:
            category = "allele"
        elif any(flag.startswith("STUTTER") for flag in flags):
            category = "stutter"
            # Now only works with stuttermark!
        elif "novel" in flags:
            if len(sequence) in marker_allele_sizes[marker]:
                category = "same_size"
            else:
                category = "novel"
        else:
            category = "other"
        processed_data[marker][category].append(mp_max)
    return processed_data
#process_sample_data

def determine_lowest_allele_highest_noise(sampledata, sample_alleles, sample, results):
    markers = set(key[0] for key in sampledata.keys()) # Get unique markers

    for marker in markers:
        if marker not in sample_alleles:
            LOG.warning("WARNING: Skipping marker %s in sample %s", marker, sample)
            continue

        allele_data = sampledata[marker]["allele"]

        lowest_allele = min(allele_data, default=0) # Percentage of highest allele. For homozygotes, this number will be 100. 
       
        scale_to_pct_allelic = 100. / sum(allele_data)      # Factor to convert mp_max to mp_allelic
        highest_stutter = max(sampledata[marker]["stutter"], 0.) * scale_to_pct_allelic
        highest_novel = max(sampledata[marker]["novel"], 0.) * scale_to_pct_allelic
        highest_samesize = max(sampledata[marker]["same_size"], 0.) * scale_to_pct_allelic
        highest_other = max(sampledata[marker]["other"], 0.) * scale_to_pct_allelic

        if marker not in results:
            results[marker] = list()

        results[marker].append((sample, lowest_allele, highest_stutter, highest_novel, highest_samesize, highest_other))

def analyse_thresholds(samples_in, outfile, allelefile, library):

    allelelist = {} if allelefile is None else parse_allelelist(
        allelefile, convert="raw", library=library)

    # Read data file of each sample and calculate for each sample the heterozygote balance per marker
    results_noise_alleles = {}
    for tag, infile in samples_in:
        if allelefile is not None and tag not in allelelist:
            continue

        alleles_out = allelelist.setdefault(tag, {}) if allelefile is None else None
        sampledata = read_sample_data_file(infile, seqformat="raw", library=library,
            drop_special_seq=True, after_correction=True, combine_strands=True,
            extra_columns={"flags": True}, alleles_out=alleles_out)
        
        processed_sampledata = process_sample_data(sampledata, allelelist)

        # Determine the pct_max for the lowest allele for heterozygotes for each marker of this sample
        determine_lowest_allele_highest_noise(processed_sampledata, allelelist[tag], tag, results_noise_alleles)

    # Wrting lowest alleles per marker to output (either file or terminal)
    outfile.write("\t".join(("sample", "marker", "lowest_allele", "highest_stutter", "highest_novel", "highest_same_size_novel", "highest_other")))
    outfile.write("\n")

    for marker in results_noise_alleles:
        for (sample, lowest_allele, highest_stutter, highest_novel, highest_samesize, highest_other) in results_noise_alleles[marker]:
            outfile.write(f"{sample}\t{marker}\t{lowest_allele}\t{highest_stutter}\t"
                          f"{highest_novel}\t{highest_samesize}\t{highest_other}\n")
#analyse_thresholds


def add_arguments(parser):
    add_bulk_input_output_args(parser)
    parser.add_argument("-a", "--allelelist", metavar="ALLELEFILE",
        help="file containing a list of the true alleles of each sample (e.g., obtained from "
             "allelefinder); if not given, the 'allele' flag from the 'flags' column is used")
    add_sequence_format_args(parser, default_format="raw", force=True)
#add_arguments


def run(args):
    samples_in = iterate_input_samples(args.infiles, args.tag_expr, args.tag_format)
    with std_or_open(args.outfile, "wt") as outfile:
        analyse_thresholds(samples_in, outfile, args.allelelist, args.library)
#run