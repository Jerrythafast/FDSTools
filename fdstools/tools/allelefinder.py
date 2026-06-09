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
Find true alleles in reference samples and detect possible
contaminations.

In each sample, the sequences with the highest read counts of each
marker are called alleles, with a user-defined maximum number of alleles
per marker.  The allele balance is kept within given bounds.  If the
highest non-allelic sequence exceeds a given limit, no alleles are
called for this marker.  If this happens for multiple markers in one
sample, no alleles are called for this sample at all.

If the input file contains a 'flags' column, any sequences with a flag
starting with 'STUTTER' will be ignored.  Therefore, it is highly
recommended to run Allelefinder on the output of Stuttermark.

The allele list obtained from allelefinder should always be checked
carefully before using it as the input of various other tools operating
on reference samples.  These tools rely heavily on the correctness of
this file to do their job.  One may use the allelefinder report
(-R/--report output argument) and the bganalyse tool to get a quick
overview of what might be wrong.
"""
import logging

from ..lib.cli import add_bulk_input_output_args, add_sequence_format_args, \
                      iterate_input_samples, pos_int_arg
from ..lib.io import parse_flags, read_sample_data_file, std_or_open
from ..lib.library import get_max_expected_alleles
from ..lib.seq import ensure_sequence_format

__version__ = "1.1.1"

LOG = logging.getLogger(__name__)


# Default values for parameters are specified below.

# Default minimum number of reads required for the highest allele.
# This value can be overridden by the -n command line option.
_DEF_MIN_READS = 50

# Default minimum number of reads required for the lowest allele.
# This value can be overridden by the -N command line option.
_DEF_MIN_READS_LOWEST = 15

# Default minimum number of reads required for an allele to be called,
# as a percentage of the number of reads of the highest allele.
# This value can be overridden by the -m command line option.
_DEF_MIN_ALLELE_PCT = 30.0

# Default maximum amount of noise to allow, as a percentage of the
# number of reads of the highest allele of each marker.  If any noise
# (i.e., non-allelic sequences) above this threshold are detected, the
# sample is considered 'noisy' for this marker.
# This value can be overridden by the -M command line option.
_DEF_MAX_NOISE_PCT = 10.0

# Default maximum number or fraction of noisy markers allowed per sample.
# This value can be overridden by the -x command line option.
_DEF_MAX_NOISY = 0.1


def report_or_log(reportfile, message):
    """Write message to reportfile, or to LOG.warning if reportfile is None."""
    if reportfile:
        reportfile.write(message)
    else:
        LOG.warning(message.strip().replace("\n", " "))
#report_or_log


def find_alleles(samples_in, outfile, reportfile, min_reads, min_reads_lowest, min_allele_pct,
                 max_noise_pct, max_alleles, max_noisy, seqformat, library):
    outfile.write("\t".join(["sample", "marker", "total", "allele"]) + "\n")
    for tag, infile in samples_in:
        data = read_sample_data_file(
            infile, drop_special_seq=True, after_correction=True, combine_strands=True,
            extra_columns={"flags": True})
        find_alleles_sample(
            data, outfile, reportfile, tag, min_reads, min_reads_lowest, min_allele_pct,
            max_noise_pct, max_alleles, max_noisy, seqformat, library)
#find_alleles


def find_alleles_sample(data, outfile, reportfile, tag, min_reads, min_reads_lowest,
                        min_allele_pct, max_noise_pct, max_alleles, max_noisy, seqformat, library):
    top_noise = {}
    top_allele = {}
    alleles = {}
    for marker, sequence in data:
        # Skip lines that have been marked as stutter artefacts.
        try:
            flags = parse_flags(data[marker, sequence].pop()["flags"])
            if any(flag.startswith("STUTTER") for flag in flags):
                continue
        except KeyError:
            # Flags column wasn't present.
            pass

        reads = data[marker, sequence][0]

        if marker not in alleles:
            alleles[marker] = {}
            top_allele[marker] = 0
            top_noise[marker] = ["-", 0]

        if reads > top_allele[marker]:
            # New highest allele!
            top_allele[marker] = reads
            for allele in tuple(alleles[marker]):
                if alleles[marker][allele] < top_allele[marker] * (min_allele_pct / 100):
                    if alleles[marker][allele] > top_noise[marker][1]:
                        top_noise[marker] = [allele, alleles[marker][allele]]
                    del alleles[marker][allele]
            alleles[marker][sequence] = reads
        elif reads >= top_allele[marker] * (min_allele_pct / 100):
            # New secundary allele!
            alleles[marker][sequence] = reads
        elif reads >= top_noise[marker][1]:
            # New highest noise!
            top_noise[marker] = [sequence, reads]

    # Find and eliminate noisy markers in this sample first.
    noisy_markers = 0
    for marker in alleles:
        if top_allele[marker] < min_reads:
            report_or_log(reportfile,
                f"Sample {tag} is not suitable for marker {marker}:\n"
                f"highest allele has only {top_allele[marker]} reads\n\n")
            alleles[marker] = {}
            continue
        lowest_allele_reads = min(alleles[marker].values())
        if lowest_allele_reads < min_reads_lowest:
            report_or_log(reportfile,
                f"Sample {tag} is not suitable for marker {marker}:\n"
                f"lowest allele has only {lowest_allele_reads} reads\n\n")
            alleles[marker] = {}
            continue
        expect = get_max_expected_alleles(max_alleles, marker, library)
        if len(alleles[marker]) > expect:
            allele_order = sorted(alleles[marker], key=lambda x: -alleles[marker][x])
            top_noise[marker] = [allele_order[expect], alleles[marker][allele_order[expect]]]
            alleles[marker] = {x: alleles[marker][x] for x in allele_order[:expect]}
        if top_noise[marker][1] > top_allele[marker] * (max_noise_pct / 100):
            highest_nonallele_pct = 100 * top_noise[marker][1] / top_allele[marker]
            report_or_log(reportfile,
                f"Sample {tag} is not suitable for marker {marker}:\n"
                f"highest non-allele is {highest_nonallele_pct:.1f}% of the highest allele\n")
            if reportfile:
                for allele in sorted(alleles[marker], key=lambda x: -alleles[marker][x]):
                    seq = allele if seqformat is None else ensure_sequence_format(
                            allele, seqformat, library=library, marker=marker)
                    reportfile.write(f"{alleles[marker][allele]}\tALLELE\t{seq}\n")
                seq = top_noise[marker][0] if seqformat is None else ensure_sequence_format(
                        top_noise[marker][0], seqformat, library=library, marker=marker)
                reportfile.write(f"{top_noise[marker][1]}\tNOISE\t{seq}\n")
            noisy_markers += 1
            alleles[marker] = {}

    # Drop this sample completely if it has too many noisy markers.
    if noisy_markers > (len(alleles) * max_noisy if max_noisy < 1 else max_noisy):
        report_or_log(reportfile, f"Sample {tag} appears to be contaminated!\n\n")
        return

    # The sample is OK, write out its alleles.
    for marker in alleles:
        for allele in sorted(alleles[marker], key=lambda x: -alleles[marker][x]):
            seq = allele if seqformat is None else ensure_sequence_format(
                allele, seqformat, library=library, marker=marker)
            outfile.write("\t".join((tag, marker, str(alleles[marker][allele]), seq)) + "\n")
#find_alleles_sample


def add_arguments(parser):
    outputgroup = add_bulk_input_output_args(parser, output_group=True)
    outputgroup.add_argument("-R", "--report", metavar="FILE",
        help="file to write a detailed report to (default: log briefly on stderr)")

    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--min-allele-pct", metavar="PCT",
        type=float, default=_DEF_MIN_ALLELE_PCT,
        help="call heterozygous if the second allele is at least this "
             "percentage of the highest allele of a marker (default: %(default)s)")
    filtergroup.add_argument("-M", "--max-noise-pct", metavar="PCT",
        type=float, default=_DEF_MAX_NOISE_PCT,
        help="a sample is considered contaminated/unsuitable for a marker if "
             "the highest non-allelic sequence is at least this percentage of "
             "the highest allele of that marker (default: %(default)s)")
    filtergroup.add_argument("-n", "--min-reads", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_READS,
        help="require at least this number of reads for the highest allele "
             "of each marker (default: %(default)s)")
    filtergroup.add_argument("-N", "--min-reads-lowest", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_READS_LOWEST,
        help="require at least this number of reads for the lowest allele "
             "of each marker (default: %(default)s)")
    filtergroup.add_argument("-a", "--max-alleles", metavar="N", type=pos_int_arg,
        help="allow no more than this number of alleles per marker; if unspecified, the amounts "
             "given in the library file are used, which have a default value of 1 for markers on "
             "the mitochondrial genome and Y chromosome, or 2 otherwise")
    filtergroup.add_argument("-x", "--max-noisy", metavar="X",
        type=float, default=_DEF_MAX_NOISY,
        help="entirely reject a sample if more than this fraction of markers (if less than 1) or "
             "absolute number of markers (if 1 or more) have a high non-allelic sequence "
             "(default: %(default)s)")
    add_sequence_format_args(parser)
#add_arguments


def run(args):
    samples_in = iterate_input_samples(args.infiles, args.tag_expr, args.tag_format)
    with std_or_open(args.outfile, "wt") as outfile, \
         std_or_open(args.report, "wt") as reportfile:
        find_alleles(samples_in, outfile, reportfile, args.min_reads, args.min_reads_lowest,
                     args.min_allele_pct, args.max_noise_pct, args.max_alleles, args.max_noisy,
                     args.sequence_format, args.library)
#run
