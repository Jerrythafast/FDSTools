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
Compute noise ratios for all noise detected in homozygous reference
samples.

With this tool, separate data points are produced for each sample, which
can be visualised using 'fdstools vis bgraw'.  Use bghomstats or
bgestimate to compute aggregate statistics on noise instead.
"""
from ..lib.cli import add_bulk_input_output_args, add_sequence_format_args, \
                      iterate_input_samples, pos_int_arg, std_or_open
from ..lib.io import parse_allelelist, read_sample_data_file
from ..lib.seq import ensure_sequence_format

__version__ = "1.1.0"


# Default values for parameters are specified below.

# Default minimum amount of background to consider, as a percentage of
# the highest allele.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT = 0.5

# Default minimum number of reads to consider.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS = 5

# Default minimum number of samples for each true allele.
# This value can be overridden by the -s command line option.
_DEF_MIN_SAMPLES = 2

# Default minimum number of samples required for each background product
# to be included in the analysis, as a percentage of the number of
# samples with a certain true allele.
# This value can be overridden by the -S command line option.
_DEF_MIN_SAMPLE_PCT = 80.


def add_sample_data(data, sample_data, sample_alleles, min_pct, min_abs, tag, combine_strands):
    # Check presence of all alleles.
    for marker in sample_alleles:
        allele = sample_alleles[marker]
        try:
            reads = sample_data[marker, allele]
        except KeyError:
            raise ValueError(
                "Missing allele %s of marker %s in sample %s!" % (allele, marker, tag))
        if 0 in reads:
            raise ValueError(
                "Allele %s of marker %s has 0 reads in sample %s!" % (allele, marker, tag))

    # Enter the read counts into data and check the thresholds.
    for marker, sequence in sample_data:
        if marker not in sample_alleles:
            # Sample does not participate in this marker.
            continue
        allele = sample_alleles[marker]

        reads = sample_data[marker, sequence]
        if not combine_strands:
            reads.append(sum(reads))
        factors = [100 / (x or 1) for x in sample_data[marker, allele]]
        if (marker, allele) not in data:
            data[marker, allele] = {}
        if sequence not in data[marker, allele]:
            data[marker, allele][sequence] = {
                "tag": [],
                "total": [],
                "tnoise": [],
                "passed_filter": 0}
            if not combine_strands:
                data[marker, allele][sequence]["forward"] = []
                data[marker, allele][sequence]["reverse"] = []
                data[marker, allele][sequence]["fnoise"] = []
                data[marker, allele][sequence]["rnoise"] = []

        this_data = data[marker, allele][sequence]
        this_data["tag"].append(tag)
        this_data["total"].append(reads[-1])
        this_data["tnoise"].append(reads[-1] * factors[-1])

        if not combine_strands:
            this_data["forward"].append(reads[0])
            this_data["reverse"].append(reads[1])
            this_data["fnoise"].append(reads[0] * factors[0])
            this_data["rnoise"].append(reads[1] * factors[1])

        if combine_strands:
            if reads[0] >= min_abs and reads[0] * factors[0] >= min_pct:
                this_data["passed_filter"] += 1
        elif any(count >= min_abs and count * factor >= min_pct
               for count, factor in zip(reads[:2], factors[:2])):
            this_data["passed_filter"] += 1
#add_sample_data


def filter_data(data, min_samples, min_sample_pct):
    """
    Remove all alleles from data that have less than min_samples samples
    and remove all data of sequences that don't pass the detection
    thresholds in at least min_sample_pct percent of the samples with a
    particular allele.
    """
    for marker, allele in tuple(data):
        if data[marker, allele][allele]["passed_filter"] < min_samples:
            del data[marker, allele]
            continue
        factor = 100 / data[marker, allele][allele]["passed_filter"]
        for sequence in tuple(data[marker, allele]):
            if data[marker, allele][sequence]["passed_filter"] * factor < min_sample_pct:
                del data[marker, allele][sequence]
                continue
#filter_data


def compute_ratios(samples_in, outfile, allelefile, min_pct, min_abs, min_samples, min_sample_pct,
                   no_alleles, seqformat, library, only_marker, combine_strands):
    if allelefile is None:
        # No need to convert input sequence data.
        allelelist = {}
        input_format = None
    else:
        # Parse allele list, convert all input to raw sequences.
        input_format = "raw"
        allelelist = parse_allelelist(allelefile, convert="raw", library=library)

    # Read sample data.
    data = {}
    for tag, infile in samples_in:
        # Fill allelelist using "allele" flag if no allelelist file was given.
        alleles_out = allelelist.setdefault(tag, {}) if allelefile is None else None
        if tag not in allelelist:
            continue  # No need to read this sample.
        sample_data = read_sample_data_file(infile, "raw", library, drop_special_seq=True,
            combine_strands=combine_strands, only_marker=only_marker, alleles_out=alleles_out)
        sample_alleles = {
            marker: alleles.pop()
            for marker, alleles in allelelist[tag].items()
            if len(alleles) == 1
            and (only_marker is None or marker == only_marker)}
        add_sample_data(data, sample_data, sample_alleles, min_pct, min_abs, tag, combine_strands)

    # Ensure minimum number of samples per allele and filter
    # insignificant background products.
    filter_data(data, min_samples, min_sample_pct)

    if combine_strands:
        outfile.write(
            "\t".join(("sample", "marker", "allele", "sequence", "total", "tnoise")) + "\n")
    else:
        outfile.write("\t".join(("sample", "marker", "allele", "sequence",
            "forward", "reverse", "total", "fnoise", "rnoise", "tnoise")) + "\n")
    for marker, raw_allele in data:
        allele = ensure_sequence_format(
            raw_allele, seqformat, from_format=input_format, library=library, marker=marker)
        for raw_sequence, values in data[marker, raw_allele].items():
            if no_alleles and raw_sequence == raw_allele:
                continue
            sequence = ensure_sequence_format(
                raw_sequence, seqformat, from_format=input_format, library=library, marker=marker)
            for i, tag in enumerate(values["tag"]):
                if combine_strands:
                    noise = values["tnoise"][i]
                    outfile.write("\t".join((tag, marker, allele, sequence,
                        str(values["total"][i]),
                        "%.3g" % noise if abs(noise) > 0.0000000001 else "0")) + "\n")
                else:
                    outfile.write("\t".join([tag, marker, allele, sequence] + list(map(str, (
                            values["forward"][i],
                            values["reverse"][i],
                            values["total"][i]))) + [
                        "%.3g" % x if abs(x) > 0.0000000001 else "0" for x in (
                            values["fnoise"][i],
                            values["rnoise"][i],
                            values["tnoise"][i])]) + "\n")
#compute_ratios


def add_arguments(parser):
    add_bulk_input_output_args(parser)
    parser.add_argument("-a", "--allelelist", metavar="ALLELEFILE",
        help="file containing a list of the true alleles of each sample (e.g., obtained from "
             "allelefinder); if not given, the 'allele' flag from the 'flags' column is used")
    parser.add_argument("-C", "--combine-strands", action="store_true",
        help="if specified, noise ratios will be calculated for the total number of reads, "
             "instead of separately for either strand")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--min-pct", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to consider, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument("-n", "--min-abs", metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="minimum amount of background to consider, as an absolute "
             "number of reads (default: %(default)s)")
    filtergroup.add_argument("-s", "--min-samples", metavar="N", type=pos_int_arg,
        default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele "
             "(default: %(default)s)")
    filtergroup.add_argument("-S", "--min-sample-pct", metavar="PCT", type=float,
        default=_DEF_MIN_SAMPLE_PCT,
        help="require this minimum number of samples for each background "
             "product, as a percentage of the number of samples with a "
             "particular true allele (default: %(default)s)")
    filtergroup.add_argument("-A", "--no-alleles", action="store_true",
        help="if specified, lines corresponding to the alleles are omitted from the output")
    filtergroup.add_argument("-M", "--marker", metavar="MARKER", help="work only on MARKER")
    add_sequence_format_args(parser, default_format="raw")
#add_arguments


def run(args):
    samples_in = iterate_input_samples(args.infiles, args.tag_expr, args.tag_format)
    with std_or_open(args.outfile, "wt") as outfile:
        compute_ratios(samples_in, outfile, args.allelelist, args.min_pct,
                       args.min_abs, args.min_samples, args.min_sample_pct, args.no_alleles,
                       args.sequence_format, args.library, args.marker, args.combine_strands)
#run
