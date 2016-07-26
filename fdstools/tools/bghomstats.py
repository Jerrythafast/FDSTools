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
Compute allele-centric statistics for background noise in homozygous
reference samples (min, max, mean, sample variance).

Compute a profile of recurring background noise for each unique allele
in the database of reference samples.  The profiles obtained can be used
by bgcorrect to filter background noise from samples.  If many reference
samples are heterozygous (as is usually the case with forensic STR
markers), it is preferable to use bgestimate instead, since it can
handle heterozygous samples as well.
"""
from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_allele_detection_args, parse_allelelist,\
                  get_sample_data, add_sequence_format_args, adjust_stats,\
                  add_random_subsampling_args

__version__ = "1.0.0"


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


def add_sample_data(data, sample_data, sample_alleles, min_pct, min_abs):
    # Check presence of all alleles.
    for marker in sample_alleles:
        allele = sample_alleles[marker]
        if (marker, allele) not in sample_data:
            raise ValueError(
                "Missing allele %s of marker %s!" % (allele, marker))
        elif 0 in sample_data[marker, allele]:
            raise ValueError(
                "Allele %s of marker %s has 0 reads!" % (allele, marker))

    # Enter the read counts into data and check the thresholds.
    for marker, sequence in sample_data:
        if marker not in sample_alleles:
            # Sample does not participate in this marker.
            continue
        allele = sample_alleles[marker]
        factors = [100./x for x in sample_data[marker, allele]]
        if (marker, allele) not in data:
            data[marker, allele] = {}
        if sequence not in data[marker, allele]:
            data[marker, allele][sequence] = [None, None, 0]
        for direction in (0, 1):
            data[marker, allele][sequence][direction] = adjust_stats(
                sample_data[marker, sequence][direction] * factors[direction],
                data[marker, allele][sequence][direction])
        if sum(count >= min_abs and count*factor >= min_pct
               for count, factor in
               zip(sample_data[marker, sequence], factors)):
            data[marker, allele][sequence][2] += 1
#add_sample_data


def filter_data(data, min_samples, min_sample_pct):
    """
    Remove all alleles from data that have less than min_samples samples
    and remove all stats of sequences that don't pass the detection
    thresholds in at least min_sample_pct per cent of the samples with a
    particular allele.  Also add explicit zeros to the stats of the
    sequences that were not seen in all samples with a given allele.
    """
    for marker, allele in data.keys():
        if data[marker, allele][allele][2] < min_samples:
            del data[marker, allele]
            continue
        factor = 100./data[marker, allele][allele][2]
        for sequence in data[marker, allele].keys():
            if data[marker, allele][sequence][2] * factor < min_sample_pct:
                del data[marker, allele][sequence]
                continue
            for i in range(data[marker, allele][sequence][0]["n"],
                           data[marker, allele][allele][2]):
                for direction in (0, 1):
                    adjust_stats(0, data[marker, allele][sequence][direction])
#filter_data


def compute_stats(samples_in, outfile, allelefile, annotation_column, min_pct,
                  min_abs, min_samples, min_sample_pct, seqformat, library,
                  marker, limit_reads, drop_samples):

    # Parse allele list.
    allelelist = {} if allelefile is None \
                    else parse_allelelist(allelefile, seqformat, library)

    # Read sample data.
    data = {}
    get_sample_data(
        samples_in,
        lambda tag, sample_data: add_sample_data(
            data, sample_data,
            {m: allelelist[tag][m].pop() for m in allelelist[tag]},
            min_pct, min_abs),
        allelelist, annotation_column, seqformat, library, marker, True,
        limit_reads, drop_samples, True)

    # Ensure minimum number of samples per allele and filter
    # insignificant background products.
    filter_data(data, min_samples, min_sample_pct)

    outfile.write("\t".join(["marker", "allele", "sequence", "n", "fmin",
                     "fmax", "fmean", "fvariance", "rmin", "rmax", "rmean",
                     "rvariance", "tool"]) + "\n")
    for marker, allele in data:
        for sequence in data[marker, allele]:
            outfile.write("\t".join([marker, allele, sequence] + [
                "%.3g" % x if abs(x) > 0.0000000001 else "0" for x in (
                    data[marker, allele][sequence][0]["n"],
                    data[marker, allele][sequence][0]["min"],
                    data[marker, allele][sequence][0]["max"],
                    data[marker, allele][sequence][0]["mean"],
                    data[marker, allele][sequence][0]["variance"],
                    data[marker, allele][sequence][1]["min"],
                    data[marker, allele][sequence][1]["max"],
                    data[marker, allele][sequence][1]["mean"],
                    data[marker, allele][sequence][1]["variance"])] +
                ["bghomstats"]) + "\n")
#compute_stats


def add_arguments(parser):
    add_input_output_args(parser)
    add_allele_detection_args(parser)
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-m', '--min-pct', metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to consider, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument('-n', '--min-abs', metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="minimum amount of background to consider, as an absolute "
             "number of reads (default: %(default)s)")
    filtergroup.add_argument('-s', '--min-samples', metavar="N",
        type=pos_int_arg,
        default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele "
             "(default: %(default)s)")
    filtergroup.add_argument('-S', '--min-sample-pct', metavar="PCT",
        type=float,
        default=_DEF_MIN_SAMPLE_PCT,
        help="require this minimum number of samples for each background "
             "product, as a percentage of the number of samples with a "
             "particular true allele (default: %(default)s)")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    add_sequence_format_args(parser)
    add_random_subsampling_args(parser)
#add_arguments


def run(args):
    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    compute_stats(files[0], files[1], args.allelelist, args.annotation_column,
                  args.min_pct, args.min_abs, args.min_samples,
                  args.min_sample_pct, args.sequence_format, args.library,
                  args.marker, args.limit_reads, args.drop_samples)
#run
