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
Estimate allele-centric background noise profiles (means) from reference
samples.

Compute a profile of recurring background noise for each unique allele
in the database of reference samples.  The profiles obtained can be used
by bgcorrect to filter background noise from samples.
"""
import argparse
import math
import sys
import time
#import numpy as np  # Only imported when actually running this tool.

from errno import EPIPE

from ..lib.cli import add_sequence_format_args, add_input_output_args, get_input_output_files,\
                      add_allele_detection_args, pos_int_arg
from ..lib.io import get_sample_data, parse_allelelist, try_write_pipe
from ..lib.noise import load_profiles
from ..lib.util import nnls

__version__ = "1.2.0"


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

# Default minimum number of heterozygous genotypes for each true allele.
# This value can be overridden by the -g command line option.
_DEF_MIN_GENOTYPES = 3


def solve_profile_mixture_single(samples, genotypes, n, starting_values={},
                                 variance=False, reportfile=None):
    """
    Solve the single-allele profiles of n true alleles:
      Profile of A:  [100,   5,  50, 25, ...]
      Profile of B:  [ 10, 100,  20,  0, ...]
      Profile of C:  [ 10,  30, 100, 30, ...]

    Given a list of observed profiles and known genotypes.
    The first n elements in each profile should be the true alleles.

    For each sample i, set genotypes[i] to a list of indices < n of the
    true allele.

    If variance=True, return a tuple (profile_means, profile_variances).
    This is an experimental feature and currently does not calculate the
    actual variance but rather a non-standard measure of variation
    inspired on variance.

    If reportfile is a writable handle, diagnostic/progress information
    is written to it.
    """
    num_samples = len(samples)
    profile_size = len(samples[0])

    A = np.empty((n, n))  # Will zero it at loop start.
    P = np.zeros((n, profile_size))
    C = np.zeros((n, profile_size))

    # Assume the true alleles do not have cross contributions at first.
    np.fill_diagonal(P, 100.)

    # Enter starting values into P.
    for x, y in starting_values:
        P[x, y] = starting_values[x, y]

    # Enter the samples into C.
    for i in range(num_samples):
        for j in genotypes[i]:
            try:
                # Compute factor to rescale such that true allele j has
                # 100 reads, and then divide by the number of true
                # alleles to make sure heterozygotes are not 'counted
                # twice' w.r.t. homozygotes.
                scale_factor = 100 / samples[i][j] / len(genotypes[i])
            except ZeroDivisionError:
                if reportfile:
                    try_write_pipe(reportfile, "Sample %i does not have allele %i\n" % (i, j))
                continue
            C[j, :] += [x * scale_factor for x in samples[i]]


    # Iteratively refine the goodness of fit to the data.
    prev_score = cur_score = sys.float_info.max
    for v in range(200):  # max 200 iterations here

        # Fill in A.
        A[:, :] = 0
        for i in range(num_samples):
            if len(genotypes[i]) == 1:
                # Shortcut for homozygotes.
                A[genotypes[i][0], genotypes[i][0]] += 1
                continue

            # Estimate allele balance in this sample based on the
            # current profiles.
            Px = P[genotypes[i], :][:, genotypes[i]]
            Cx = np.array(([samples[i][x] for x in genotypes[i]],))
            Ax = (100 / Cx / len(genotypes[i])).T @ nnls(Px.T, Cx.T).T

            # Update A with the values in Ax.
            for j in range(len(genotypes[i])):
                for k in range(len(genotypes[i])):
                    A[genotypes[i][j], genotypes[i][k]] += Ax[j, k]

        # Compute best-fitting profiles.
        # Doing this with the same nonnegative least squares method.
        E = A.T @ A
        F = A.T @ C
        prev_scorex = cur_scorex = sys.float_info.max
        for w in range(200):
            for p in range(n):
                if not E[p, p]:
                    # This would be an utter failure, but let's be safe.
                    if reportfile:
                        try_write_pipe(reportfile,
                            "%4i - No samples appear to have allele %i\n" % (v, p))
                    P[p, :p] = 0
                    P[p, p+1:] = 0
                    continue
                tmp = (F[p, :] - E[p:p+1, :p] @ P[:p, :] - E[p:p+1, p+1:] @ P[p+1:, :]) / E[p, p]
                tmp[tmp < 0] = 0
                tmp[0, p] = 100
                P[p, :] = tmp
            prev_scorex = cur_scorex
            cur_scorex = np.square(C - A @ P).sum()
            score_changex = (prev_scorex - cur_scorex) / prev_scorex
            if not cur_scorex or score_changex < 0.0001:
                break

        # Check whether profiles have converged.
        prev_score = cur_score
        cur_score = np.square(C - A @ P).sum()
        score_change = (prev_score - cur_score) / prev_score
        if v and reportfile:
            try_write_pipe(reportfile, "%4i %15.6f %15.6f %6.2f\n" %
                (v, cur_score, prev_score - cur_score, 100 * score_change))
        elif reportfile:
            try_write_pipe(reportfile, "%4i %15.6f\n" % (v, cur_score))
        if not cur_score or score_change < 0.0001:
            break

    if variance:
        # Variance estimation...
        # Going to solve A * V = C for V.  This time with C the
        # piecewise squared deviations of the samples from P.
        # We will include one row in A and C per allele per sample.
        if reportfile:
            if sum(map(len, genotypes)) - len(genotypes):
                try_write_pipe(reportfile,
                    "Computing variances...\n"
                    "EXPERIMENAL feature! The values produced may give a "
                    "sense of the amount of variation,\nbut should not be "
                    "used in further computations that expect true variances!")
            else:
                # With only homozygotes, the variances thus computed are
                # population variances.  It is more appropriate to
                # compute the sample variance, but how?
                try_write_pipe(reportfile,
                    "Computing variances...\n"
                    "EXPERIMENAL feature! The values produced are population "
                    "variances.\nThis may (or may not) change to sample "
                    "variance in a future version. Use with care.")

        # Fill in A and C.
        A = []
        C = []
        for i in range(num_samples):
            # Get relevant profiles for this sample.
            Px = P[genotypes[i], :]
            Cx = np.array((samples[i],))
            scale_factors = (100 / Cx[:, genotypes[i]] / len(genotypes[i])).T

            # Estimate allele balance in this sample based on the
            # current profiles.
            if len(genotypes[i]) == 1:
                # Shortcut for homozygotes.
                Ax = np.array(1., ndmin=2)
            else:
                Ax = scale_factors @ nnls(Px[:, genotypes[i]].T, Cx[:, genotypes[i]].T).T
            Cx = np.square(scale_factors @ Cx - Ax @ Px)

            # Update A and C with the values in Ax and the squared
            # deviations of this sample, respectively.
            for j in range(len(genotypes[i])):
                C.append(Cx[j, :].tolist())
                A.append([0.0] * n)
                for k in range(len(genotypes[i])):
                    A[-1][genotypes[i][k]] = Ax[j, k] ** 2
                    #A[-1][genotypes[i][k]] = Ax[j, k]
        A = np.array(A)
        C = np.array(C)

        # Compute best-fitting profiles.
        # Doing this with the same nonnegative least squares method.
        V = np.zeros(P.shape)
        E = A.T @ A
        #E = np.diagflat(A.sum(0)-1)
        F = A.T @ C
        prev_scorex = cur_scorex = sys.float_info.max
        for w in range(200):
            for p in range(n):
                if not E[p, p]:
                    V[p, :] = 0
                    continue
                tmp = (F[p, :] - E[p:p+1, :p] @ V[:p, :] - E[p:p+1, p+1:] @ V[p+1:, :]) / E[p, p]
                tmp[tmp < 0] = 0
                tmp[0, p] = 0  # No variance for actual allele.
                tmp[P[p, :] == 0] = 0  # No variance for zero means.
                V[p, :] = tmp
            prev_scorex = cur_scorex
            cur_scorex = np.square(C - A @ V).sum()
            score_changex = (prev_scorex - cur_scorex) / prev_scorex
            if not cur_scorex or score_changex < 0.0001:
                break
        return P.tolist(), V.tolist()
    else:
        return P.tolist()
#solve_profile_mixture_single


def filter_data(allelelist, min_samples, min_genotypes):
    if min_samples <= 1 and min_genotypes <= 1:
        return

    marker_names = set()
    for tag in allelelist:
        marker_names.update(allelelist[tag])

    for marker in marker_names:
        # Get a sample count of each true allele of this marker
        # and a sample count of each unique genotype of each allele.
        true_alleles = {}
        for tag in allelelist:
            if marker not in allelelist[tag]:
                continue
            for true_allele in allelelist[tag][marker]:
                try:
                    true_alleles[true_allele][0] += 1
                except KeyError:
                    true_alleles[true_allele] = [1, {}]
                if len(allelelist[tag][marker]) == 1:  # Got homozygote!
                    true_alleles[true_allele][1] = None
                elif true_alleles[true_allele][1] is not None:
                    genotype = "\t".join(sorted(allelelist[tag][marker]))
                    try:
                        true_alleles[true_allele][1][genotype] += 1
                    except KeyError:
                        true_alleles[true_allele][1][genotype] = 1

        # Drop any alleles that occur in less than min_samples samples
        # or in less than min_genotypes unique heterozygous genotypes
        # (by dropping the sample for this marker completely).
        repeat = True
        while repeat:
            repeat = False
            for true_allele in true_alleles:
                if 0 < true_alleles[true_allele][0] < min_samples or (
                        true_alleles[true_allele][1] is not None and
                        0 < len(true_alleles[true_allele][1]) < min_genotypes):
                    repeat = True
                    for tag in allelelist:
                        if marker not in allelelist[tag]:
                            continue
                        if true_allele in allelelist[tag][marker]:
                            genotype = "\t".join(sorted(allelelist[tag][marker]))
                            for allele in allelelist[tag][marker]:
                                true_alleles[allele][0] -= 1
                                if true_alleles[allele][1] is not None:
                                    true_alleles[allele][1][genotype] -= 1
                                    if not true_alleles[allele][1][genotype]:
                                        del true_alleles[allele][1][genotype]
                            del allelelist[tag][marker]
#filter_data


def add_sample_data(data, sample_data, sample_alleles, min_pct, min_abs, tag, combine_strands):
    # Make sure the true alleles of this sample are added to data.
    # Also compute the allele-specific inclusion thresholds for noise.
    thresholds = {}
    for marker in sample_alleles:
        if not sample_alleles[marker]:
            continue
        if marker not in data:
            data[marker] = {
                "profiles": {
                    "true alleles": 0,
                    "alleles": []},
                "allele_counts": {},
                "genotypes": []}
            if combine_strands:
                data[marker]["profiles"]["profiles_total"] = []
            else:
                data[marker]["profiles"]["profiles_forward"] = []
                data[marker]["profiles"]["profiles_reverse"] = []
        p = data[marker]["profiles"]
        if combine_strands:
            p["profiles_total"].append([0] * len(p["alleles"]))
        else:
            p["profiles_forward"].append([0] * len(p["alleles"]))
            p["profiles_reverse"].append([0] * len(p["alleles"]))
        data[marker]["genotypes"].append([])
        thresholds[marker] = {}
        for allele in sample_alleles[marker]:
            try:
                if 0 in sample_data[marker, allele]:
                    raise ValueError("Allele %s of marker %s has 0 reads in sample %s!" %
                        (allele, marker, tag))
            except KeyError:
                raise ValueError(
                    "Missing allele %s of marker %s in sample %s!" % (allele, marker, tag))
            try:
                i = p["alleles"].index(allele)
            except ValueError:
                i = len(p["alleles"])
                p["alleles"].append(allele)
                if combine_strands:
                    for profile in p["profiles_total"]:
                        profile.append(0)
                else:
                    for profile in p["profiles_forward"]:
                        profile.append(0)
                    for profile in p["profiles_reverse"]:
                        profile.append(0)
                for gi in data[marker]["allele_counts"]:
                    data[marker]["allele_counts"][gi].append(0)
            if i not in data[marker]["allele_counts"]:
                data[marker]["allele_counts"][i] = [0] * len(p["alleles"])
                p["true alleles"] += 1
            data[marker]["genotypes"][-1].append(i)
            thresholds[marker][i] = [math.ceil(x * min_pct / 100) for x in
                                     sample_data[marker, allele]]

    # Now enter the read counts into data and check the thresholds.
    for marker, allele in sample_data:
        if marker not in sample_alleles or not sample_alleles[marker]:
            # Sample does not participate in this marker (no alleles).
            continue

        p = data[marker]["profiles"]
        try:
            i = p["alleles"].index(allele)
        except ValueError:
            p["alleles"].append(allele)
            if combine_strands:
                for profile in p["profiles_total"]:
                    profile.append(0)
            else:
                for profile in p["profiles_forward"]:
                    profile.append(0)
                for profile in p["profiles_reverse"]:
                    profile.append(0)
            for gi in data[marker]["allele_counts"]:
                data[marker]["allele_counts"][gi].append(0)
            i = -1
        if combine_strands:
            p["profiles_total"][-1][i] = sample_data[marker, allele][0]
        else:
            p["profiles_forward"][-1][i] = sample_data[marker, allele][0]
            p["profiles_reverse"][-1][i] = sample_data[marker, allele][1]

        for gi in thresholds[marker]:
            if sum(count >= max(min_abs, threshold)
                   for count, threshold in
                   zip(sample_data[marker, allele], thresholds[marker][gi])):
                data[marker]["allele_counts"][gi][i] += 1
#add_sample_data


def preprocess_data(data, min_sample_pct, combine_strands):
    """
    Drop any sequence that is less than threshold_pct percent of the
    highest allele in more than 100-min_sample_pct of the samples with
    any particular true allele.

    The data is re-ordered as well, to ensure that the true alleles are
    placed before any other sequences in the profiles.
    """
    for marker in data:
        p = data[marker]["profiles"]
        counts = data[marker]["allele_counts"]
        thresholds = {i: min_sample_pct * counts[i][i] / 100 for i in counts}
        order = sorted(counts)
        for i in range(len(p["alleles"])):
            if i in counts:
                continue
            for gi in counts:
                if counts[gi][i] >= thresholds[gi]:
                    order.append(i)
                    break
        p["alleles"] = [p["alleles"][i] for i in order]
        if combine_strands:
            p["profiles_total"] = [[x[i] for i in order] for x in p["profiles_total"]]
        else:
            p["profiles_forward"] = [[x[i] for i in order] for x in p["profiles_forward"]]
            p["profiles_reverse"] = [[x[i] for i in order] for x in p["profiles_reverse"]]
        data[marker]["genotypes"] = [
            [order.index(y) for y in x] for x in data[marker]["genotypes"]]
        del data[marker]["allele_counts"]
#preprocess_data


def generate_profiles(samples_in, outfile, reportfile, allelefile, annotation_column, min_pct,
                      min_abs, min_samples, min_sample_pct, min_genotypes, library, profiles_in,
                      marker, homozygotes, combine_strands):
    if reportfile:
        t0 = time.time()

    # Parse allele list.
    allelelist = {} if allelefile is None else parse_allelelist(
        allelefile, convert="raw", library=library)

    # Read sample data.
    sample_data = {}
    get_sample_data(samples_in, lambda tag, data: sample_data.update({tag: data}),
                    allelelist, annotation_column, "raw", library, marker, homozygotes,
                    drop_special_seq=True, combine_strands=combine_strands)

    # Ensure minimum number of samples and genotypes per allele.
    allelelist = {tag: allelelist[tag] for tag in sample_data}
    filter_data(allelelist, min_samples, min_genotypes)

    # Combine data from all samples.  This takes most time.
    data = {}
    for tag in tuple(sample_data):
        add_sample_data(data, sample_data[tag], allelelist[tag], min_pct, min_abs, tag,
            combine_strands)
        del sample_data[tag]

    # Filter insignificant background products.
    preprocess_data(data, min_sample_pct, combine_strands)

    # Load starting profiles.
    start_profiles = {} if profiles_in is None else load_profiles(profiles_in, library)

    if reportfile:
        t1 = time.time()
        try_write_pipe(reportfile, "Data loading and filtering took %f seconds\n" % (t1 - t0))

    if combine_strands:
        outfile.write("\t".join(("marker", "allele", "sequence", "tmean", "tools")) + "\n")
    else:
        outfile.write("\t".join(
            ("marker", "allele", "sequence", "fmean", "rmean", "tools")) + "\n")
    for marker in tuple(data):
        p = data[marker]["profiles"]
        profile_size = len(p["alleles"])

        # Read relevent starting values from the starting profiles.
        init = {
            (p["alleles"].index(allele), p["alleles"].index(seq)):
                start_profiles[marker][allele][seq]["total"] if combine_strands else (
                start_profiles[marker][allele][seq]["forward"],
                start_profiles[marker][allele][seq]["reverse"])
            for allele in start_profiles[marker] if allele in p["alleles"]
            for seq in start_profiles[marker][allele] if seq in p["alleles"]
        } if marker in start_profiles else {}

        # Solve for the profiles of the true alleles.
        if reportfile:
            try_write_pipe(reportfile, "Solving marker %s with n=%i, m=%i, k=%i\n" % (
                marker, p["true alleles"], profile_size,
                len(p["profiles_total" if combine_strands else "profiles_forward"])))
            t0 = time.time()

        if combine_strands:
            p["profiles_total"] = solve_profile_mixture_single(
                p["profiles_total"], data[marker]["genotypes"], p["true alleles"],
                starting_values=init, reportfile=reportfile)
        else:
            for i, direction in enumerate(("forward", "reverse")):
                if reportfile:
                    try_write_pipe(reportfile, "Solving %s read profiles\n" % direction)
                p["profiles_" + direction] = solve_profile_mixture_single(
                    p["profiles_" + direction], data[marker]["genotypes"], p["true alleles"],
                    starting_values={x: init[x][i] for x in init}, reportfile=reportfile)

        if reportfile:
            t1 = time.time()
            try_write_pipe(reportfile, "Solved marker %s in %f seconds\n" % (marker, t1 - t0))

        # Round to 3 digits to get rid of funny rounding effects.
        # This method is not that precise anyway.
        if combine_strands:
            for profile in p["profiles_total"]:
                for i in range(profile_size):
                    profile[i] = round(profile[i], 3)
        else:
            for profile in p["profiles_forward"]:
                for i in range(profile_size):
                    profile[i] = round(profile[i], 3)
            for profile in p["profiles_reverse"]:
                for i in range(profile_size):
                    profile[i] = round(profile[i], 3)

        for i in range(p["true alleles"]):
            for j in range(len(p["alleles"])):
                if combine_strands:
                    if not p["profiles_total"][i][j]:
                        continue
                    outfile.write("\t".join((marker, p["alleles"][i], p["alleles"][j],
                        str(p["profiles_total"][i][j]), "bgestimate")) + "\n")
                else:
                    if not p["profiles_forward"][i][j] + p["profiles_reverse"][i][j]:
                        continue
                    outfile.write("\t".join(
                        [marker, p["alleles"][i], p["alleles"][j]] +
                        list(map(str, [p["profiles_forward"][i][j],
                                  p["profiles_reverse"][i][j]])) +
                        ["bgestimate"]) + "\n")
        del data[marker]
#generate_profiles


def add_arguments(parser):
    add_input_output_args(parser, single_in=False, batch_support=False, report_out=True)
    add_allele_detection_args(parser)
    parser.add_argument("-C", "--combine-strands", action="store_true",
        help="if specified, noise profiles will be calculated for the total number of reads, "
             "instead of separately for either strand")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--min-pct", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to consider, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument("-n", "--min-abs", metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="minimum amount of background to consider, as an absolute "
             "number of reads for at least one orientation (default: %(default)s)")
    filtergroup.add_argument("-s", "--min-samples", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele (default: %(default)s)")
    filtergroup.add_argument("-S", "--min-sample-pct", metavar="PCT",
        type=float, default=_DEF_MIN_SAMPLE_PCT,
        help="require this minimum number of samples for each background "
             "product, as a percentage of the number of samples with a "
             "particular true allele (default: %(default)s)")
    filtergroup.add_argument("-g", "--min-genotypes", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_GENOTYPES,
        help="require this minimum number of unique heterozygous genotypes "
             "for each allele for which no homozygous samples are available "
             "(default: %(default)s)")
    filtergroup.add_argument("-p", "--profiles", metavar="FILE",
        type=argparse.FileType("tr", encoding="UTF-8"),
        help="use the given noise profiles file as a starting point")
    filtergroup.add_argument("-M", "--marker", metavar="MARKER", help="work only on MARKER")
    filtergroup.add_argument("-H", "--homozygotes", action="store_true",
        help="if specified, only homozygous samples will be considered")
    add_sequence_format_args(parser, default_format="raw", force=True)
#add_arguments


def run(args):
    # Import numpy now.
    global np
    import numpy as np

    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output of another program")
    try:
        generate_profiles(files[0], files[1], args.report, args.allelelist, args.annotation_column,
                          args.min_pct, args.min_abs, args.min_samples, args.min_sample_pct,
                          args.min_genotypes, args.library, args.profiles, args.marker,
                          args.homozygotes, args.combine_strands)
    except IOError as e:
        if e.errno == EPIPE:
            try:
                try_write_pipe(args.report, "Stopped early because the output was closed.\n")
            except:
                pass
            return
        raise
#run
