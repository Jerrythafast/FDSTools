#!/usr/bin/env python
"""
Estimate allele-centric background noise profiles (means) from reference
samples.

Compute a profile of recurring background noise for each unique allele
in the database of reference samples.  The profiles obtained can be used
by bgcorrect to filter background noise from samples.
"""
import argparse
import sys
import time
import math
#import numpy as np  # Only imported when actually running this tool.

from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_allele_detection_args, nnls, add_sequence_format_args,\
                  parse_allelelist, parse_library, get_sample_data, \
                  add_random_subsampling_args

__version__ = "0.1dev"


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


def solve_profile_mixture(forward, reverse, genotypes, n, variance=False,
                          reportfile=None):
    if reportfile:
        reportfile.write("Solving forward read profiles\n")
    f = solve_profile_mixture_single(forward, genotypes, n,
            variance=variance, reportfile=reportfile)
    if reportfile:
        reportfile.write("Solving reverse read profiles\n")
    r = solve_profile_mixture_single(reverse, genotypes, n,
            variance=variance, reportfile=reportfile)
    if variance:
        return f[0], r[0], f[1], r[1]
    return f, r
#solve_profile_mixture


def solve_profile_mixture_single(samples, genotypes, n, variance=False,
                                 reportfile=None):
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
    import numpy as np
    num_samples = len(samples)
    profile_size = len(samples[0])

    A = np.matrix(np.empty([n, n]))  # Will zero it at loop start.
    P = np.matrix(np.zeros([n, profile_size]))
    C = np.matrix(np.zeros([n, profile_size]))

    # Assume the true alleles do not have cross contributions at first.
    np.fill_diagonal(P, 100.)

    # Enter the samples into C.
    for i in range(num_samples):
        for j in genotypes[i]:
            try:
                # Compute factor to rescale such that true allele j has
                # 100 reads, and then divide by the number of true
                # alleles to make sure heterozygotes are not 'counted
                # twice' w.r.t. homozygotes.
                scale_factor = (100.0/samples[i][j])/len(genotypes[i])
            except ZeroDivisionError:
                if reportfile:
                    reportfile.write(
                        "Sample %i does not have allele %i\n" % (i, j))
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
            Cx = np.matrix([[samples[i][x] for x in genotypes[i]]])
            Ax = (100./Cx/len(genotypes[i])).T * nnls(Px.T, Cx.T).T

            # Update A with the values in Ax.
            for j in range(len(genotypes[i])):
                for k in range(len(genotypes[i])):
                    A[genotypes[i][j], genotypes[i][k]] += Ax[j, k]

        # Compute best-fitting profiles.
        # Doing this with the same nonnegative least squares method.
        E = A.T * A
        F = A.T * C
        prev_scorex = cur_scorex = sys.float_info.max
        for w in range(200):
            for p in range(n):
                nn = list(range(p)) + list(range(p + 1, n))
                if not E[p, p]:
                    # This would be an utter failure, but let's be safe.
                    if reportfile:
                        reportfile.write(
                            "%4i - No samples appear to have allele %i\n" %
                            (v, p))
                    P[p, nn] = 0
                else:
                    tmp = (F[p, :] - E[p, nn] * P[nn, :]) / E[p, p]
                    tmp[tmp < 0] = 0
                    tmp[0, p] = 100
                    P[p, :] = tmp
            prev_scorex = cur_scorex
            cur_scorex = np.square(C - A * P).sum()
            score_changex = (prev_scorex-cur_scorex)/prev_scorex
            if not cur_scorex or score_changex < 0.0001:
                break

        # Check whether profiles have converged.
        prev_score = cur_score
        cur_score = np.square(C - A * P).sum()
        score_change = (prev_score-cur_score)/prev_score
        if v and reportfile:
            reportfile.write("%4i %15.6f %15.6f %6.2f\n" %
                        (v, cur_score, prev_score-cur_score, 100*score_change))
        elif reportfile:
            reportfile.write("%4i %15.6f\n" % (v, cur_score))
        if not cur_score or score_change < 0.0001:
            break

    if variance:
        # Variance estimation...
        # Going to solve A * V = C for V.  This time with C the
        # piecewise squared deviations of the samples from P.
        # We will include one row in A and C per allele per sample.

        # The variances thus computed are population variances.  It is
        # more appropriate to compute the sample variance, but how?
        if reportfile:
            if sum(map(len, genotypes))-len(genotypes):
                reportfile.write(
                    "Computing variances...\n"
                    "EXPERIMENAL feature! The values produced may give a "
                    "sense of the amount of variation,\nbut should not be "
                    "used in further computations that expect true variances!")
            else:
                reportfile.write(
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
            Cx = np.matrix(samples[i])
            scale_factors = (100./Cx[:, genotypes[i]]/len(genotypes[i])).T

            # Estimate allele balance in this sample based on the
            # current profiles.
            if len(genotypes[i]) == 1:
                # Shortcut for homozygotes.
                Ax = np.matrix([[1.]])
            else:
                Ax = scale_factors * \
                     nnls(Px[:, genotypes[i]].T, Cx[:, genotypes[i]].T).T
            Cx = np.square(scale_factors * Cx - Ax * Px)

            # Update A and C with the values in Ax and the squared
            # deviations of this sample, respectively.
            for j in range(len(genotypes[i])):
                C.append(Cx[j, :].tolist()[0])
                A.append([0.0] * n)
                for k in range(len(genotypes[i])):
                    A[-1][genotypes[i][k]] = Ax[j, k] ** 2
                    #A[-1][genotypes[i][k]] = Ax[j, k]
        A = np.matrix(A)
        C = np.matrix(C)

        # Compute best-fitting profiles.
        # Doing this with the same nonnegative least squares method.
        V = np.matrix(np.zeros(P.shape))
        E = A.T * A
        #E = np.diagflat(A.sum(0)-1)
        F = A.T * C
        prev_scorex = cur_scorex = sys.float_info.max
        for w in range(200):
            for p in range(n):
                nn = list(range(p)) + list(range(p + 1, n))
                if not E[p, p]:
                    V[p, :] = 0
                else:
                    tmp = (F[p, :] - E[p, nn] * V[nn, :]) / E[p, p]
                    tmp[tmp < 0] = 0
                    tmp[0, p] = 0  # No variance for actual allele.
                    tmp[P[p, :] == 0] = 0  # No variance for zero means.
                    V[p, :] = tmp
            prev_scorex = cur_scorex
            cur_scorex = np.square(C - A * V).sum()
            score_changex = (prev_scorex-cur_scorex)/prev_scorex
            if not cur_scorex or score_changex < 0.0001:
                break
        return P.tolist(), V.tolist()
    else:
        return P.tolist()
#solve_profile_mixture


def ensure_min_samples(allelelist, min_samples):
    if min_samples <= 1:
        return

    marker_names = set()
    for tag in allelelist:
        marker_names.update(allelelist[tag].keys())

    for marker in marker_names:
        # Get a sample count of each true allele of this marker.
        true_alleles = {}
        for tag in allelelist:
            if marker not in allelelist[tag]:
                continue
            for true_allele in allelelist[tag][marker]:
                try:
                    true_alleles[true_allele] += 1
                except KeyError:
                    true_alleles[true_allele] = 1

        # Drop any alleles that occur in less than min_samples samples
        # (by dropping the sample for this marker completely).
        repeat = True
        while repeat:
            repeat = False
            for true_allele in true_alleles:
                if 0 < true_alleles[true_allele] < min_samples:
                    repeat = True
                    for tag in allelelist:
                        if marker not in allelelist[tag]:
                            continue
                        if true_allele in allelelist[tag][marker]:
                            for allele in allelelist[tag][marker]:
                                true_alleles[allele] -= 1
                            del allelelist[tag][marker]
#ensure_min_samples


def add_sample_data(data, sample_data, sample_alleles, min_pct, min_abs):
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
                    "alleles": [],
                    "profiles_forward": [],
                    "profiles_reverse": []},
                "allele_counts": {},
                "genotypes": []}
        p = data[marker]["profiles"]
        p["profiles_forward"].append([0] * len(p["alleles"]))
        p["profiles_reverse"].append([0] * len(p["alleles"]))
        data[marker]["genotypes"].append([])
        thresholds[marker] = {}
        for allele in sample_alleles[marker]:
            if (marker, allele) not in sample_data:
                raise ValueError(
                    "Missing allele %s of marker %s!" % (allele, marker))
            elif 0 in sample_data[marker, allele]:
                raise ValueError(
                    "Allele %s of marker %s has 0 reads!" % (allele, marker))
            try:
                i = p["alleles"].index(allele)
            except ValueError:
                i = len(p["alleles"])
                p["alleles"].append(allele)
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
            thresholds[marker][i] = [math.ceil(x*min_pct/100.) for x in
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
            for profile in p["profiles_forward"]:
                profile.append(0)
            for profile in p["profiles_reverse"]:
                profile.append(0)
            for gi in data[marker]["allele_counts"]:
                data[marker]["allele_counts"][gi].append(0)
            i = -1
        p["profiles_forward"][-1][i] = sample_data[marker, allele][0]
        p["profiles_reverse"][-1][i] = sample_data[marker, allele][1]

        for gi in thresholds[marker]:
            if sum(count >= max(min_abs, threshold)
                   for count, threshold in
                   zip(sample_data[marker, allele], thresholds[marker][gi])):
                data[marker]["allele_counts"][gi][i] += 1
#add_sample_data


def preprocess_data(data, min_sample_pct):
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
        thresholds = {i: min_sample_pct * counts[i][i] / 100. for i in counts}
        order = sorted(counts)
        for i in range(len(p["alleles"])):
            if i in counts:
                continue
            for gi in counts:
                if counts[gi][i] >= thresholds[gi]:
                    order.append(i)
                    break
        p["alleles"] = [p["alleles"][i] for i in order]
        p["profiles_forward"] = [
            [x[i] for i in order] for x in p["profiles_forward"]]
        p["profiles_reverse"] = [
            [x[i] for i in order] for x in p["profiles_reverse"]]
        data[marker]["genotypes"] = [
            [order.index(y) for y in x] for x in data[marker]["genotypes"]]
        del data[marker]["allele_counts"]
#preprocess_data


def generate_profiles(samples_in, outfile, reportfile, allelefile,
                      annotation_column, min_pct, min_abs, min_samples,
                      min_sample_pct, seqformat, library, crosstab, marker,
                      homozygotes, limit_reads, drop_samples):
    if reportfile:
        t0 = time.time()

    # Parse library and allele list.
    library = parse_library(library) if library is not None else None
    allelelist = {} if allelefile is None \
                    else parse_allelelist(allelefile, seqformat, library)

    # Read sample data.
    sample_data = {}
    get_sample_data(samples_in,
                    lambda tag, data: sample_data.update({tag: data}),
                    allelelist, annotation_column, seqformat, library, marker,
                    homozygotes, limit_reads, drop_samples)

    # Ensure minimum number of samples per allele.
    allelelist = {tag: allelelist[tag] for tag in sample_data}
    ensure_min_samples(allelelist, min_samples)

    # Combine data from all samples.
    data = {}
    for tag in sample_data.keys():
        add_sample_data(data, sample_data[tag], allelelist[tag], min_pct,
                        min_abs)
        del sample_data[tag]

    # Filter insignificant background products.
    preprocess_data(data, min_sample_pct)

    if reportfile:
        t1 = time.time()
        reportfile.write("Data loading and filtering took %f seconds\n" %
                         (t1-t0))

    if not crosstab:
        outfile.write("\t".join(
            ["marker", "allele", "sequence", "fmean", "rmean"]) + "\n")
    for marker in data.keys():
        p = data[marker]["profiles"]
        profile_size = len(p["alleles"])

        # Solve for the profiles of the true alleles.
        if reportfile:
            reportfile.write("Solving marker %s with n=%i, m=%i, k=%i\n" %
                             (marker, p["true alleles"], profile_size,
                              len(p["profiles_forward"])))
            t0 = time.time()
        p["profiles_forward"], p["profiles_reverse"] = solve_profile_mixture(
                p["profiles_forward"],
                p["profiles_reverse"],
                data[marker]["genotypes"],
                p["true alleles"],
                reportfile=reportfile)
        if reportfile:
            t1 = time.time()
            reportfile.write("Solved marker %s in %f seconds\n" %
                             (marker, t1-t0))

        # Round to 3 digits to get rid of funny rounding effects.
        # This method is not that precise anyway.
        for profile in p["profiles_forward"]:
            for i in range(profile_size):
                profile[i] = round(profile[i], 3)
        for profile in p["profiles_reverse"]:
            for i in range(profile_size):
                profile[i] = round(profile[i], 3)

        if crosstab:
            # Cross-tabular output (profiles in rows)
            outfile.write("\t".join([marker, "0"] + p["alleles"]) + "\n")
            for i in range(p["true alleles"]):
                outfile.write("\t".join(
                    [marker, str(i+1)] + map(str, p["profiles_forward"][i])) +
                    "\n")
                outfile.write("\t".join(
                    [marker, str(-i-1)] + map(str, p["profiles_reverse"][i])) +
                    "\n")
        else:
            # Tab-separated columns format.
            for i in range(p["true alleles"]):
                for j in range(len(p["alleles"])):
                    if not (p["profiles_forward"][i][j] +
                            p["profiles_reverse"][i][j]):
                        continue
                    outfile.write("\t".join(
                        [marker, p["alleles"][i], p["alleles"][j]] +
                        map(str, [p["profiles_forward"][i][j],
                                  p["profiles_reverse"][i][j]])) + "\n")
        del data[marker]
#generate_profiles


def add_arguments(parser):
    add_input_output_args(parser, False, False, True)
    parser.add_argument('-C', '--cross-tabular', action="store_true",
        help="if specified, a space-efficient cross-tabular output format is "
             "used instead of the default tab-separated columns format")
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
        type=pos_int_arg, default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele "
             "(default: %(default)s)")
    filtergroup.add_argument('-S', '--min-sample-pct', metavar="PCT",
        type=float, default=_DEF_MIN_SAMPLE_PCT,
        help="require this minimum number of samples for each background "
             "product, as a percentage of the number of samples with a "
             "particular true allele (default: %(default)s)")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    filtergroup.add_argument('-H', '--homozygotes', action="store_true",
        help="if specified, only homozygous samples will be considered")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
    add_random_subsampling_args(parser)
#add_arguments


def run(args):
    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    generate_profiles(files[0], files[1], args.report, args.allelelist,
                      args.annotation_column, args.min_pct, args.min_abs,
                      args.min_samples, args.min_sample_pct,
                      args.sequence_format, args.library, args.cross_tabular,
                      args.marker, args.homozygotes, args.limit_reads,
                      args.drop_samples)
#run


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    try:
        add_arguments(parser)
        run(parser.parse_args())
    except OSError as error:
        parser.error(error)
#main


if __name__ == "__main__":
    main()
