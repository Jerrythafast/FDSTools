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
Train a stutter prediction model using homozygous reference samples.

The model obtained from this tool can be used by bgpredict to predict
background noise profiles of alleles for which no reference samples are
available.
"""

from  itertools import combinations
import re
#import numpy as np  # Only imported when actually running this tool.

from ..lib.noise import iter_overlapping_repeats, get_all_combinations, get_stutter_sequence
from ..lib.cli import add_bulk_input_output_args, add_sequence_format_args, \
                      iterate_input_samples, pos_int_arg
from ..lib.io import parse_allelelist, read_sample_data_file, std_or_open
from ..lib.seq import reverse_complement, get_repeat_pattern
from ..lib.util import nnls

__version__ = "1.2.1"


# Default values for parameters are specified below.

# Default minimum amount of background to consider, as a percentage of
# the highest allele.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT = 0

# Default minimum number of reads to consider.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS = 1

# Default minimum stutter percentage to fit a model through.
# This value can be overridden by the -p command line option.
_DEF_THRESHOLD_STUTTER_PCT = 0.5

# Default minimum number of samples for each true allele.
# This value can be overridden by the -s command line option.
_DEF_MIN_SAMPLES = 1

# Default minimum number of different repeat lengths per fit.
# This value can be overridden by the -L command line option.
_DEF_MIN_LENGTHS = 5

# Default degree of polynomials to fit.
# This value can be overridden by the -D command line option.
_DEF_DEGREE = 2

# Default maximum repeat unit length to investigate.
# This value can be overridden by the -u command line option.
_DEF_MAX_UNIT_LENGTH = 6

# Default minimum R2 score.
# This value can be overridden by the -t command line option.
_DEF_MIN_R2 = 0.5


def is_repeated_sequence(sequence):
    """Test whether a sequence consists of one repeated element."""
    for i in range(1, len(sequence)):
        if not len(sequence) % i and sequence == sequence[:i] * (len(sequence) // i):
            return True
    return False
#is_repeated_sequence


def lowest_cyclic_variant(sequence):
    """Return the lowest-sorting cyclic variant of sequence."""
    bestseq = sequence
    seqseq = sequence * 2
    qesqes = reverse_complement(sequence) * 2
    for i in range(len(sequence)):
        if seqseq[i : i + len(sequence)] < bestseq:
            bestseq = seqseq[i : i + len(sequence)]
        if qesqes[i : i + len(sequence)] < bestseq:
            bestseq = qesqes[i : i + len(sequence)]
    return bestseq
#lowest_cyclic_variant


def get_unique_repeats(maxlen, out=None, prefix=""):
    """Return the set of unique repeat sequences of length <= maxlen."""
    if out == None:
        out = set()
    for base in "ACGT":
        sequence = prefix + base
        if is_repeated_sequence(sequence):
            continue
        out.add(lowest_cyclic_variant(sequence))
        if len(sequence) < maxlen:
            get_unique_repeats(maxlen, out, sequence)
    return out
#get_unique_repeats


def compute_fit(lengths, observed_amounts, degree, nonnegative):
    if nonnegative:
        """
        The following Least Squares setting is used:
        The xi are the lengths of the i observations, the yi the amounts.
        The a and b are the quadratic and linear factors of the polynomial.
        The c is the shift amount. The example below is for degree=2.
            --            --               --  --
            | x1**2  x1  1 |               | y1 |
            | x2**2  x2  1 |     -- --     | y2 |
            | x3**2  x3  1 |     | a |     | y3 |
            | x4**2  x4  1 |  *  | b |  =  | y4 |
            | x5**2  x5  1 |     | c |     | y5 |
            | x6**2  x6  1 |     -- --     | y6 |
            | x7**2  x7  1 |               | y7 |
            --            --               --  --
        """
        A = np.vander(lengths, degree + 1)
        if np.linalg.matrix_rank(A) == degree + 1:
            return nnls(A, observed_amounts.reshape(-1, 1)).reshape(-1)
    else:
        fit = np.polyfit(lengths, observed_amounts, degree, full=True)
        if fit[2] == degree + 1:
            return fit[0]
    return None
#compute_fit


def test_fit(fit, lengths, observed_amounts):
    p = np.poly1d(fit)
    pp = p.deriv()
    max_x = lengths.max()

    # Find lowest nonnegative x for which y is nonnegative and
    # nondecreasing.
    lower_bound = 0
    while min(p(lower_bound), pp(lower_bound)) < 0 and lower_bound < max_x:
        lower_bound += 1

    # Return R2 = 1 - SSres / SStot and the lower_bound.
    predicted_amounts = np.maximum(p(lengths), 0)  # Nonnegative y.
    predicted_amounts[lengths < lower_bound] = 0  # Nondecreasing low x.
    SSres = np.square(observed_amounts - predicted_amounts).sum()
    SStot = np.square(observed_amounts - observed_amounts.mean()).sum()
    if SStot == 0:
        SStot = 1
        if SSres != 0:
            SSres = 1
    return 1 - SSres / SStot, lower_bound
#test_fit


def print_fit(outfile, fit, lengths, seq, marker, stutter_fold, direction, lower_bound, r2):
    outfile.write("%s\t%s\t%+i\t%i\t%i\t%i\t%s\t%0.3f\t" %
        (seq, marker, stutter_fold, lower_bound, min(lengths), max(lengths), direction, r2))
    outfile.write("\t".join("%.3e" % x for x in reversed(fit.tolist())) + "\n")
#print_fit


def add_sample_data(data, sample_data, sample_alleles, tag, min_pct, min_abs):
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

        if marker not in data["alleles"]:
            data["alleles"][marker] = {}
        try:
            data["alleles"][marker][allele].add(tag)
        except KeyError:
            data["alleles"][marker][allele] = set([tag])

    # Enter the read counts into data and check the thresholds.
    for marker, sequence in sample_data:
        if marker not in sample_alleles:
            # Sample does not participate in this marker.
            continue
        allele = sample_alleles[marker]

        reads = sample_data[marker, sequence]
        factors = [100 / x for x in sample_data[marker, allele]]

        if (tag, marker) not in data["samples"]:
            data["samples"][tag, marker] = {}
        amounts = [count * factor for count, factor in zip(reads, factors)]
        if any(abscount >= min_abs and relcount >= min_pct for abscount, relcount in
               zip(reads, amounts)):
            data["samples"][tag, marker][sequence] = amounts
#add_sample_data


def filter_data(data, min_samples):
    """
    Remove all alleles from data that have less than min_samples
    samples.
    """
    for marker in tuple(data["alleles"]):
        for allele in tuple(data["alleles"][marker]):
            if len(data["alleles"][marker][allele]) < min_samples:
                del data["alleles"][marker][allele]
                continue
        if not data["alleles"][marker]:
            del data["alleles"][marker]
#filter_data


def find_stutter_positions(allele, patterns):
    for seq in sorted(patterns, key=lambda seq: (len(seq), seq)):
        unit_length = len(seq)
        positions = ((m, i == 1)
            for i, pattern in enumerate(patterns[seq])  # One or two.
            for m in iter_overlapping_repeats(allele, pattern, unit_length))
        for m, is_reverse_complement in positions:
            start = m.start()
            end = m.end()
            length = end - start

            # Require at least 3 full repeats, because the pattern
            # returned by get_repeat_pattern() is not guaranteed to
            # report all repeats shorter than 3 repeats minus 2 nt.
            if length >= unit_length * 3:
                yield {
                    "seq": seq,
                    "start": start,
                    "end": end,
                    "length": length,  # revcomp is not yet used, but it should be more efficient to reverse logic and then we can use all.
                    "revcomp": is_reverse_complement}
#find_stutter_positions


def get_all_existing_stutters(data, patterns):
    return {
        marker: {
            allele: get_all_existing_stutters_of_allele(data, patterns, marker, allele, samples)
            for allele, samples in alleles.items()}
        for marker, alleles in data["alleles"].items()}
#get_all_existing_stutters


def get_all_existing_stutters_of_allele(data, patterns, marker, allele, samples):
    stutters1 = sorted(
        get_existing_stutters(data, patterns, marker, allele, samples),
        key=lambda s: (s[1][0]["start"], s[1][0]["end"]))
    all_stutters = dict(stutters1)
    stutters = [s[1][0] for s in stutters1] # TODO: It's a bit icky that I need to unpack it again here

    stuttersn = None
    n = 1
    while True:
        stuttersn = get_possible_larger_combinations(data, patterns, marker, allele, samples, stutters, stuttersn, n)
        if not stuttersn:
            break
        all_stutters.update(stuttersn.values())
        n += 1
    return all_stutters
#get_all_existing_stutters_of_allele


def get_existing_stutters(data, patterns, marker, allele, samples):
    """
    Create a list of (stutter_seq, [stutter]) pairs for all possible
    stutters of the given allele that exist in any of the samples with
    this allele.
    """
    stutters = []
    for stutterpos in find_stutter_positions(allele, patterns):
        unit_length = len(stutterpos["seq"])
        stutter_fold = -1
        while True:
            stutlen = unit_length * abs(stutter_fold)
            success = False
            if stutterpos["length"] >= stutlen:
                stutter = stutterpos.copy()
                stutter["stutlen"] = stutlen
                stutter["fold"] = stutter_fold
                stutter_seq = get_stutter_sequence(allele, [stutter])
                if any(stutter_seq in seqs for sample in samples for seqs in data["samples"][sample, marker]):
                    stutters.append((stutter_seq, [stutter]))
                    success = True
            if success:
                stutter_fold += 1 if stutter_fold > 0 else -1
            elif stutter_fold < 0:
                stutter_fold = 1
            else:
                break
    return stutters
#get_existing_stutters


def get_possible_larger_combinations(data, patterns, marker, allele, samples, stutters, stuttersn, n):
    """
    Create all combinations of n+1 stutters. For n > 1, subject to the
    requirement that all sub-combinations of n stutters exist in
    suttersn (stuttersn may be None for n <= 1).
    The input stutters must be sorted by increasing (start, end).
    """
    stuttersm = {}
    for key in get_all_combinations(stutters, minlen=n+1, maxlen=n+1, indexes_out=True):
        stutter_combination = [stutters[i] for i in key]
        if n > 1 and any(subset not in stuttersn for subset in combinations(key, n)):
            continue
        stutter_seq = get_stutter_sequence(allele, stutter_combination)
        if any(stutter_seq in seqs for sample in samples for seqs in data["samples"][sample, marker]):
            stuttersm[key] = (stutter_seq, stutter_combination)
    return stuttersm
#get_possible_larger_combinations


def fit_stutter(samples_in, outfile, allelefile, min_pct, min_abs,
                min_stutter_pct, min_lengths, min_samples, library, min_r2, orphans, degree,
                same_shape, nonnegative, ignore_zeros, max_unit_length, raw_outfile, only_marker,
                combine_strands, workers):
    # Parse allele list.
    allelelist = {} if allelefile is None else parse_allelelist(
        allelefile, convert="raw", library=library)

    # Read sample data.
    data = {"alleles": {}, "samples": {}}
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
        add_sample_data(data, sample_data, sample_alleles, tag, min_pct, min_abs)

    # Ensure minimum number of samples per allele.
    filter_data(data, min_samples)

    # Compile 1 or 2 regular expressions for each unique repeat sequence.
    patterns = {seq: [get_repeat_pattern(unit) for unit in sorted({seq, reverse_complement(seq)})]
        for seq in get_unique_repeats(max_unit_length)}

    # Calculate all possible stutters from each allele,
    # that actually exist in any samples with that allele.
    possible_stutters = get_all_existing_stutters(data, patterns)

    if raw_outfile != outfile:
        outfile.write("\t".join(
            ["unit", "marker", "stutter", "lbound", "min", "max", "direction", "r2"] +
            list(map(chr, list(range(ord("a"), ord("a") + degree + 1))))) + "\n")
    if raw_outfile is not None:
        raw_outfile.write(
            "\t".join(["unit", "marker", "stutter", "length"] +
                (["total"] if combine_strands else ["forward", "reverse"])) + "\n")

    for seq in sorted(patterns, key=lambda seq: (len(seq), seq)):
        stutter_fold = -1
        while True:
            if fit_stutter_model(outfile, raw_outfile, data, seq, patterns[seq], possible_stutters,
                    min_r2, min_stutter_pct, min_lengths, degree, same_shape, nonnegative, ignore_zeros, stutter_fold, orphans,
                    combine_strands):
                stutter_fold += 1 if stutter_fold > 0 else -1
            elif stutter_fold < 0:
                stutter_fold = 1
            else:
                break
#fit_stutter


def fit_stutter_model(outfile, raw_outfile, data, seq, patterns, possible_stutters, min_r2, min_stutter_pct,
                      min_lengths, degree, same_shape, nonnegative, ignore_zeros, stutter_fold,
                      orphans, combine_strands):
    palindromic = seq == reverse_complement(seq)
    num_amounts = 1 if combine_strands else 2
    num_fits = 1 if palindromic else num_amounts
    stutlen = len(seq) * abs(stutter_fold)

    all_lengths = []
    all_observed_amounts = []
    found_fits = {}
    if same_shape:
        markers = []
        from_markers = []
    for marker in data["alleles"]:
        lengths, observed_amounts = collect_raw_data_marker(
            data, marker, patterns, possible_stutters, palindromic, stutlen, seq, stutter_fold, num_amounts, raw_outfile)
        all_lengths.extend(lengths)
        all_observed_amounts.extend(observed_amounts)
        if lengths:
            if not same_shape:
                found_fits_marker = fit_stutter_model_marker(
                    lengths, observed_amounts, ignore_zeros, min_lengths, min_r2, degree, nonnegative)
                if any(found_fits_marker):
                    found_fits[marker] = found_fits_marker
            else:
                # Prepare for same_shape fit below.
                marker_i = len(markers)
                markers.append(marker)
                from_markers.extend(marker_i for x in lengths)

    all_observed_amounts = np.array(all_observed_amounts)
    if not all_observed_amounts.any() or len(set(all_lengths)) < min_lengths:
        return False

    # Compute fit to all data.
    found_fits_all = fit_stutter_model_marker(
        all_lengths, all_observed_amounts, ignore_zeros, min_lengths, min_r2, degree, nonnegative)
    if any(found_fits_all):
        found_fits["All data"] = found_fits_all

    """
    With same_shape=True, the following Least Squares setting is used:
    The xi are the lengths of the i observations, the yi the amounts.
    The a and b are the quadratic and linear factors of all polynomials.
    The cj are the shift amounts of the different markers j.
    Each '.' is either 0 or 1 depending on which marker the observation
    is from.  The example below is for degree=2.
        --                  --                --  --
        | x1**2  x1  .  .  . |     --  --     | y1 |
        | x2**2  x2  .  .  . |     | a  |     | y2 |
        | x3**2  x3  .  .  . |     | b  |     | y3 |
        | x4**2  x4  .  .  . |  *  | c1 |  =  | y4 |
        | x5**2  x5  .  .  . |     | c2 |     | y5 |
        | x6**2  x6  .  .  . |     | c3 |     | y6 |
        | x7**2  x7  .  .  . |     --  --     | y7 |
        --                  --                --  --
    """

    # Compute same shape fit.
    if same_shape:
        all_lengths = np.array(all_lengths)
        markers = np.array(markers)
        from_markers = np.array(from_markers)
        for i in range(num_fits):
            if not (all_observed_amounts[:, i] > min_stutter_pct).any():
                continue  # All below the threshold.
            if ignore_zeros:
                nonzero_rows = all_observed_amounts[:, i] > 0
                this_lengths = all_lengths[nonzero_rows]
                this_amounts = all_observed_amounts[nonzero_rows, i]
                if len(set(this_lengths)) < min_lengths:
                    continue
                this_markers = []
                this_from_markers = []
                for marker_i in from_markers[nonzero_rows]:
                    if markers[marker_i] not in this_markers:
                        this_markers.append(markers[marker_i])
                    this_from_markers.append(this_markers.index(markers[marker_i]))
                this_from_markers = np.array(this_from_markers)
            else:
                this_lengths = all_lengths
                this_amounts = all_observed_amounts[:, i]
                this_markers = markers
                this_from_markers = from_markers
            A = np.hstack([
                    np.vander(this_lengths, degree + 1)[:, :degree],
                    np.zeros([len(this_lengths), len(this_markers)])])
            for j in range(len(this_from_markers)):
                # TODO I could refactor this loop alike onehot_initialization_v2 from:
                # https://stackoverflow.com/a/36960495
                A[j, this_from_markers[j] + degree] = 1

            if nonnegative:
                if A.shape[1] == np.linalg.matrix_rank(A):
                    fit = nnls(A, this_amounts).reshape(-1)
                else:
                    fit = None
            else:
                fit = np.linalg.lstsq(A, this_amounts, rcond=None) # TODO: Remove rcond=None after dropping Py3.12
                fit = fit[0] if A.shape[1] == fit[2] else None

            if fit is not None:
                for marker_i, marker in enumerate(this_markers):
                    marker_rows = (markers[from_markers] == marker)
                    marker_lengths = this_lengths[this_from_markers == marker_i]
                    marker_f = fit[list(range(degree)) + [marker_i + degree]]
                    r2, lower_bound = test_fit(marker_f, all_lengths[marker_rows],
                        all_observed_amounts[marker_rows, i])
                    if r2 >= min_r2 and lower_bound < max(marker_lengths):
                        if marker not in found_fits:
                            found_fits[marker] = [False] * num_fits
                        found_fits[marker][i] = (marker_f,
                            marker_lengths, lower_bound, r2)

    if raw_outfile != outfile:
        for marker in found_fits:
            if not orphans and not all(found_fits[marker]):
                continue
            for i in range(num_fits):
                if not found_fits[marker][i]:
                    continue
                direction = "total" if combine_strands else "reverse" if i else "forward"
                print_fit(outfile, found_fits[marker][i][0], found_fits[marker][i][1], seq, marker,
                          stutter_fold, direction, found_fits[marker][i][2], found_fits[marker][i][3])

    return any(found_fits)
#fit_stutter_model


def collect_raw_data_marker(data, marker, patterns, possible_stutters, palindromic, stutlen, seq, stutter_fold, num_amounts, raw_outfile):
    lengths = []
    observed_amounts = []
    for allele in data["alleles"][marker]:

        # Get all possible stutter positions in this allele.
        positions = [(m, False) for m in patterns[0].finditer(allele)]
        if not palindromic:
            positions.extend((m, True) for m in patterns[1].finditer(allele))
        for m, is_reverse_complement in positions:
            start = m.start()
            end = m.end()
            length = end - start
            if length <= stutlen:
                continue  # Not repeated.

            for sample in data["alleles"][marker][allele]:
                amount = [0.] * num_amounts  # Reads per 100 reads of allele.
                for stutter_seq, stutter_combination in possible_stutters[marker][allele].items():
                    if not any(s["seq"] == seq and s["fold"] == stutter_fold and s["start"] == start and s["end"] == end for s in stutter_combination):
                        continue  # This stutter combination does not include the stutter we are currently analysing.
                    if stutter_seq in data["samples"][sample, marker]:
                        for i in range(num_amounts):
                            amount[i] += data["samples"][sample, marker][stutter_seq][i]
                if is_reverse_complement:
                    amount = amount[::-1]
                if palindromic and num_amounts == 2:
                    # The forward and reverse reads are now
                    # side-by-side in amounts.  But for palindromic
                    # repeat units, they must be added separately.
                    observed_amounts.append([amount[0]])
                    lengths.append(length)
                    amount = [amount[1]]
                observed_amounts.append(amount)
                lengths.append(length)

    # Write raw data for this marker.
    if raw_outfile is not None:
        for i in range(len(lengths)):
            raw_outfile.write("\t".join(map(str,
                [seq, marker, stutter_fold, lengths[i]] +
                    observed_amounts[i])) + "\n")
    return lengths, observed_amounts
#collect_raw_data_marker


def fit_stutter_model_marker(lengths, observed_amounts, ignore_zeros, min_lengths, min_r2, degree, nonnegative):
    observed_amounts = np.asarray(observed_amounts)
    num_fits = observed_amounts.shape[1]
    found_fits = [False] * num_fits

    # Compute per-marker fit for this marker.
    if len(set(lengths)) >= min_lengths:
        lengths = np.asarray(lengths)
        for i in range(num_fits):
            if not observed_amounts[:, i].any():
                continue  # All zero.
            if ignore_zeros:
                nonzero_rows = observed_amounts[:, i] > 0
                this_lengths = lengths[nonzero_rows]
                this_amounts = observed_amounts[nonzero_rows, i]
                if len(set(this_lengths)) < min_lengths:
                    continue
            else:
                this_lengths = lengths
                this_amounts = observed_amounts[:, i]
            fit = compute_fit(this_lengths, this_amounts, degree, nonnegative)
            if fit is not None:
                r2, lower_bound = test_fit(fit, lengths, observed_amounts[:, i])
                if r2 >= min_r2 and lower_bound < this_lengths.max():
                    found_fits[i] = (fit, this_lengths.tolist(), lower_bound, r2)
    return found_fits
#fit_stutter_model_marker


def add_arguments(parser):
    outputgroup = add_bulk_input_output_args(parser, output_group=True)
    outputgroup.add_argument("-r", "--raw-outfile", metavar="RAWFILE",
        help="write raw data points to this file, for use in stuttermodel "
             "visualisations (specify '-' to write to stdout; normal output on "
             "stdout is then suppressed)")

    parser.add_argument("-a", "--allelelist", metavar="ALLELEFILE",
        help="file containing a list of the true alleles of each sample (e.g., obtained from "
             "allelefinder); if not given, the 'allele' flag from the 'flags' column is used")
    parser.add_argument("-C", "--combine-strands", action="store_true",
        help="if specified, stutter will be modeled for the total number of reads, "
             "instead of separately for either strand")
    parser.add_argument("-T", "--num-threads", metavar="THREADS", type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    parser.add_argument("-D", "--degree", type=pos_int_arg, default=_DEF_DEGREE, metavar="N",
        help="degree of polynomials to fit (default: %(default)s)")
    parser.add_argument("-S", "--same-shape", action="store_true",
        help="if specified, the polynomials of all markers will have equal "
             "coefficients, except for a vertical shift")
    parser.add_argument("-N", "--nonnegative", action="store_true",
        help="if specified, force coefficients to be nonnegative (best used with degree 1 or 2)")
    parser.add_argument("-z", "--ignore-zeros", action="store_true",
        help="if specified, samples exhibiting no stutter are ignored")
    parser.add_argument("-u", "--max-unit-length", type=pos_int_arg,
        default=_DEF_MAX_UNIT_LENGTH, metavar="N",
        help="investigate stutter of repeats of units of up to this number of "
             "nucleotides in length (default: %(default)s)")

    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--min-pct", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to treat as nonzero, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument("-n", "--min-abs", metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="minimum amount of background to treat as nonzero, as an absolute "
             "number of reads (default: %(default)s)")
    filtergroup.add_argument("-p", "--min-stutter-pct", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_STUTTER_PCT,
        help="only fit stutter models where the highest observation has at least this stutter "
             "percentage (default: %4.2f)" % _DEF_THRESHOLD_STUTTER_PCT)
    filtergroup.add_argument("-L", "--min-lengths", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_LENGTHS,
        help="require this minimum number of unique repeat lengths (default: %(default)s)")
    filtergroup.add_argument("-s", "--min-samples", metavar="N",
        type=pos_int_arg, default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele (default: %(default)s)")
    filtergroup.add_argument("-M", "--marker", metavar="MARKER", help="work only on MARKER")
    filtergroup.add_argument("-t", "--min-r2", type=float, default=_DEF_MIN_R2, metavar="N",
        help="minimum required r-squared score (default: %(default)s)")
    filtergroup.add_argument("-O", "--orphans", action="store_true",
        help="if specified, a fit on one strand is reported even if no fit "
             "was obtained on the other strand for the same marker, unit, and stutter depth")

    add_sequence_format_args(parser, default_format="raw", force=True)
#add_arguments


def run(args):
    # Import numpy now.
    global np
    import numpy as np

    samples_in = iterate_input_samples(args.infiles, args.tag_expr, args.tag_format)
    with std_or_open(args.outfile, "wt") as outfile, \
         std_or_open(args.raw_outfile, "wt") as raw_outfile:
        fit_stutter(samples_in, outfile, args.allelelist, args.min_pct,
                    args.min_abs, args.min_stutter_pct, args.min_lengths, args.min_samples,
                    args.library, args.min_r2, args.orphans, args.degree, args.same_shape,
                    args.nonnegative, args.ignore_zeros, args.max_unit_length, raw_outfile,
                    args.marker, args.combine_strands, args.num_threads)
#run
