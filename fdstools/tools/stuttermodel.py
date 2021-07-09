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
Train a stutter prediction model using homozygous reference samples.

The model obtained from this tool can be used by bgpredict to predict
background noise profiles of alleles for which no reference samples are
available.
"""
import argparse
import re
#import numpy as np  # Only imported when actually running this tool.

from concurrent.futures import ProcessPoolExecutor
from errno import EPIPE

from ..lib.cli import add_sequence_format_args, add_input_output_args, add_allele_detection_args,\
                      pos_int_arg, get_input_output_files
from ..lib.io import parse_allelelist, get_sample_data
from ..lib.seq import reverse_complement, get_repeat_pattern, call_variants

__version__ = "1.2.0"


# Default values for parameters are specified below.

# Default minimum amount of background to consider, as a percentage of
# the highest allele.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT = 0

# Default minimum number of reads to consider.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS = 1

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


def compute_fit(lengths, observed_amounts, degree):
    fit = np.polyfit(lengths, observed_amounts, degree, None, True)
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


def get_variants_async(data, max_unit_length, workers):
    pat_minimal_repeat = re.compile(
        r"(.).{,%i}\1" % (max_unit_length-1) if max_unit_length > 1 else r"(.)\1")
    variants = {}
    pool = ProcessPoolExecutor(workers)
    for marker in data["alleles"]:
        for allele in data["alleles"][marker]:
            if pat_minimal_repeat.search(allele) is not None:
                sequences = set(sequence for sample in data["alleles"][marker][allele]
                    for sequence in data["samples"][sample, marker])
                for sequence in sequences:
                    variants[allele, sequence] = pool.submit(call_variants, allele, sequence, cache=False)
    pool.shutdown(wait=False)
    return variants
#get_variants_async


def fit_stutter(samples_in, outfile, allelefile, annotation_column, min_pct, min_abs, min_lengths,
                min_samples, library, min_r2, orphans, degree, same_shape, ignore_zeros,
                max_unit_length, raw_outfile, marker, combine_strands, workers):

    # Parse allele list.
    allelelist = {} if allelefile is None else parse_allelelist(allelefile,
        convert="raw", library=library)

    # Read sample data.
    data = {"alleles": {}, "samples": {}}
    get_sample_data(
        samples_in,
        lambda tag, sample_data: add_sample_data(
            data, sample_data,
            {m: allelelist[tag][m].pop() for m in allelelist[tag]},
            tag, min_pct, min_abs),
        allelelist, annotation_column, "raw", library, marker,
        homozygotes=True, drop_special_seq=True, combine_strands=combine_strands)

    # Ensure minimum number of samples per allele.
    filter_data(data, min_samples)

    # Get started on variant calling in separate worker processes.
    variants = get_variants_async(data, max_unit_length, workers) if workers > 1 else None

    # Compile 2 regular expressions for each unique repeat sequence.
    patterns = {seq: [get_repeat_pattern(seq), get_repeat_pattern(reverse_complement(seq))]
                for seq in get_unique_repeats(max_unit_length)}

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
            if fit_stutter_model(outfile, raw_outfile, data, library, seq, patterns[seq], min_r2,
                    min_lengths, degree, same_shape, ignore_zeros, stutter_fold, orphans,
                    combine_strands, variants):
                stutter_fold += 1 if stutter_fold > 0 else -1
            elif stutter_fold < 0:
                stutter_fold = 1
            else:
                break
#fit_stutter


def fit_stutter_model(outfile, raw_outfile, data, library, seq, patterns, min_r2, min_lengths,
                      degree, same_shape, ignore_zeros, stutter_fold, orphans, combine_strands,
                      variants):
    palindromic = seq == reverse_complement(seq)
    num_amounts = 1 if combine_strands else 2
    num_fits = 1 if palindromic else num_amounts
    success = False
    found_fits = {marker: [False] * num_fits for marker in data["alleles"]}
    found_fits["All data"] = [False] * num_fits

    stutlen = len(seq) * abs(stutter_fold)
    all_lengths = []
    all_observed_amounts = []
    if same_shape:
        marker_i = -1
        markers = []
        from_markers = []
    for marker in data["alleles"]:
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
                position = end - stutlen
                if length <= stutlen:
                    continue  # Not repeated.

                # Go via variants to allow variant combinations.
                # NOTE: Beware variant clashes.  When looking for
                # e.g., "13AGAT>-" with allele "AGATAGACAGATAGAT",
                # to go from this allele to AGATAGATAGAT could be
                # "8C>T_+13AGAT>-" but optimal is "8CAGA>-".  It
                # should be included in the analysis but it is not.
                # FIXME: Adjust variant caller settings to fix this.
                if stutter_fold > 0:
                    variant = "%i.1->%s" % (end, allele[position : end])
                else:
                    variant = "%i%s>-" % (position + 1, allele[position : end])
                for sample in data["alleles"][marker][allele]:
                    amount = [0.] * num_amounts  # Reads per 100 reads of allele.
                    for sequence in data["samples"][sample, marker]:
                        if variants:
                            this_variants = variants[allele, sequence].result()
                        else:
                            this_variants = call_variants(allele, sequence)
                        if variant in this_variants:
                            for i in range(num_amounts):
                                amount[i] += data["samples"][sample, marker][sequence][i]
                    if is_reverse_complement:
                        amount = amount[::-1]
                    if palindromic and not combine_strands:
                        # The forward and reverse reads are now
                        # side-by-side in amounts.  But for palindromic
                        # repeat units, they must be added separately.
                        all_observed_amounts.append([amount[0]])
                        observed_amounts.append([amount[0]])
                        all_lengths.append(length)
                        lengths.append(length)
                        amount = [amount[1]]
                    all_observed_amounts.append(amount)
                    observed_amounts.append(amount)
                    all_lengths.append(length)
                    lengths.append(length)
                    if same_shape:
                        if not markers or markers[-1] != marker:
                            markers.append(marker)
                            marker_i += 1
                        from_markers.append(marker_i)
                        if palindromic and not combine_strands:
                            from_markers.append(marker_i)

        # Write raw data for this marker.
        observed_amounts = np.array(observed_amounts)
        if not observed_amounts.any():
            continue
        if raw_outfile is not None:
            for i in range(len(lengths)):
                raw_outfile.write("\t".join(map(str,
                    [seq, marker, stutter_fold, lengths[i]] +
                        observed_amounts[i].tolist())) + "\n")

        # Compute per-marker fit for this marker.
        if same_shape or len(set(lengths)) < min_lengths:
            continue
        lengths = np.array(lengths)
        for i in range(num_fits):
            if not observed_amounts[:, i].any():
                # All zero.
                continue
            if ignore_zeros:
                this_lengths = lengths[observed_amounts[:, i] > 0]
                this_amounts = observed_amounts[observed_amounts[:, i] > 0, i]
                if len(set(this_lengths)) < min_lengths:
                    continue
            else:
                this_lengths = lengths
                this_amounts = observed_amounts[:, i]
            fit = compute_fit(this_lengths, this_amounts, degree)
            if fit is not None:
                r2, lower_bound = test_fit(fit, lengths, observed_amounts[:, i])
                if r2 < min_r2 or lower_bound >= this_lengths.max():
                    continue
                success = True
                if raw_outfile != outfile:
                    found_fits[marker][i] = (fit, this_lengths.tolist(), lower_bound, r2)

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

    # Compute same shape fit and the fit to all data at once.
    all_observed_amounts = np.array(all_observed_amounts)
    if not all_observed_amounts.any() or len(set(all_lengths)) < min_lengths:
        return success
    all_lengths = np.array(all_lengths)
    if same_shape:
        markers = np.array(markers)
        from_markers = np.array(from_markers)
    for i in range(num_fits):
        if not all_observed_amounts[:, i].any():
            # All zero.
            continue
        if ignore_zeros:
            this_lengths = all_lengths[all_observed_amounts[:, i] > 0]
            this_amounts = all_observed_amounts[all_observed_amounts[:, i] > 0, i]
            if len(set(this_lengths)) < min_lengths:
                continue
            if same_shape:
                this_markers = []
                this_from_markers = []
                for marker_i in from_markers[all_observed_amounts[:, i] > 0]:
                    if markers[marker_i] not in this_markers:
                        this_markers.append(markers[marker_i])
                    this_from_markers.append(this_markers.index(markers[marker_i]))
                this_from_markers = np.array(this_from_markers)
        else:
            this_lengths = all_lengths
            this_amounts = all_observed_amounts[:, i]
            if same_shape:
                this_markers = markers
                this_from_markers = from_markers
        if same_shape:
            A = np.hstack([
                    np.vander(this_lengths, degree + 1)[:, :degree],
                    np.zeros([len(this_lengths), len(this_markers)])])
            for j in range(len(this_from_markers)):
                A[j, this_from_markers[j] + degree] = 1
            fit = np.linalg.lstsq(A, this_amounts)
            if A.shape[1] == fit[2]:
                for marker_i in range(len(this_markers)):
                    marker_rows = (markers[from_markers] == this_markers[marker_i])
                    marker_lengths = this_lengths[this_from_markers == marker_i]
                    marker_f = fit[0][list(range(degree)) + [marker_i + degree]]
                    r2, lower_bound = test_fit(marker_f, all_lengths[marker_rows],
                        all_observed_amounts[marker_rows, i])
                    if r2 < min_r2 or lower_bound >= max(marker_lengths):
                        continue
                    success = True
                    if raw_outfile != outfile:
                        found_fits[this_markers[marker_i]][i] = (marker_f,
                        marker_lengths, lower_bound, r2)

        # Fit a polynomial to all data at once as well.
        fit = compute_fit(this_lengths, this_amounts, degree)
        if fit is not None:
            r2, lower_bound = test_fit(fit, all_lengths, all_observed_amounts[:, i])
            if r2 < min_r2 or lower_bound >= this_lengths.max():
                continue
            success = True
            if raw_outfile != outfile:
                found_fits["All data"][i] = (fit, this_lengths.tolist(), lower_bound, r2)

    for marker in found_fits:
        if not orphans and not all(found_fits[marker]):
            continue
        for i in range(num_fits):
            if not found_fits[marker][i]:
                continue
            direction = "total" if combine_strands else "reverse" if i else "forward"
            print_fit(outfile, found_fits[marker][i][0], found_fits[marker][i][1], seq, marker,
                      stutter_fold, direction, found_fits[marker][i][2], found_fits[marker][i][3])

    return success
#fit_stutter_model


def add_arguments(parser):
    add_input_output_args(parser)
    add_allele_detection_args(parser)
    parser.add_argument("-C", "--combine-strands", action="store_true",
        help="if specified, stutter will be modeled for the total number of reads, "
             "instead of separately for either strand")
    parser.add_argument("-T", "--num-threads", metavar="THREADS", type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--min-pct", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to consider, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument("-n", "--min-abs", metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="minimum amount of background to consider, as an absolute "
             "number of reads (default: %(default)s)")
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
    parser.add_argument("-D", "--degree", type=pos_int_arg, default=_DEF_DEGREE, metavar="N",
        help="degree of polynomials to fit (default: %(default)s)")
    parser.add_argument("-S", "--same-shape", action="store_true",
        help="if specified, the polynomials of all markers will have equal "
             "coefficients, except for a vertical shift")
    parser.add_argument("-z", "--ignore-zeros", action="store_true",
        help="if specified, samples exhibiting no stutter are ignored")
    parser.add_argument("-u", "--max-unit-length", type=pos_int_arg,
        default=_DEF_MAX_UNIT_LENGTH, metavar="N",
        help="investigate stutter of repeats of units of up to this number of "
             "nucleotides in length (default: %(default)s)")
    parser.add_argument("-r", "--raw-outfile", type=argparse.FileType("tw", encoding="UTF-8"),
        metavar="RAWFILE",
        help="write raw data points to this file, for use in stuttermodel "
             "visualisations (specify '-' to write to stdout; normal output on "
             "stdout is then suppressed)")
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
        fit_stutter(files[0], files[1], args.allelelist, args.annotation_column, args.min_pct,
                    args.min_abs, args.min_lengths, args.min_samples, args.library, args.min_r2,
                    args.orphans, args.degree, args.same_shape, args.ignore_zeros,
                    args.max_unit_length, args.raw_outfile, args.marker, args.combine_strands,
                    args.num_threads)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#run
