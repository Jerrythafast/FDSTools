#!/usr/bin/env python
"""
Train a stutter prediction model using homozygous samples.
"""
import argparse
import re
#import numpy as np  # Only imported when actually running this tool.

from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_allele_detection_args, parse_allelelist, parse_library,\
                  get_sample_data, add_sequence_format_args, call_variants,\
                  add_random_subsampling_args, reverse_complement,\
                  get_repeat_pattern

__version__ = "0.1dev"


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
_DEF_MIN_R2 = 0.8


def is_repeated_sequence(sequence):
    """Test whether a sequence consists of one repeated element."""
    for i in range(1, len(sequence)):
        if not len(sequence)%i and sequence == sequence[:i]*(len(sequence)/i):
            return True
    return False
#is_repeated_sequence


def lowest_cyclic_variant(sequence):
    """Return the lowest-sorting cyclic variant of sequence."""
    bestseq = sequence
    seqseq = sequence * 2
    qesqes = reverse_complement(sequence) * 2
    for i in range(len(sequence)):
        if seqseq[i:i+len(sequence)] < bestseq:
            bestseq = seqseq[i:i+len(sequence)]
        if qesqes[i:i+len(sequence)] < bestseq:
            bestseq = qesqes[i:i+len(sequence)]
    return bestseq
#lowest_cyclic_variant


def get_unique_repeats(maxlen, out=None, prefix=""):
    """Return the set of unique repeat sequences of length <= maxlen."""
    if out == None:
        out = set()
    for base in 'ACGT':
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

    # Return R2 = 1 - SSres/SStot and the lower_bound.
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


def print_fit(outfile, fit, lengths, seq, marker, stutter_fold, direction,
              lower_bound, r2):
    if lower_bound < max(lengths):
        outfile.write("%s\t%s\t%+i\t%i\t%i\t%i\t%s\t%0.3f\t" %
            (seq, marker, stutter_fold, lower_bound, min(lengths),
             max(lengths), direction, r2))
        outfile.write("\t".join("%.3e" % x for x in fit.tolist()) + "\n")
#print_fit


def add_sample_data(data, sample_data, sample_alleles, tag, min_pct, min_abs):
    # Check presence of all alleles.
    for marker in sample_alleles:
        allele = sample_alleles[marker]
        if (marker, allele) not in sample_data:
            raise ValueError(
                "Missing allele %s of marker %s in sample %s!" %
                (allele, marker, tag))
        elif 0 in sample_data[marker, allele]:
            raise ValueError(
                "Allele %s of marker %s has 0 reads in sample %s!" %
                (allele, marker, tag))
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
        if (tag, marker) not in data["samples"]:
            data["samples"][tag, marker] = {}
        amounts = [count*factor for count, factor in zip(
            sample_data[marker, sequence],
            (100./x for x in sample_data[marker, sample_alleles[marker]]))]
        if sum(abscount >= min_abs and relcount >= min_pct
               for abscount, relcount in
               zip(sample_data[marker, sequence], amounts)):
            data["samples"][tag, marker][sequence] = amounts
#add_sample_data


def filter_data(data, min_samples):
    """
    Remove all alleles from data that have less than min_samples
    samples.
    """
    for marker in data["alleles"].keys():
        for allele in data["alleles"][marker].keys():
            if len(data["alleles"][marker][allele]) < min_samples:
                del data["alleles"][marker][allele]
                continue
        if not data["alleles"][marker]:
            del data["alleles"][marker]
#filter_data


def fit_stutter(samples_in, outfile, allelefile, annotation_column, min_pct,
                min_abs, min_lengths, min_samples, library, min_r2, degree,
                same_shape, ignore_zeros, max_unit_length, raw_outfile, marker,
                limit_reads, drop_samples):

    # Parse library and allele list.
    library = parse_library(library) if library is not None else None
    allelelist = {} if allelefile is None \
                    else parse_allelelist(allelefile, "raw", library)

    # Read sample data.
    data = {"alleles": {}, "samples": {}}
    get_sample_data(
        samples_in,
        lambda tag, sample_data: add_sample_data(
            data, sample_data,
            {m: allelelist[tag][m].pop() for m in allelelist[tag]},
            tag, min_pct, min_abs),
        allelelist, annotation_column, "raw", library, marker, True,
        limit_reads, drop_samples)

    # Ensure minimum number of samples per allele.
    filter_data(data, min_samples)

    # Compile 2 regular expressions for each unique repeat sequence.
    patterns = {seq: [get_repeat_pattern(seq),
                      get_repeat_pattern(reverse_complement(seq))]
                for seq in get_unique_repeats(max_unit_length)}

    if raw_outfile != outfile:
        outfile.write("\t".join(
            ["unit", "marker", "stutter", "lbound", "min", "max", "direction",
             "r2"] +
            map(chr, list(range(ord("a"), ord("a") + degree + 1)))) + "\n")
    if raw_outfile is not None:
        raw_outfile.write("\t".join(
            ["unit", "marker", "stutter", "length", "forward", "reverse"]) +
            "\n")

    for seq in sorted(patterns, key=lambda seq: (len(seq), seq)):
        stutter_fold = -1
        while True:
            if fit_stutter_model(outfile, raw_outfile, data, library, seq,
                    patterns[seq], min_r2, min_lengths, degree, same_shape,
                    ignore_zeros, stutter_fold):
                stutter_fold += 1 if stutter_fold > 0 else -1
            elif stutter_fold < 0:
                stutter_fold = 1
            else:
                break
#fit_stutter


def fit_stutter_model(outfile, raw_outfile, data, library, seq, patterns,
                      min_r2, min_lengths, degree, same_shape, ignore_zeros,
                      stutter_fold):
    success = False

    stutlen = len(seq) * abs(stutter_fold)
    all_lengths = []
    all_observed_amounts = []
    if same_shape:
        marker_i = -1
        markers = []
        from_markers = []
    for marker in data["alleles"]:
        if library and "flanks" in library and marker in library["flanks"]:
            flanks = library["flanks"][marker]
        else:
            flanks = ["", ""]
        lengths = []
        observed_amounts = []
        for allele in data["alleles"][marker]:
            # Include flanks in case the allele starts in a repeat.
            full_allele = flanks[0] + allele + flanks[1]

            # Get all possible stutter positions in this allele.
            positions = reduce(lambda positions, y:
                positions + map(lambda m: (m.end() - stutlen,
                    m.end() - m.start(), y[0]), y[1]),
                [(False, patterns[0].finditer(full_allele)),
                 (True, patterns[1].finditer(full_allele))],
                [])
            positions=[(m, False) for m in patterns[0].finditer(full_allele)]+\
                      [(m, True) for m in patterns[1].finditer(full_allele)]
            for m, is_reverse_complement in positions:
                start = m.start()
                end = m.end()
                length = end - start
                position = end - stutlen
                if length <= stutlen:
                    continue  # Not repeated.
                if stutter_fold > 0:
                    if (end < len(flanks[0]) or
                            start > len(full_allele)-len(flanks[1])):
                        continue  # Repeat is embedded in flank.
                elif (len(flanks[0]) > position or
                        start+stutlen > len(full_allele)-len(flanks[1])):
                    continue  # Shortening repeat disrupts flank.

                # Go via variants to allow variant combinations.
                # NOTE: Beware variant clashes.  When looking for
                # e.g., "+13AGAT>-" with allele "AGATAGACAGATAGAT",
                # to go from this allele to AGATAGATAGAT could be
                # "+8C>T_+13AGAT>-" but optimal is "+8CAGA>-".  It
                # should be included in the analysis but it is not.
                if stutter_fold > 0:
                    variant = "%+i.1->%s" % (
                        end - len(flanks[0]),
                        full_allele[position:end])
                else:
                    variant = "%+i%s>-" % (
                        position + 1 - len(flanks[0]),
                        full_allele[position:end])
                for sample in data["alleles"][marker][allele]:
                    amount = [0., 0.]  # Reads per 100 reads of allele.
                    for sequence in data["samples"][sample, marker]:
                        if variant in call_variants(allele, sequence):
                            amount[0] += \
                               data["samples"][sample, marker][sequence][0]
                            amount[1] += \
                               data["samples"][sample, marker][sequence][1]
                    if is_reverse_complement:
                        amount = amount[::-1]
                    all_observed_amounts.append(amount)
                    observed_amounts.append(amount)
                    all_lengths.append(length)
                    lengths.append(length)
                    if same_shape:
                        if not markers or markers[-1] != marker:
                            markers.append(marker)
                            marker_i += 1
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
        if raw_outfile == outfile or same_shape:
            continue
        if len(set(lengths)) < min_lengths:
            continue
        lengths = np.array(lengths)
        for i in (0, 1):
            if not observed_amounts[:, i].any():
                # All zero.
                continue
            if ignore_zeros:
                this_lengths = lengths[observed_amounts[:, i] > 0]
                this_amounts = observed_amounts[
                    observed_amounts[:, i] > 0, i]
                if len(set(this_lengths)) < min_lengths:
                    continue
            else:
                this_lengths = lengths
                this_amounts = observed_amounts[:, i]
            fit = compute_fit(this_lengths, this_amounts, degree)
            if fit is not None:
                r2, lower_bound = test_fit(fit, lengths,
                                           observed_amounts[:, i])
                if r2 < min_r2:
                    continue
                print_fit(outfile, fit, this_lengths.tolist(), seq, marker,
                          stutter_fold,
                          "reverse" if i else "forward", lower_bound, r2)
                success = True

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
    if raw_outfile == outfile or not all_observed_amounts.any():
        return success
    if len(set(all_lengths)) < min_lengths:
        return success
    all_lengths = np.array(all_lengths)
    if same_shape:
        markers = np.array(markers)
        from_markers = np.array(from_markers)
    for i in (0, 1):
        if not all_observed_amounts[:, i].any():
            # All zero.
            continue
        if ignore_zeros:
            this_lengths = all_lengths[all_observed_amounts[:, i] > 0]
            this_amounts = all_observed_amounts[
                all_observed_amounts[:, i] > 0, i]
            if len(set(this_lengths)) < min_lengths:
                continue
            if same_shape:
                this_markers = []
                this_from_markers = []
                for marker_i in from_markers[all_observed_amounts[:, i] > 0]:
                    if markers[marker_i] not in this_markers:
                        this_markers.append(markers[marker_i])
                    this_from_markers.append(
                        this_markers.index(markers[marker_i]))
                this_from_markers = np.array(this_from_markers)
        else:
            this_lengths = all_lengths
            this_amounts = all_observed_amounts[:, i]
            if same_shape:
                this_markers = markers
                this_from_markers = from_markers
        if same_shape:
            A = np.hstack([
                    np.vander(this_lengths, degree+1)[:,:degree],
                    np.zeros([len(this_lengths), len(this_markers)])])
            for j in range(len(this_from_markers)):
                A[j, this_from_markers[j]+degree] = 1
            fit = np.linalg.lstsq(A, this_amounts)
            if A.shape[1] == fit[2]:
                for marker_i in range(len(this_markers)):
                    marker_rows = (markers[from_markers] ==
                                   this_markers[marker_i])
                    marker_lengths = this_lengths[this_from_markers==marker_i]
                    marker_f = fit[0][list(range(degree)) + [marker_i+degree]]
                    r2, lower_bound = test_fit(marker_f,
                        all_lengths[marker_rows],
                        all_observed_amounts[marker_rows, i])
                    if r2 < min_r2:
                        continue
                    success = True
                    print_fit(outfile, marker_f, marker_lengths, seq,
                              this_markers[marker_i], stutter_fold,
                              "reverse" if i else "forward", lower_bound, r2)

        # Fit a polynomial to all data at once as well.
        fit = compute_fit(this_lengths, this_amounts, degree)
        if fit is not None:
            r2, lower_bound = test_fit(fit, all_lengths,
                                       all_observed_amounts[:, i])
            if r2 < min_r2:
                continue
            success = True
            print_fit(outfile, fit, this_lengths.tolist(), seq, "All data",
                      stutter_fold,
                      "reverse" if i else "forward", lower_bound, r2)

    return success
#fit_stutter_model


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
    filtergroup.add_argument('-L', '--min-lengths', metavar="N",
        type=pos_int_arg, default=_DEF_MIN_LENGTHS,
        help="require this minimum number of unique repeat lengths "
             "(default: %(default)s)")
    filtergroup.add_argument('-s', '--min-samples', metavar="N",
        type=pos_int_arg, default=_DEF_MIN_SAMPLES,
        help="require this minimum number of samples for each true allele "
             "(default: %(default)s)")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    filtergroup.add_argument('-t', '--min-r2', type=float,
        default=_DEF_MIN_R2, metavar="N",
        help="minimum required r-squared score (default: %(default)s)")


    parser.add_argument('-D', '--degree', type=pos_int_arg,
        default=_DEF_DEGREE, metavar="N",
        help="degree of polynomials to fit (default: %(default)s)")
    parser.add_argument('-S', '--same-shape', action="store_true",
        help="if specified, the polynomials of all markers will have equal "
             "coefficients, except for a vertical shift")
    parser.add_argument('-z', '--ignore-zeros', action="store_true",
        help="if specified, samples exhibiting no stutter are ignored")
    parser.add_argument('-u', '--max-unit-length', type=pos_int_arg,
        default=_DEF_MAX_UNIT_LENGTH, metavar="N",
        help="investigate stutter of repeats of units of up to this number of "
             "nucleotides in length (default: %(default)s)")
    parser.add_argument('-r', '--raw-outfile', type=argparse.FileType('w'),
        metavar="RAWFILE",
        help="write raw data points to this file (specify '-' to write to "
             "stdout; normal output on stdout is then supressed)")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
    add_random_subsampling_args(parser)
#add_arguments


def run(args):
    # Import numpy now.
    import numpy as np
    global np

    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    fit_stutter(files[0], files[1], args.allelelist, args.annotation_column,
                args.min_pct, args.min_abs, args.min_lengths, args.min_samples,
                args.library, args.min_r2, args.degree, args.same_shape,
                args.ignore_zeros, args.max_unit_length, args.raw_outfile,
                args.marker, args.limit_reads, args.drop_samples)
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