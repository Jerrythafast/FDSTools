#!/usr/bin/env python
"""
Predict background profiles of new alleles based on a model of stutter
occurrence obtained from stuttermodel.

This tool can be used to compute background noise profiles for alleles
for which no reference samples are available.  The profiles are
predicted using a model of stutter occurrence that must have been
created previously using stuttermodel.  A list of sequences should be
given; bgpredict will predict a background noise profile for each of the
provided sequences separately.  The prediction is based completely on
the provided stutter model.

The predicted background noise profiles obtained from bgpredict can be
combined with the output of bgestimate and/or bghomstats using bgmerge.

It is possible to use an entire forensic case sample as the SEQS input
argument of bgpredict to obtain a predicted background noise profile for
each sequence detected in the sample.  When the background noise
profiles thus obtained are combined with those obtained from bgestimate,
bgcorrect may subsequently produce 'cleaner' results if the sample
contained alleles for which no reference samples were available.
"""
import argparse
import sys
#import numpy as np  # Only imported when actually running this tool.

from operator import mul

from ..lib import get_column_ids, reverse_complement, get_repeat_pattern,\
                  mutate_sequence, parse_library,\
                  PAT_SEQ_RAW, ensure_sequence_format, add_sequence_format_args

__version__ = "0.1dev"


# Default values for parameters are specified below.

# Default name of the column that contains the marker name.
# This value can be overridden by the -m command line option.
_DEF_COLNAME_MARKER = "name"

# Default name of the column that contains the allele.
# This value can be overridden by the -a command line option.
_DEF_COLNAME_ALLELE = "allele"

# Default name of the column to write the output to.
# This value can be overridden by the -o command line option.
_DEF_COLNAME_ALLELE_OUT = "allele"

# Default minimum amount of background to consider, as a percentage of
# the highest allele.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT = 0.5

# Default minimum R2 score.
# This value can be overridden by the -t command line option.
_DEF_MIN_R2 = 0.8


def parse_stuttermodel(stuttermodel, min_r2=0, use_all_data=False):
    column_names = stuttermodel.readline().rstrip("\r\n").split("\t")
    (colid_unit, colid_marker, colid_stutter, colid_lbound, colid_direction,
     colid_r2) = get_column_ids(column_names, "unit", "marker", "stutter",
        "lbound", "direction", "r2")
    degree = ord('a')
    colids_coefs = []
    while True:
        try:
            colids_coefs.append(get_column_ids(column_names, chr(degree)))
            degree += 1
        except ValueError:
            break
    degree -= ord('b')
    if degree < 0:
        raise ValueError("Invalid stutter model file: Unable to determine "
                         "polynomial degree!")

    repeat_patterns = {}
    model = {}
    for line in stuttermodel:
        line = line.rstrip("\r\n").split("\t")
        marker = line[colid_marker]
        seq = line[colid_unit]
        stutter_fold = int(line[colid_stutter])
        direction = line[colid_direction]
        lbound = int(line[colid_lbound])
        r2 = float(line[colid_r2])
        if r2 < min_r2:
            continue
        coefs = [float(line[colid_coef])
                 for colid_coef in reversed(colids_coefs)]
        if marker not in model:
            model[marker] = {}
        if not seq or not PAT_SEQ_RAW.match(seq):
            raise ValueError(
                "Invalid stutter model file: Encountered invalid repeat "
                "sequence '%s'!" % seq)
        if direction == "reverse":
            seq = reverse_complement(seq)
        elif direction != "forward":
            raise ValueError(
                "Invalid stutter model file: Unknown sequence strand '%s'!" %
                direction)
        if (seq, stutter_fold) in model[marker]:
            raise ValueError(
                "Invalid stutter model file: Encountered two models for %+i "
                "stutter of %s repeats in marker %s!" %
                (stutter_fold, seq, marker))
        if seq not in repeat_patterns:
            repeat_patterns[seq] = get_repeat_pattern(seq)
        model[marker][seq, stutter_fold] = {
            "lbound": lbound,
            "r2": r2,
            "pat": repeat_patterns[seq],
            "func": lambda x, lbound=lbound,coefs=coefs,degree=degree:
                0. if x < lbound else max(0.,
                    sum(coefs[i] * x**(degree-i) for i in range(len(coefs))))}

    # Extend marker-specific models with "All data" fits where possible.
    if use_all_data and "All data" in model:
        for marker in model:
            if marker == "All data":
                continue
            for key in model["All data"]:
                if key not in model[marker]:
                    model[marker][key] = model["All data"][key]

    return model
#parse_stuttermodel


def get_all_stutters(allele, flanks, model, min_pct):
    """Return a sorted list of all possible stutters."""
    # Include flanks in case the allele starts in a repeat.
    full_allele = flanks[0] + allele + flanks[1]

    stutters = []
    for seq, stutter_fold in model:
        stutlen = len(seq) * abs(stutter_fold)

        # Generate all stutter variants for this sequence.
        for m in model[seq, stutter_fold]["pat"].finditer(full_allele):
            start = m.start()
            end = m.end()
            length = end - start
            if length <= stutlen:
                continue  # Not repeated.
            if stutter_fold > 0:
                if (end < len(flanks[0]) or
                        start > len(full_allele)-len(flanks[1])):
                    continue  # Repeat is embedded in flank.
            elif (len(flanks[0]) > end-stutlen or
                    start+stutlen > len(full_allele)-len(flanks[1])):
                continue  # Shortening repeat disrupts flank.
            amount = model[seq, stutter_fold]["func"](length)
            if amount < min_pct:
                continue
            stutters.append({
                "seq": seq,
                "stutlen": stutlen,
                "start": max(0, start-len(flanks[0])),
                "end":min(len(full_allele)-len(flanks[1]), end)-len(flanks[0]),
                "fold": stutter_fold,
                "amount": amount/100.})
    return sorted(stutters, key=lambda x: (x["start"], x["end"]))
#get_all_stutters


def get_all_combinations(stutters, combinations=None, appendto=None, pos=0,
                         start=0):
    """Return a list of all non-overlapping combinations of stutters."""
    if combinations is None:
        combinations = []
    if appendto is None:
        appendto = []
    for i in range(start, len(stutters)):
        if stutters[i]["start"] < pos:
            continue
        withnewelement = [x for x in appendto]
        withnewelement.append(stutters[i])
        combinations.append(withnewelement)
        get_all_combinations(stutters, combinations, withnewelement,
                             stutters[i]["end"], i+1)
    return combinations
#get_all_combinations


def get_relative_frequencies(stutters, combinations):
    # Compute expected amount of each combination.
    A = np.array([reduce(mul, (s["amount"] for s in combo))
                  for combo in combinations])
    C = np.array([[s in c for c in combinations] for s in stutters])
    S = np.array([stutter["amount"] for stutter in stutters])
    prev_score = cur_score = sys.float_info.max
    for iterations in xrange(200):  # Max 200 iterations should suffice.
        for i in range(len(stutters)):
            A[C[i]] *= stutters[i]["amount"] / A[C[i]].sum()
        prev_score = cur_score
        cur_score = np.square(S - A.dot(C.T)).sum()
        score_change = (prev_score-cur_score)/prev_score
        if abs(cur_score) < 1e-20 or score_change < 1e-10:
            # We have converged (usually within 5 iterations).
            break
    return A.tolist()
#get_relative_frequencies


def predict_profiles(stuttermodel, seqsfile, outfile, marker_column,
                     allele_column, default_marker, use_all_data, crosstab,
                     min_pct, min_r2, seqformat, library):
    # Parse library and stutter model file.
    library = parse_library(library) if library is not None else None
    model = parse_stuttermodel(stuttermodel, min_r2, use_all_data)

    # Read list of sequences.
    seqlist = {}
    column_names = seqsfile.readline().rstrip("\r\n").split("\t")
    colid_allele = get_column_ids(column_names, "allele")
    try:
        colid_name = get_column_ids(column_names, "name")
    except:
        colid_name = None
    for line in seqsfile:
        line = line.rstrip("\r\n").split("\t")
        marker = line[colid_name] if colid_name is not None else default_marker
        if marker not in model:
            if use_all_data and "All data" in model:
                # No marker-specific model available, use "All data".
                model[marker] = model["All data"]
            else:
                continue
        sequence = ensure_sequence_format(line[colid_allele], "raw",
                                          library=library, marker=marker)
        try:
            seqlist[marker].append(sequence)
        except KeyError:
            seqlist[marker] = [sequence]

    if not crosstab:
        outfile.write("\t".join(
            ["marker", "allele", "sequence", "fmean", "rmean"]) + "\n")
    for marker in seqlist:
        p = {
            "true alleles": len(seqlist[marker]),
            "alleles": seqlist[marker],
            "profiles_forward":
                [[100 if i == j else 0
                  for i in range(len(seqlist[marker]))]
                 for j in range(len(seqlist[marker]))],
            "profiles_reverse":
                [[100 if i == j else 0
                  for i in range(len(seqlist[marker]))]
                 for j in range(len(seqlist[marker]))]}
        if library and "flanks" in library and marker in library["flanks"]:
            flanks = library["flanks"][marker]
            flanks_r = map(reverse_complement, reversed(flanks))
        else:
            flanks = ["", ""]
            flanks_r = flanks
        for ai in range(len(seqlist[marker])):
            allele = seqlist[marker][ai]
            for rc in (False, True):
                if rc:
                    allele = reverse_complement(allele)
                stutters = get_all_stutters(allele, flanks_r if rc else flanks,
                                            model[marker], min_pct)
                combinations = get_all_combinations(stutters)
                frequencies = get_relative_frequencies(stutters, combinations)
                for i in range(len(combinations)):
                    freq = frequencies[i] * 100.
                    if freq < min_pct:
                        continue
                    sequence = mutate_sequence(allele, [
                        "%+i.1->%s" %
                            (s["end"], allele[s["end"]-s["stutlen"]:s["end"]])
                        if s["fold"] > 0 else
                            "%+i%s>-" % (s["end"]-s["stutlen"]+1,
                                allele[s["end"]-s["stutlen"]:s["end"]])
                        for s in combinations[i]])
                    if rc:
                        sequence = reverse_complement(sequence)
                    try:
                        si = p["alleles"].index(sequence)
                    except ValueError:
                        p["alleles"].append(sequence)
                        for profile in p["profiles_forward"]:
                            profile.append(0)
                        for profile in p["profiles_reverse"]:
                            profile.append(0)
                        si = -1
                    if rc:
                        p["profiles_reverse"][ai][si] = freq
                    else:
                        p["profiles_forward"][ai][si] = freq

        # Round to 3 digits to get rid of funny rounding effects.
        # This method is not that precise anyway.
        for profile in p["profiles_forward"]:
            for i in range(len(profile)):
                profile[i] = round(profile[i], 3)
        for profile in p["profiles_reverse"]:
            for i in range(len(profile)):
                profile[i] = round(profile[i], 3)

        if crosstab:
            # Cross-tabular output (profiles in rows).
            outfile.write("\t".join([marker, "0"] +
                [ensure_sequence_format(seq, seqformat, library=library,
                    marker=marker) for seq in p["alleles"]]) + "\n")
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
                    outfile.write("\t".join([marker] +
                        [ensure_sequence_format(seq, seqformat,
                            library=library, marker=marker)
                         for seq in (p["alleles"][i], p["alleles"][j])] +
                        map(str, [p["profiles_forward"][i][j],
                                  p["profiles_reverse"][i][j]])) + "\n")
#predict_profiles


def add_arguments(parser):
    parser.add_argument('stuttermodel', metavar="STUT",
        type=argparse.FileType("r"),
        help="file containing a trained stutter model")
    parser.add_argument('seqs', metavar="SEQS",
        type=argparse.FileType("r"),
        help="file containing the sequences for which a profile should be "
             "predicted")
    parser.add_argument('outfile', metavar="OUT", nargs="?", default=sys.stdout,
        type=argparse.FileType("w"),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-m', '--marker-column', metavar="COLNAME",
        default=_DEF_COLNAME_MARKER,
        help="name of the column that contains the marker name "
             "(default: '%(default)s')")
    parser.add_argument('-a', '--allele-column', metavar="COLNAME",
        default=_DEF_COLNAME_ALLELE,
        help="name of the column that contains the allele "
             "(default: '%(default)s')")
    parser.add_argument('-M', '--marker', metavar="MARKER",
        help="assume the specified marker for all sequences")
    parser.add_argument('-A', '--use-all-data', action="store_true",
        help="if specified, the 'All data' model is used to predict stutter "
             "whenever no marker-specific model is available for a certain "
             "repeat unit")
    parser.add_argument('-C', '--cross-tabular', action="store_true",
        help="if specified, a space-efficient cross-tabular output format is "
             "used instead of the default tab-separated columns format")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-n', '--min-pct', metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="minimum amount of background to consider, as a percentage "
             "of the highest allele (default: %4.2f)" % _DEF_THRESHOLD_PCT)
    filtergroup.add_argument('-t', '--min-r2', type=float,
        default=_DEF_MIN_R2, metavar="N",
        help="minimum required r-squared score (default: %(default)s)")
    add_sequence_format_args(parser, "raw", True)  # Force raw seqs.
#add_arguments


def run(args):
    # Import numpy now.
    import numpy as np
    global np

    predict_profiles(args.stuttermodel, args.seqs, args.outfile,
                     args.marker_column, args.allele_column, args.marker,
                     args.use_all_data, args.cross_tabular, args.min_pct,
                     args.min_r2, args.sequence_format, args.library)
#run
