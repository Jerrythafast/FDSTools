#!/usr/bin/env python3

#
# Copyright (C) 2020 Jerry Hoogenboom
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

import sys
#import numpy as np  # Only imported when it gets used.


def adjust_stats(value, stats=None):
    """
    Adjust the given stats in place with the given observed value and
    return the adjusted stats as well.  If no stats dict is given,
    create a new stats dict with the following initial values:
    {"n": 1, "min": value, "max": value, "mean": value, "m2": 0.0,
     "variance": 0.0}
    """
    value += 0.0
    if not stats:
        return {"n": 1, "min": value, "max": value, "mean": value, "m2": 0.0, "variance": 0.0}
    stats["n"] += 1
    delta = value - stats["mean"]
    stats["mean"] += delta / stats["n"]
    stats["m2"] += delta * (value - stats["mean"])
    stats["variance"] = stats["m2"] / (stats["n"] - 1)
    stats["min"] = min(stats["min"], value)
    stats["max"] = max(stats["max"], value)
    return stats
#adjust_stats


def nnls(A, C, B=None, max_iter=200, min_change=0.0001, debug=False):
    """
    Solve for B in A * B = C in the least squares sense, s.t. B >= 0.

    Hint: call nnls(B.T, C.T).T to solve for A.

    Algorithm has converged if the sum of squared error has decreased
    by less than a factor of min_change in one iteration.  If debug is
    True, print the sum of squared error to sys.stdout after each
    iteration.

    This code was partially adopted from nimfa.methods.factorization.bd.
    """
    import numpy as np
    if B is None:
        B = np.zeros((A.shape[1], C.shape[1]))
    E = A.T @ A
    F = A.T @ C
    prev_score = cur_score = sys.float_info.max
    for i in range(max_iter):
        for n in range(B.shape[0]):
            if E[n, n] == 0:
                B[n, :] = 0
            else:
                tmp = (F[n, :] - E[None, n, :n] @ B[:n, :] - E[None, n, n+1:] @ B[n+1:, :])/E[n, n]
                tmp[tmp < 0] = 0
                B[n, :] = tmp
        prev_score = cur_score
        cur_score = np.square(C - A @ B).sum()
        score_change = (prev_score - cur_score) / prev_score

        if debug:
            if i:
                print("%4i %15.6f %15.6f %6.2f" % (i, cur_score,
                    prev_score - cur_score, 100 * score_change))
            else:
                print("%4i %15.6f" % (i, cur_score))

        if not cur_score or score_change < min_change:
            # We have converged.
            break

    return B
#nnls
