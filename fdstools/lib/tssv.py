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

from ..lib.seq import PAT_SEQ_RAW


def seq_pass_filt(sequence, reads, threshold, explen=None):
    """Return False if the sequence does not meet the criteria."""
    return (reads >= threshold and PAT_SEQ_RAW.match(sequence) is not None and
        (explen is None or explen[0] <= len(sequence) <= explen[1]))
#seq_pass_filt


def filter_sequences(sequences, reported_range, minimum, aggregate_filtered, include_missing):
    """
    Filter out sequences with low read counts or invalid bases.
    If aggregate_filtered is True, aggregate them as "Other sequences".
    If include_missing is True, add "No data" if sequences would
    otherwise be empty.
    Return False if sequences is empty, True otherwise.
    """
    # Filter out sequences with low read counts and invalid bases.
    aggregates = [0, 0]
    expected_length = reported_range.get_option("expected_allele_length")
    for sequence, (forward, reverse, *quals) in tuple(sequences.items()):
        if not seq_pass_filt(sequence, forward + reverse, minimum, expected_length):
            aggregates[0] += forward
            aggregates[1] += reverse
            del sequences[sequence]

    # Add aggregate row if the user requested so.
    if aggregate_filtered and any(aggregates):
        sequences["Other sequences"] = aggregates

    # Check presence of any data.
    if include_missing and not sequences:
        sequences["No data"] = [0, 0]

    return bool(sequences)
#filter_sequences
