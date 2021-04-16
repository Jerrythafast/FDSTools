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

import re

from strnaming import libsequence

# Patterns that match entire sequences.
PAT_SEQ_RAW = re.compile("^[ACGT]*$")
PAT_SEQ_IUPAC = re.compile("^[ACGTUWSMKRYBDHVN]*$")
PAT_SEQ_TSSV = re.compile("^(?:[ACGT]+\(\d+\))*$")
PAT_SEQ_ALLELENAME_STR = re.compile(
    "^(?:CE)?-?\d+(?:\.\d+)?_"  # Second line: ACG[n][qA>C]GT[m]
    "(?:[ACGT]+\[\d+\]|\[(?: ?\d+(?:\.1)?[ACGT-]+>[ACGT-]+)*\])*"
    "(?:_[-+]\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # _+3A>
        "(?!(?P=a))(?:[ACGT]+|-))*$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_SNP = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
        "(?!(?P=a))(?:[ACGT]+|-))+$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_MT = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "(?:-?\d+\.\d+[ACGT]|(?P<a>[ACGT])?\d+(?(a)(?!(?P=a)))(?:[ACGT-]|DEL)))+$")

# Special values that may appear in the place of a sequence.
SEQ_SPECIAL_VALUES = ("No data", "Other sequences")


# The following functions are moved to STRNaming.
reverse_complement = libsequence.reverse_complement
call_variants = libsequence.call_variants
mutate_sequence = libsequence.mutate_sequence
PAT_TSSV_BLOCK = libsequence.PAT_TSSV_BLOCK


def detect_sequence_format(seq):
    """Return format of seq.  One of 'raw', 'tssv', or 'allelename'."""
    if not seq:
        raise ValueError("Empty sequence")
    if seq in SEQ_SPECIAL_VALUES:
        # Special case.
        return False
    if PAT_SEQ_RAW.match(seq):
        return "raw"
    if PAT_SEQ_TSSV.match(seq):
        return "tssv"
    if PAT_SEQ_ALLELENAME_STR.match(seq) or PAT_SEQ_ALLELENAME_MT.match(seq) \
            or PAT_SEQ_ALLELENAME_SNP.match(seq):
        return "allelename"
    raise ValueError("Unrecognised sequence format")
#detect_sequence_format


def ensure_sequence_format(seq, to_format, *, from_format=None, library=None, marker=None):
    """Convert seq to 'raw', 'tssv', or 'allelename' format."""
    if seq in SEQ_SPECIAL_VALUES:
        # Special case.
        return seq
    known_formats = ("raw", "tssv", "allelename")
    if to_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" % (to_format, known_formats))
    if from_format is None:
        from_format = detect_sequence_format(seq) if seq else "raw"
    elif from_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" % (from_format, known_formats))

    # No conversion needed?
    if to_format == from_format:
        return seq

    if library is None or marker is None:
        raise ValueError("Sequence needs to be converted from %s to %s, this "
                          "conversion requires a library file" % (from_format, to_format))

    # Perform conversions.
    reported_range = library.get_range(marker)
    if from_format == "allelename":
        seq = reported_range.from_name(seq)
    elif from_format == "tssv":
        seq = reported_range.from_tssv(seq)
    if to_format == "tssv":
        return reported_range.get_tssv(seq)
    if to_format == "allelename":
        return reported_range.get_name(seq)
    return seq
#ensure_sequence_format


def get_repeat_pattern(seq):
    """Return compiled regular expression that finds repeats of seq."""
    return re.compile("".join(             # For AGAT, one obtains:
        ["(?:" * (len(seq)-1)] +           # (?:(?:(?:
        ["%s)?" % x for x in seq[1:]] +    # G)?A)?T)?
        ["(?:", seq, ")+"] +               # (?AGAT)+
        ["(?:%s" % x for x in seq[:-1]] +  # (?:A(?:G(?:A
        [")?" * (len(seq)-1)]))            # )?)?)?
#get_repeat_pattern
