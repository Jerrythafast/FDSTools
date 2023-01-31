#!/usr/bin/env python3

#
# Copyright (C) 2022 Jerry Hoogenboom
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
PAT_SEQ_ALLELENAME_MH = re.compile(
    "^MH_[ACGTN]*"  # Next line: variants preceded by '_', ref may have N.
    "(?:_\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGTN]+)>"
        "(?!(?P=a))(?:[ACGT]+|-))*$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_MT = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "(?:-?\d+\.\d+[ACGT]|(?P<a>[ACGT])?\d+(?(a)(?!(?P=a)))(?:[ACGT-]|DEL)))+$")
PAT_VARIANT_MH = re.compile("^(\d+)([ACGTN]*N[ACGTN]*)>([ACGT-]+)$")
PAT_SEQ_OR_N = re.compile("[ACGT]+|N")

# Special values that may appear in the place of a sequence.
SEQ_SPECIAL_VALUES = ("No data", "Other sequences")


# The following functions are moved to STRNaming.
reverse_complement = libsequence.reverse_complement
call_variants = libsequence.call_variants
mutate_sequence = libsequence.mutate_sequence
PAT_TSSV_BLOCK = libsequence.PAT_TSSV_BLOCK
PAT_VARIANT_SNP = libsequence.PAT_VARIANT_SNP


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
            or PAT_SEQ_ALLELENAME_SNP.match(seq) or PAT_SEQ_ALLELENAME_MH.match(seq):
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
        if reported_range.has_option("microhaplotype_positions"):
            seq = microhaplotype_to_variants(
                seq, reported_range.get_option("microhaplotype_positions"))
        seq = reported_range.from_name(seq)
    elif from_format == "tssv":
        seq = reported_range.from_tssv(seq)
    if to_format == "tssv":
        return reported_range.get_tssv(seq)
    if to_format == "allelename":
        seq = reported_range.get_name(seq)
        if reported_range.has_option("microhaplotype_positions"):
            seq = name_microhaplotype(seq, reported_range.get_option("microhaplotype_positions"))
    return seq
#ensure_sequence_format


def name_microhaplotype(seq, positions):
    """Convert seq from space-separated variants to a microhap name."""
    microhaplotype = {}
    variants = []
    for variant in seq.split(" "):
        match = PAT_VARIANT_MH.match(variant)
        if match and len(match.group(2)) == len(match.group(3)):
            # For now, disallowing most length-changing variants.
            # Only 'N>-' is allowed (since bases_to can be '-').
            offset = int(match.group(1))
            bases_from = match.group(2)
            bases_to = match.group(3)
            for part in PAT_SEQ_OR_N.finditer(bases_from):
                # The part_from is either [ACGT]+ or a single N.
                pos = part.start() + offset
                part_slice = slice(*part.span())
                part_from = bases_from[part_slice]
                part_to = bases_to[part_slice]
                if part_from == "N":
                    # This is one of the defined microhaplotype positions.
                    microhaplotype[pos] = part_to
                else:
                    # Variants adjacent to the microhaplotype positions.
                    variants.append("%i%s>%s" % (pos, part_from, part_to))
        else:
            variants.append(variant)
    return "_".join(
        ["MH", "".join(microhaplotype.get(position, "N") for position in positions)] + variants)
#name_microhaplotype


def microhaplotype_to_variants(seq, positions):
    """Convert seq from a microhap name to space-separated variants."""
    _, microhaplotype, *variants = seq.split("_")
    if len(positions) != len(microhaplotype):
        raise ValueError("Unable to convert microhaplotype name back to sequence: " + seq)
    microhaplotype = {
        position: base for position, base in zip(sorted(positions), microhaplotype) if base != "N"}

    # Merge consecutive microhaplotype positions.
    # Note: adjacent non-microhaplotype SNPs are currently not merged.
    for position in list(sorted(microhaplotype, reverse=True)):
        if position - 1 in microhaplotype and "-" not in (
                microhaplotype[x] for x in (position, position - 1)):
            microhaplotype[position - 1] += microhaplotype[position]
            microhaplotype[position] = ""
    return " ".join(mh_variant for pos, mh_variant in sorted(
        [(position, "%i%s>%s" % (position, "N" * len(base), base))
            for position, base in microhaplotype.items() if base] +
        [(int(PAT_VARIANT_SNP.match(variant).group("pos")), variant)
            for variant in variants]))
#microhaplotype_to_variants


def get_repeat_pattern(seq):
    """Return compiled regular expression that finds repeats of seq."""
    return re.compile("".join(             # For AGAT, one obtains:
        ["(?:" * (len(seq)-1)] +           # (?:(?:(?:
        ["%s)?" % x for x in seq[1:]] +    # G)?A)?T)?
        ["(?:", seq, ")+"] +               # (?:AGAT)+
        ["(?:%s" % x for x in seq[:-1]] +  # (?:A(?:G(?:A
        [")?" * (len(seq)-1)]))            # )?)?)?
#get_repeat_pattern
