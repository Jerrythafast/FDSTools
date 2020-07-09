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

import itertools
import re
import sys

from functools import reduce

# Patterns that match entire sequences.
PAT_SEQ_RAW = re.compile("^[ACGT]*$")
PAT_SEQ_IUPAC = re.compile("^[ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]*$")
PAT_SEQ_TSSV = re.compile("^(?:[ACGT]+\(\d+\))*$")
PAT_SEQ_ALLELENAME_STR = re.compile(  # First line: n_ACT[m] or alias.
    "^(?:(?:(?:CE)?-?\d+(?:\.\d+)?_(?:[ACGT]+\[\d+\])*)|((?!_).+?))"
    "(?:_[-+]\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # _+3A>
        "(?!(?P=a))(?:[ACGT]+|-))*$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_SNP = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
        "(?!(?P=a))(?:[ACGT]+|-))+$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_MT = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "(?:-?\d+\.\d+[ACGT]|(?P<a>[ACGT])?\d+(?(a)(?!(?P=a)))(?:[ACGT-]|del)))+$")

# Patterns that match blocks of TSSV-style sequences and allele names.
PAT_TSSV_BLOCK = re.compile("([ACGT]+)\((\d+)\)")
PAT_ALLELENAME_BLOCK = re.compile("([ACGT]+)\[(\d+)\]")

# Patterns that match a single variant.
PAT_VARIANT_STR = re.compile(
    "^(?P<pos>[-+]\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_SNP = re.compile(
    "^(?P<pos>\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_MT = re.compile(
    "^(?P<old>(?P<a>[ACGT])|-?)"
    "(?P<pos>\d+)(?(a)|(?:\.(?P<ins>\d+))?)"
    "(?P<new>[ACGT-]|del)$")

# IUPAC Table of complementary bases.
COMPL = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G", "R": "Y", "Y": "R",
         "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D",
         "a": "t", "t": "a", "u": "a", "g": "c", "c": "g", "r": "y", "y": "r",
         "k": "m", "m": "k", "b": "v", "v": "b", "d": "h", "h": "d"}

# IUPAC Table of ambiguous bases.
IUPAC = {"A": "A",   "C": "C",   "G": "G",   "T": "T",  "U": "T",  "W": "AT",
         "S": "CG",  "M": "AC",  "K": "GT",  "R": "AG", "Y": "CT", "B": "CGT",
         "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
         "a": ("A", ""), "c": ("C", ""), "g": ("G", ""), "t": ("T", ""),
         "u": ("T", ""), "w": ("A", "T", ""), "s": ("C", "G", ""),
         "m": ("A", "C", ""), "k": ("G", "T", ""), "r": ("A", "G", ""),
         "y": ("C", "T", ""), "b": ("C", "G", "T", ""),
         "d": ("A", "G", "T", ""), "h": ("A", "C", "T", ""),
         "v": ("A", "C", "G", ""), "n": ("A", "C", "G", "T", "")}

# Special values that may appear in the place of a sequence.
SEQ_SPECIAL_VALUES = ("No data", "Other sequences")


def get_genome_pos(location, x, *, invert=False):
    """Get the genome position of the x-th base in a sequence."""
    if invert:
        offset = 0
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif pos <= x <= location[i]:
                # x is in the current range
                break
            else:
                offset += location[i] - pos + 1
        else:
            if len(location) % 2:
                raise ValueError("Position %i is outside sequence range" % x)
        return offset + x - pos
    else:
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif location[i]-pos < x:
                # x is after this ending position
                x -= location[i] - pos + 1
            else:
                # x is before this ending position
                break
        return pos + x
#get_genome_pos


def call_variants(template, sequence, *, location=("?", 1), cache=True, debug=False):
    """
    Perform a global alignment of sequence to template and return a
    list of variants detected.  The format (nomenclature) of the
    returned variants depends on the location argument.

    If location is a tuple ("chromosome name", position), with any
    integer for the position, all variants are given as substitutions in
    the form posX>Y.  Insertions and deletions are written as pos.1->Y
    and posX>-, respectively.  The given position is that of the first
    base in the template.  With the location set to "suffix", a plus
    sign is prepended to position numbers and the first base in the
    template is pos=1.  With location set to "prefix", a minus sign is
    prepended and bases are counted from right to left instead.

    If location is a tuple ("M", position) with any integer for the
    position, variants are written following the mtDNA nomenclature
    guidelines.  The given position is that of the first base in the
    template.

    By default, the results of this function are cached.  Set cache to
    False to suppress caching the result and reduce memory usage.

    Setting debug to True will cause the alignment matrices to be
    printed to sys.stdout.  Be aware that they can be quite large.
    """
    # Saving the results in a cache to avoid repeating alignments.
    try:
        return call_variants.cache[template, sequence, location]
    except KeyError:
        cache_key = location

    row_offset = len(template) + 1
    matrix_match = [0] * row_offset * (len(sequence) + 1)
    matrix_gap1 = [-sys.maxsize - 1] * row_offset * (len(sequence) + 1)
    matrix_gap2 = [-sys.maxsize - 1] * row_offset * (len(sequence) + 1)
    matrix_direction = [0] * row_offset * (len(sequence) + 1)

    # Matrix and arrow enum constants.
    M_MATCH = 0
    M_GAP1 = 1
    M_GAP2 = 2
    A_MATCH  = 0b00001
    A_HORZ_O = 0b00010
    A_HORZ_E = 0b00100
    A_VERT_O = 0b01000
    A_VERT_E = 0b10000

    # Settings.
    MATCH_SCORE = 1
    MISMATCH_SCORE = -3
    GAP_OPEN_SCORE = -7
    GAP_EXTEND_SCORE = -2
    variant_format = "%i%s>%s"

    if location == "prefix":
        location = ("prefix", -len(template))
    elif location == "suffix":
        # Include plus signs for position numbers.
        variant_format = "%+i%s>%s"
        location = ("suffix", 1)
    elif not isinstance(location, tuple) or len(location) < 2:
        raise ValueError("Unknown location %r. It should be 'prefix', "
            "'suffix', or a tuple (chromosome, position [, endpos])" % location)
    elif location[0] == "M":
        MATCH_SCORE = 1
        MISMATCH_SCORE = -1
        GAP_OPEN_SCORE = -2
        GAP_EXTEND_SCORE = -1

    for i in range(len(matrix_match)):
        x = i % row_offset
        y = i // row_offset

        # Initialisation of first row and column.
        if x == 0 or y == 0:
            if x != 0:
                # Top row.
                matrix_gap1[i] = GAP_OPEN_SCORE + GAP_EXTEND_SCORE * (x - 1)
                matrix_match[i] = matrix_gap1[i]
                matrix_direction[i] = A_HORZ_E | (A_HORZ_O if x == 1 else 0)
            elif y != 0:
                # Left column.
                matrix_gap2[i] = GAP_OPEN_SCORE + GAP_EXTEND_SCORE * (y - 1)
                matrix_match[i] = matrix_gap2[i]
                matrix_direction[i] = A_VERT_E | (A_VERT_O if y == 1 else 0)
            else:
                # Top left corner.
                matrix_direction[i] = A_MATCH
            continue

        if template[x - 1] == sequence[y - 1]:
            match = MATCH_SCORE
        else:
            match = MISMATCH_SCORE

        options_gap1 = (
            matrix_match[i - 1] + GAP_OPEN_SCORE,
            matrix_gap1[i - 1] + GAP_EXTEND_SCORE)
        matrix_gap1[i] = max(options_gap1)
        if options_gap1[0] > options_gap1[1]:
            matrix_direction[i] |= A_HORZ_O  # Must exit M_GAP1 here.

        options_gap2 = (
            matrix_match[i - row_offset] + GAP_OPEN_SCORE,
            matrix_gap2[i - row_offset] + GAP_EXTEND_SCORE)
        matrix_gap2[i] = max(options_gap2)
        if options_gap2[0] > options_gap2[1]:
            matrix_direction[i] |= A_VERT_O  # Must exit M_GAP2 here.

        options = (
            matrix_match[i - 1 - row_offset] + match,
            matrix_gap1[i],
            matrix_gap2[i])
        matrix_match[i] = max(options)
        if options[0] == matrix_match[i]:
            matrix_direction[i] |= A_MATCH  # Can stay in M_MATCH here.
        if options[1] == matrix_match[i]:
            matrix_direction[i] |= A_HORZ_E  # Can enter M_GAP1 here.
        if options[2] == matrix_match[i]:
            matrix_direction[i] |= A_VERT_E  # Can enter M_GAP2 here.

    if debug:
        print("GAP1")
        for i in range(len(matrix_gap1), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap1[i : i + row_offset]))
        print("GAP2")
        for i in range(len(matrix_gap2), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap2[i : i + row_offset]))
        print("Match")
        for i in range(len(matrix_match), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_match[i : i + row_offset]))
        print("FLAGS")
        for i in range(len(matrix_direction), row_offset):
            print(("%5s|" * row_offset) % tuple("".join((
                "h" if x & A_HORZ_O else " ",
                "H" if x & A_HORZ_E else " ",
                "D" if x & A_MATCH  else " ",
                "V" if x & A_VERT_E else " ",
                "v" if x & A_VERT_O else " "
            )) for x in matrix_direction[i : i + row_offset]))
        print("Traceback")


    # Backtracking.
    variants = []
    variant_template = 0
    variant_sequence = 0
    i = len(matrix_match) - 1
    in_matrix = M_MATCH  # May change before first step.
    while i >= 0:
        x = i % row_offset
        y = i // row_offset
        if debug:
            print("(%i, %i)" % (x, y))

        if in_matrix == M_MATCH:
            # Make gaps as soon as possible (pushed right).
            if matrix_direction[i] & A_HORZ_E:
                in_matrix = M_GAP1
            elif matrix_direction[i] & A_VERT_E:
                in_matrix = M_GAP2
            elif not (matrix_direction[i] & A_MATCH):
                raise ValueError(
                    "Alignment error: Dead route! (This is a bug.) [%s,%s]" % (template, sequence))

        if in_matrix == M_GAP1:
            # Go horizontally.  Deletion.
            variant_template += 1
            if matrix_direction[i] & A_HORZ_O:
                # End of gap, go diagonally after this.
                in_matrix = M_MATCH
            i -= 1
            continue

        if in_matrix == M_GAP2:
            # Go vertically.  Insertion.
            variant_sequence += 1
            if matrix_direction[i] & A_VERT_O:
                # End of gap, go diagonally after this.
                in_matrix = M_MATCH
            i -= row_offset
            continue

        # Go diagonally.  Either match or mismatch.
        if i != 0 and template[x - 1] != sequence[y - 1]:
            # Start/extend mismatch.
            variant_template += 1
            variant_sequence += 1

        else:
            # Match.  Flush variants.
            if variant_template or variant_sequence:
                if location[0] == "M":
                    # MtDNA variants are one-base-at-a-time.
                    for j in range(max(variant_template, variant_sequence)-1, -1, -1):
                        variants.append("%s%i%s%s" % (
                            template[x + j] if j < variant_template else "",#"-",
                            get_genome_pos(location, x + min(j, variant_template - 1)),
                            ".%i" % (j - variant_template + 1) if j >= variant_template else "",
                            sequence[y + j] if j < variant_sequence else "del"))
                elif variant_template == 0:
                    # Insertions: "-131.1->C" instead of "-130->C".
                    variants.append(variant_format % (
                        get_genome_pos(location, x - 1),
                        ".1-",
                        sequence[y : y + variant_sequence]))
                else:
                    variants.append(variant_format % (
                        get_genome_pos(location, x),
                        template[x : x + variant_template],
                        sequence[y : y + variant_sequence] or "-"))
                variant_template = 0
                variant_sequence = 0
        i -= 1 + row_offset

    # Variants were called from right to left.  Reverse their order.
    if location[0] != "prefix":
        variants.reverse()

    # Store the result in the cache.
    if cache:
        call_variants.cache[template, sequence, cache_key] = variants
    return variants
#call_variants
call_variants.cache = {}


def mutate_sequence(seq, variants, *, location=None):
    """Apply the given variants to the given sequence."""
    if not isinstance(location, tuple) or len(location) < 2:
        pattern = PAT_VARIANT_STR
        location = (None, 0)
    elif location[0] == "M":
        pattern = PAT_VARIANT_MT
        location = (location[0], location[1] - 1) + tuple(location[2:])
    else:
        pattern = PAT_VARIANT_SNP
        location = (location[0], location[1] - 1) + tuple(location[2:])

    seq = [[]] + [[base] for base in seq]
    for variant in variants:
        vm = pattern.match(variant)
        if vm is None:
            raise ValueError("Unrecognised variant '%s'" % variant)
        pos = int(vm.group("pos"))
        ins = int(vm.group("ins") or 0)
        old = vm.group("old")
        new = vm.group("new")
        if old == "-":
            old = ""
        if new == "-" or new == "del":
            new = ""
        if pos < 0:
            pos += len(seq)
        pos = get_genome_pos(location, pos, invert=True)
        if pos < 0 or (pos == 0 and not ins) or pos >= len(seq):
            raise ValueError("Position of variant '%s' is outside sequence range" % variant)
        if not ins and old and old != "".join("".join(x[:1]) for x in seq[pos : pos + len(old)]):
            raise ValueError(
                "Incorrect original sequence in variant '%s'; should be '%s'!"
                % (variant, "".join("".join(x[:1]) for x in seq[pos : pos + len(old)])))
        elif not ins and not old:
            # MtDNA substitution with reference base omitted.
            old = "".join("".join(x[:1]) for x in seq[pos : pos + len(new)])
        if not ins:
            # Remove old bases, retaining those inserted between/after.
            seq[pos : pos + len(old)] = [[""] + x[1:] for x in seq[pos : pos + len(old)]]
            # Place new entirely in the position of the first old base.
            seq[pos][0] = new
        else:
            # Insert new exactly ins positions after pos.
            while len(seq[pos]) <= ins:
                seq[pos].append("")
            seq[pos][ins] = new
    return "".join("".join(x) for x in seq)
#mutate_sequence


def reverse_complement(sequence):
    """Return the reverse complement of the given DNA sequence."""
    return "".join(COMPL[x] if x in COMPL else x for x in reversed(sequence))
#reverse_complement


def iupac_expand_ambiguous(seq):
    """Return generator for all possible values of ambiguous seq."""
    return ("".join(x) for x in itertools.product(*((b for b in IUPAC[a]) for a in seq)))
#iupac_expand_ambiguous


def pattern_longest_match(pattern, subject):
    """Return the longest match of the pattern in the subject string."""
    # FIXME, this function tries only one match at each position in the
    # sequence, which is not neccessarily the longest match at that
    # position. For now, we'll search the reverse sequence as well.
    # Re-implement this without regular expressions to test all options.
    reverse = False
    match = None
    pos = 0
    pat = re.compile("".join("(%s){%i,%i}" % x for x in pattern))
    while pos < len(subject):
        m = pat.search(subject, pos)
        if m is None:
            break
        if match is None or m.end() - m.start() > match.end() - match.start():
            match = m
        pos = m.start() + 1

    # Try to find a longer match from the other end.
    if match is not None:
        subject = reverse_complement(subject)
        pos = 0
        pat = re.compile("".join(
            "(%s){%i,%i}" % (reverse_complement(x[0]), x[1], x[2])
            for x in reversed(pattern)))
        while pos < len(subject):
            m = pat.search(subject, pos)
            if m is None:
                break
            if m.end() - m.start() > match.end() - match.start():
                match = m
                reverse = True
            pos = m.start() + 1

    # Extract the blocks from the match.
    match = [] if match is None or not match.group() else reduce(
        lambda x, i: (
            x[0] + [match.group(i)] * ((match.end(i) - x[1]) // len(match.group(i))),
            match.end(i)) if match.group(i) else x,
        range(1, match.lastindex + 1),
        ([], match.start()))[0]

    # Return the match in the same sequence orientation as the input.
    return list(map(reverse_complement, reversed(match))) if reverse else match
#pattern_longest_match


def pattern_longest_match_veryslow(pattern, subject):
    """Return the longest match of the pattern in the subject string."""
    longest = 0
    the_match = []
    # Generate all possible matching sequences for this pattern.
    #print("Finding match of pattern %r to sequence %s" % (pattern, subject))
    for matching_blocks in itertools.product(*(
            [[block[0]] * i for i in range(block[1], block[2] + 1)]
            for block in pattern)):
        matching = itertools.chain.from_iterable(matching_blocks)
        matching_seq = "".join(matching)
        matching_len = len(matching_seq)
        if matching_len <= longest:
            continue
        if matching_seq in subject:
            longest = matching_len
            the_match = matching
    #print("Found match covering %i/%i bases" % (longest, len(subject)))
    return the_match
#pattern_longest_match_veryslow


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


def convert_sequence_tssv_raw(seq):
    """Convert TSSV-style sequence seq to a raw sequence."""
    return "".join(block[0] * int(block[1]) for block in PAT_TSSV_BLOCK.findall(seq))
#convert_sequence_tssv_raw


def convert_sequence_raw_tssv(seq, library, marker, *, return_alias=False):
    """Convert raw sequence seq to shortened TSSV-style notation."""
    # Try to match this marker's pattern, or any of its aliases.
    match = None
    if "aliases" in library:
        for alias in library["aliases"]:
            if library["aliases"][alias]["marker"] == marker:
                match = library["regex"][alias].match(seq)
                if match is not None:
                    marker = alias
                    break
    if match is None and marker in library["regex"]:
        match = library["regex"][marker].match(seq)
    if match is not None:
        parts = ((match.group(i), match.end(i)) for i in range(1, 1 if
            match.lastindex is None else match.lastindex + 1) if match.group(i))

    # Use heuristics if the sequence does not match the pattern.
    else:

        # Find explictily provided prefix and/or suffix if present.
        pre_suf = ["", ""]
        if "prefix" in library and marker in library["prefix"]:
            for prefix in (x for x in library["prefix"][marker]
                    for x in iupac_expand_ambiguous(x)):
                if seq.startswith(prefix):
                    pre_suf[0] = prefix
                    seq = seq[len(prefix):]
                    break
        if "suffix" in library and marker in library["suffix"]:
            for suffix in (x for x in library["suffix"][marker]
                    for x in iupac_expand_ambiguous(x)):
                if seq.endswith(suffix):
                    pre_suf[1] = suffix
                    seq = seq[:-len(suffix)]
                    break

        # Find longest match of middle pattern.
        middle = [(seq, len(pre_suf[0]) + len(seq))] if seq else []
        if middle and marker in library["blocks_middle"]:
            match = pattern_longest_match(library["blocks_middle"][marker], seq)
            matched = "".join(match)
            if matched:

                # If this allele does not match the prefix of this
                # marker, but the canonical prefix of the marker ends
                # with the same sequence as the start of our match, we
                # move that portion of the match into the prefix.
                # Then, we do the same thing with the suffix.
                middle = match
                match_start = seq.index(matched)
                match_end = match_start + len(matched)
                start = match_start
                end = match_end
                modified = False
                if not pre_suf[0] and "prefix" in library and marker in library["prefix"]:
                    ref = library["prefix"][marker][0]
                    if start != len(ref):
                        i = min(len(ref), len(matched))
                        while i > 0:
                            if ref.endswith(matched[:i]):
                                start += i
                                matched = matched[i:]
                                modified = True
                                break
                            i -= 1
                if not pre_suf[1] and "suffix" in library and marker in library["suffix"]:
                    ref = library["suffix"][marker][0]
                    if len(seq) - end != len(ref):
                        i = min(len(ref), len(matched))
                        while i > 0:
                            if ref.startswith(matched[-i:]):
                                end -= i
                                matched = matched[:-i]
                                modified = True
                                break
                            i -= 1
                if modified:
                    from_start = start - match_start
                    from_end = match_end - end
                    while from_start:
                        if from_start < len(middle[0]):
                            middle[0] = middle[0][from_start:]
                            break
                        else:
                            from_start -= len(middle[0])
                            middle = middle[1:]
                    while from_end:
                        if from_end < len(middle[-1]):
                            middle[-1] = middle[-1][:-from_end]
                            break
                        else:
                            from_end -= len(middle[-1])
                            middle = middle[:-1]
                if start:
                    if pre_suf[0]:
                        middle.insert(0, seq[:start])
                        matched = seq[:start] + matched
                    else:
                        pre_suf[0] = seq[:start]
                if end < len(seq):
                    if pre_suf[1]:
                        middle.append(seq[end:])
                        matched = matched + seq[end:]
                    else:
                        pre_suf[1] = seq[end:]
                if middle:
                    middle = reduce(
                        lambda x, y: (x[:-1] if x[-1][0] == y else x) +
                            [(y, x[-1][1] + len(y))], middle[1:],
                            [(middle[0],
                              len(middle[0]) + len(pre_suf[0]))])
                seq = matched

        # Now construct parts.
        parts = []
        if pre_suf[0]:
            parts.append((pre_suf[0], len(pre_suf[0])))
        parts += middle
        if pre_suf[1]:
            parts.append((pre_suf[1], sum(map(len, pre_suf)) + len(seq)))

    seq = reduce(
        lambda a, b: (a[0] + "%s(%i)" % (b[0], (b[1] - a[1]) // len(b[0])), b[1]),
        reduce(
            lambda x, y: x[:-1] + [y] if x[-1][0] == y[0] else x + [y],
            parts,
            [("", 0)]))[0]
    return (seq, marker) if return_alias else seq
#convert_sequence_raw_tssv


def convert_sequence_allelename_tssv(seq, library, marker):
    """Convert allele name seq to TSSV-style sequence notation."""
    # Check whether there is an alias for this sequence.
    alias_of = None
    if "aliases" in library:
        for alias in library["aliases"]:
            if library["aliases"][alias]["marker"] == marker and (
                    seq == library["aliases"][alias]["name"] or
                    seq.startswith(library["aliases"][alias]["name"] + "_")):
                alias_of = marker
                marker = alias
                seq = "".join((
                    "0_", library["aliases"][alias]["sequence"], "[1]",
                    seq[len(library["aliases"][alias]["name"]):]))
                break

    nameformat = None
    if PAT_SEQ_ALLELENAME_MT.match(seq) is not None:
        nameformat = "MtDNA"
    elif PAT_SEQ_ALLELENAME_SNP.match(seq) is not None:
        nameformat = "SNP"
    if nameformat is not None:
        # MtDNA and SNP markers.
        try:
            reference = library["nostr_reference"][marker]
        except KeyError:
            raise ValueError(
                "%s allele '%s' found for marker %s, but "
                "no reference sequence was found in the library" % (nameformat, seq, marker))
        if seq == "REF":
            return reference + "(1)"
        return mutate_sequence(reference, seq.split(),
            location=library["genome_position"].get(marker,
                ("M" if nameformat == "MtDNA" else "", 1))) + "(1)"

    # Note: aliases of mtDNA and SNP markers end up here as well.
    # It should NOT look like an alias now, however.
    match = PAT_SEQ_ALLELENAME_STR.match(seq)
    if match is None or match.group(1) is not None:
        raise ValueError("Invalid allele name '%s' encountered!" % seq)

    allele = seq.split("_")

    # Get and mutate prefix and suffix.
    prefix = ""
    suffix = ""
    if "prefix" in library:
        if marker in library["prefix"]:
            prefix = library["prefix"][marker][0]
        elif alias_of is not None and alias_of in library["prefix"]:
            prefix = library["prefix"][alias_of][0]
    if "suffix" in library:
        if marker in library["suffix"]:
            suffix = library["suffix"][marker][0]
        elif alias_of is not None and alias_of in library["suffix"]:
            suffix = library["suffix"][alias_of][0]
    variants = [[], []]
    for variant in allele[2:]:
        if variant[0] == "-":
            if not prefix:
                raise ValueError("Encountered prefix variant '%s', but marker "
                                 "'%s' has no prefix!" % (variant, marker))
            variants[0].append(variant)
        elif variant[0] == "+":
            if not suffix:
                raise ValueError("Encountered suffix variant '%s', but marker "
                                 "'%s' has no suffix!" % (variant, marker))
            variants[1].append(variant)
        else:
            raise ValueError("Unrecognised variant '%s'" % variant)
    if variants[0]:
        prefix = mutate_sequence(prefix, variants[0])
    if variants[1]:
        suffix = mutate_sequence(suffix, variants[1])

    blocks = []
    if prefix:
        blocks.append((prefix, 1))
    for block in PAT_ALLELENAME_BLOCK.findall(allele[1]):
        blocks.append((block[0], int(block[1])))
    if suffix:
        blocks.append((suffix, 1))
    return "".join("%s(%i)" % block for block in blocks)
#convert_sequence_allelename_tssv


def convert_sequence_raw_allelename(seq, library, marker):
    """Convert raw sequence seq to an allele name."""
    # We actually convert raw->allelename via TSSV format.
    seq, alias = convert_sequence_raw_tssv(seq, library, marker, return_alias=True)
    blocks = PAT_TSSV_BLOCK.findall(seq)

    if "nostr_reference" in library and marker in library["nostr_reference"]:
        # Handle non-STR markers here.
        if alias != marker:
            return library["aliases"][alias]["name"]
        if not blocks:
            # Oh dear, empty sequence... Primer dimer?
            blocks = (("",),)
        if library["nostr_reference"][marker] == blocks[0][0]:
            return "REF"
        return " ".join(
            call_variants(library["nostr_reference"][marker], blocks[0][0],
                location=library["genome_position"].get(marker, ("?", 1))))

    # Find prefix and suffix.
    prefix = suffix = this_prefix = this_suffix = ""
    remaining_blocks = len(blocks)
    if "prefix" in library:
        if alias in library["prefix"]:
            prefix = library["prefix"][alias][0]
        elif marker in library["prefix"]:
            prefix = library["prefix"][marker][0]
        if prefix and remaining_blocks > 0 and blocks[0][1] == "1":
            remaining_blocks -= 1
    if "suffix" in library:
        if alias in library["suffix"]:
            suffix = library["suffix"][alias][0]
        elif marker in library["suffix"]:
            suffix = library["suffix"][marker][0]
        if suffix and remaining_blocks > 0 and blocks[-1][1] == "1":
            remaining_blocks -= 1
    if remaining_blocks > 0 and prefix and blocks[0][1] == "1":
        this_prefix = blocks[0][0]
        blocks = blocks[1:]
    if remaining_blocks > 0 and suffix and blocks[-1][1] == "1":
        this_suffix = blocks[-1][0]
        blocks = blocks[:-1]

    # Generate prefix/suffix variants.
    length = 0
    variants = []
    if prefix != this_prefix:
        variants += call_variants(prefix, this_prefix, location="prefix")
        length += len(this_prefix) - len(prefix)
    if suffix != this_suffix:
        variants += call_variants(suffix, this_suffix, location="suffix")
        length += len(this_suffix) - len(suffix)

    # We are ready to return the allele name of aliases.
    if alias != marker:
        return "_".join([library["aliases"][alias]["name"]] + variants)

    # Compute CE allele number for the other alleles.
    # TODO: perhaps produce a more intelligent name if there is exactly
    #       1 alias with the same length, or only 1 alias sequence is
    #       contained somewhere within the allele.
    blocknames = []
    blocksize = library.get("block_length", {}).get(marker, 4)
    length -= library.get("length_adjust", {}).get(marker, 0)
    for block in blocks:
        blocknames.append("%s[%s]" % (block[0], block[1]))
        length += len(block[0]) * int(block[1])

    allelename = "CE" + str(length // blocksize)
    if length % blocksize:
        allelename += "." + str(length % blocksize)
    return "_".join([allelename, "".join(blocknames)] + variants)
#convert_sequence_raw_allelename


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

    # From TSSV to raw sequence is easy.
    # We'll need a library and marker name for anything else.
    if (library is None or marker is None) and (from_format != "tssv" or to_format != "raw"):
        raise ValueError("Sequence needs to be converted from %s to %s, this "
                          "conversion requires a library file" % (from_format, to_format))

    # Perform conversions.
    if from_format == "allelename":
        seq = convert_sequence_allelename_tssv(seq, library, marker)
    if to_format == "tssv":
        if from_format == "raw":
            return convert_sequence_raw_tssv(seq, library, marker)
        return seq
    if from_format != "raw":
        seq = convert_sequence_tssv_raw(seq)
    if to_format == "raw":
        return seq
    return convert_sequence_raw_allelename(seq, library, marker)
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
