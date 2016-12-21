#!/usr/bin/env python

#
# Copyright (C) 2016 Jerry Hoogenboom
#
# This file is part of FDSTools, data analysis tools for Next
# Generation Sequencing of forensic DNA markers.
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

import re, sys, argparse, random, itertools, textwrap, glob
#import numpy as np  # Imported only when calling nnls()

from ConfigParser import RawConfigParser, MissingSectionHeaderError
from StringIO import StringIO

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
PAT_ALIAS = re.compile("^(?!_).+$")

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

# Patterns that match (parts of) an STR definition.
PAT_STR_DEF = re.compile("^(?:(?:(?<=^)|(?<!^)\s+)[ACGT]+\s+\d+\s+\d+)*$")
PAT_STR_DEF_BLOCK = re.compile("([ACGT]+)\s+(\d+)\s+(\d+)")

# Pattern to split a comma-, semicolon-, or space-separated list.
PAT_SPLIT = re.compile("\s*[,; \t]\s*")
PAT_SPLIT_QUOTED = re.compile(
    r""""((?:\\"|[^"])*)"|'((?:\\'|[^'])*)'|(\S+)""")

# Pattern that matches a chromosome name/number.
PAT_CHROMOSOME = re.compile(
    "^(?:[Cc][Hh][Rr](?:[Oo][Mm])?)?([1-9XYM]|1\d|2[0-2])$")

# Default regular expression to capture sample tags in file names.
# This is the default of the -e command line option.
DEF_TAG_EXPR = "^(.*?)(?:\.[^.]+)?$"

# Default formatting template to write sample tags.
# This is the default of the -f command line option.
DEF_TAG_FORMAT = "\\1"

# Default formatting template to construct output file names for batch
# processing.  \1 and \2 refer to sample tag and tool name.
# This is the default for the -o command line option with batch support.
DEF_OUTFILE_FORMAT = "\\1-\\2.out"

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

# TextWrapper object for formatting help texts in generated INI files.
INI_COMMENT = textwrap.TextWrapper(width=79, initial_indent="; ",
    subsequent_indent="; ", break_on_hyphens=False)


def get_genome_pos(location, x, invert=False):
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
                offset += location[i]-pos+1
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
                x -= location[i]-pos+1
            else:
                # x is before this ending position
                break
        return pos + x
#get_genome_pos


def call_variants(template, sequence, location=("?", 1), cache=True,
                  debug=False):
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
    matrix_match = [0] * row_offset * (len(sequence)+1)
    matrix_gap1 = [-sys.maxint-1] * row_offset * (len(sequence)+1)
    matrix_gap2 = [-sys.maxint-1] * row_offset * (len(sequence)+1)
    matrix_direction = [0] * row_offset * (len(sequence)+1)

    # Matrix and arrow enum constants.
    M_MATCH = 0
    M_GAP1 = 1
    M_GAP2 = 2
    A_MATCH = 1
    A_HORZ_O = 2
    A_HORZ_E = 4
    A_VERT_O = 8
    A_VERT_E = 16

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
    elif type(location) != tuple or len(location) < 2:
        raise ValueError("Unknown location %r. It should be 'prefix', "
            "'suffix', or a tuple (chromosome, position [, endpos])" %
            location)
    elif location[0] == "M":
        MATCH_SCORE = 1
        MISMATCH_SCORE = -1
        GAP_OPEN_SCORE = -2
        GAP_EXTEND_SCORE = -1

    for i in range(len(matrix_match)):
        x = i % row_offset
        y = i / row_offset

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

        if template[x-1] == sequence[y-1]:
            match = MATCH_SCORE
        else:
            match = MISMATCH_SCORE

        options_gap1 = (
            matrix_match[i-1] + GAP_OPEN_SCORE,
            matrix_gap1[i-1] + GAP_EXTEND_SCORE)
        matrix_gap1[i] = max(options_gap1)
        if options_gap1[0] > options_gap1[1]:
            matrix_direction[i] |= A_HORZ_O  # Must exit M_GAP1 here.

        options_gap2 = (
            matrix_match[i-row_offset] + GAP_OPEN_SCORE,
            matrix_gap2[i-row_offset] + GAP_EXTEND_SCORE)
        matrix_gap2[i] = max(options_gap2)
        if options_gap2[0] > options_gap2[1]:
            matrix_direction[i] |= A_VERT_O  # Must exit M_GAP2 here.

        options = (
            matrix_match[i-1-row_offset] + match,
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
        for i in range(0, len(matrix_gap1), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap1[i:i+row_offset]))
        print("GAP2")
        for i in range(0, len(matrix_gap2), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap2[i:i+row_offset]))
        print("Match")
        for i in range(0, len(matrix_match), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_match[i:i+row_offset]))
        print("FLAGS")
        for i in range(0, len(matrix_direction), row_offset):
            print(("%5s|" * row_offset) % tuple("".join([
                "h" if x & A_HORZ_O else " ",
                "H" if x & A_HORZ_E else " ",
                "D" if x & A_MATCH else " ",
                "V" if x & A_VERT_E else " ",
                "v" if x & A_VERT_O else " "
            ]) for x in matrix_direction[i:i+row_offset]))
        print("Traceback")


    # Backtracking.
    variants = []
    variant_template = 0
    variant_sequence = 0
    i = len(matrix_match) - 1
    in_matrix = M_MATCH  # May change before first step.
    while i >= 0:
        x = i % row_offset
        y = i / row_offset
        if debug:
            print("(%i, %i)" % (x,y))

        if in_matrix == M_MATCH:
            # Make gaps as soon as possible (pushed right).
            if matrix_direction[i] & A_HORZ_E:
                in_matrix = M_GAP1
            elif matrix_direction[i] & A_VERT_E:
                in_matrix = M_GAP2
            elif not (matrix_direction[i] & A_MATCH):
                raise ValueError(
                    "Alignment error: Dead route! (This is a bug.) [%s,%s]" %
                    (template,sequence))

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
                    for j in range(
                            max(variant_template, variant_sequence)-1, -1, -1):
                        variants.append("%s%i%s%s" % (
                            template[x+j] if j < variant_template else "",#"-",
                            get_genome_pos(
                                location, x + min(j, variant_template-1)),
                            ".%i" % (j-variant_template+1)
                                if j >= variant_template else "",
                            sequence[y+j] if j < variant_sequence else "del"))
                elif variant_template == 0:
                    # Insertions: "-131.1->C" instead of "-130->C".
                    variants.append(variant_format % (
                        get_genome_pos(location, x - 1),
                        ".1-",
                        sequence[y:y+variant_sequence]))
                else:
                    variants.append(variant_format % (
                        get_genome_pos(location, x),
                        template[x:x+variant_template],
                        sequence[y:y+variant_sequence] or "-"))
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


def mutate_sequence(seq, variants, location=None):
    """Apply the given variants to the given sequence."""
    if type(location) != tuple or len(location) < 2:
        pattern = PAT_VARIANT_STR
        location = (None, 0)
    elif location[0] == "M":
        pattern = PAT_VARIANT_MT
        location = (location[0], location[1]-1) + tuple(location[2:])
    else:
        pattern = PAT_VARIANT_SNP
        location = (location[0], location[1]-1) + tuple(location[2:])

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
        pos = get_genome_pos(location, pos, True)
        if pos < 0 or (pos == 0 and not ins) or pos >= len(seq):
            raise ValueError(
                "Position of variant '%s' is outside sequence range" %
                    (variant))
        if (not ins and old and old != "".join("".join(x[:1])
                for x in seq[pos:pos+len(old)])):
            raise ValueError(
                "Incorrect original sequence in variant '%s'; should be '%s'!"
                % (variant, "".join("".join(x[:1])
                    for x in seq[pos:pos+len(old)])))
        elif not ins and not old:
            # MtDNA substitution with reference base omitted.
            old = "".join("".join(x[:1]) for x in seq[pos:pos+len(new)])
        if not ins:
            # Remove old bases, retaining those inserted between/after.
            seq[pos:pos+len(old)] = [
                [""] + x[1:] for x in seq[pos:pos+len(old)]]
            # Place new entirely in the position of the first old base.
            seq[pos][0] = new
        else:
            # Insert new exactly ins positions after pos.
            while len(seq[pos]) <= ins:
                seq[pos].append("")
            seq[pos][ins] = new
    return "".join("".join(x) for x in seq)
#mutate_sequence


def iupac_expand_ambiguous(seq):
    """Return generator for all possible values of ambiguous seq."""
    return ("".join(x) for x in itertools.product(
        *((b for b in IUPAC[a]) for a in seq)))
#iupac_expand_ambiguous


def iupac_pattern_ambiguous(seq):
    """Return regex pattern that matches the ambiguous seq."""
    return "".join(
        (("[%s]" if len(IUPAC[x.upper()]) > 1 else "%s") +
        ("?" if x.islower() else "")) % IUPAC[x.upper()] for x in seq)
#iupac_expand_ambiguous


def parse_library(libfile, stream=False):
    try:
        if not stream:
            libfile = sys.stdin if libfile == "-" else open(libfile, "r")
        if libfile == sys.stdin:
            # Can't seek on pipes, so read it into a buffer first.
            libfile = StringIO(sys.stdin.read())
        try:
            library = parse_library_ini(libfile)
            if not stream:
                libfile.close()
            return library
        except MissingSectionHeaderError:
            # Not an ini file.
            pass
        libfile.seek(0)
        library = parse_library_tsv(libfile)
        if not stream and libfile != sys.stdin:
            libfile.close()
        return library
    except ValueError as err:
        raise argparse.ArgumentTypeError(err)
#parse_library


def parse_library_tsv(handle):
    """
    Parse a TSSV library file (tab-separated values format).

    The provided file should contain at least four columns: marker name,
    left flanking sequence, right flanking sequence, and STR definition.

    Return a nested dict with top-level keys "flanks" and "regex".
    """
    library = {
      "flanks": {},
      "regex": {},
      "blocks_middle": {}
    }
    for line in handle:
        line = [x.strip() for x in line.rstrip("\r\n").split("\t")]
        if line == [""]:
            continue
        if len(line) < 4:
            raise ValueError(
                "Invalid library file: encountered line with %i columns, "
                "need at least 4" % len(line))
        marker = line[0]
        if PAT_SEQ_RAW.match(line[1]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" %
                             (line[1], marker))
        if PAT_SEQ_RAW.match(line[2]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" %
                             (line[2], marker))
        if PAT_STR_DEF.match(line[3]) is None:
            raise ValueError("STR definition '%s' of marker %s is invalid" %
                             (line[3], marker))
        library["flanks"][marker] = line[1:3]
        library["blocks_middle"][marker] = [
            (block[0], int(block[1]), int(block[2])) for block in
                PAT_STR_DEF_BLOCK.findall(line[3])]
        # NOTE: Libconvert tool expects "(seq){num,num}" blocks ONLY!
        library["regex"][marker] = re.compile(
            "".join(("^", "".join(
                "(%s){%s,%s}" % x for x in PAT_STR_DEF_BLOCK.findall(line[3])),
                "$")))
    return library
#parse_library_tsv


def parse_library_ini(handle):
    library = {
      "flanks": {},
      "prefix": {},
      "suffix": {},
      "regex": {},
      "blocks_middle": {},
      "nostr_reference": {},
      "genome_position": {},
      "length_adjust": {},
      "block_length": {},
      "max_expected_copies": {},
      "expected_length": {},
      "aliases": {}
    }
    markers = set()

    ini = RawConfigParser()
    ini.optionxform = str
    ini.readfp(handle)
    for section in ini.sections():
        for marker in ini.options(section):
            value = ini.get(section, marker)
            section_low = section.lower()
            if section_low == "flanks":
                values = PAT_SPLIT.split(value)
                if len(values) != 2:
                    raise ValueError(
                        "For marker %s, %i flanking sequences were given,"
                        "need exactly 2" % (marker, len(values)))
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Flanking sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["flanks"][marker] = values
                markers.add(marker)
            elif section_low == "prefix":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "A prefix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for i in range(len(values)):
                    if (not i and PAT_SEQ_RAW.match(values[i]) is None or
                            i and PAT_SEQ_IUPAC.match(values[i]) is None):
                        raise ValueError(
                            "Prefix sequence '%s' of marker %s is invalid" %
                            (values[i], marker))
                library["prefix"][marker] = values
                markers.add(marker)
            elif section_low == "suffix":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "A suffix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for i in range(len(values)):
                    if (not i and PAT_SEQ_RAW.match(values[i]) is None or
                            i and PAT_SEQ_IUPAC.match(values[i]) is None):
                        raise ValueError(
                            "Suffix sequence '%s' of marker %s is invalid" %
                            (values[i], marker))
                library["suffix"][marker] = values
                markers.add(marker)
            elif section_low == "genome_position":
                values = PAT_SPLIT.split(value)
                chromosome = PAT_CHROMOSOME.match(values[0])
                if chromosome is None:
                    raise ValueError(
                        "Invalid chromosome '%s' for marker %s." %
                        (values[0], marker))
                pos = [chromosome.group(1)]
                for i in range(1, len(values)):
                    try:
                        pos.append(int(values[i]))
                    except:
                        raise ValueError(
                            "Position '%s' of marker %s is not a valid integer"
                            % (values[i], marker))
                    if not i % 2 and pos[-2] >= pos[-1]:
                        raise ValueError(
                            "End position %i of marker %s must be higher than "
                            "corresponding start position %i" %
                            (pos[-1], marker, pos[-2]))
                if len(values) == 1:
                    pos.append(1)
                library["genome_position"][marker] = tuple(pos)
                markers.add(marker)
            elif section_low == "length_adjust":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Length adjustment '%s' of marker %s is not a valid "
                        "integer" % (value, marker))
                library["length_adjust"][marker] = value
                markers.add(marker)
            elif section_low == "block_length":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Block length '%s' of marker %s is not a valid integer"
                        % (value, marker))
                library["block_length"][marker] = value
                markers.add(marker)
            elif section_low == "max_expected_copies":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Maximum number of expected copies '%s' of marker %s "
                        "is not a valid integer" % (value, marker))
                library["max_expected_copies"][marker] = value
                markers.add(marker)
            elif section_low == "aliases":
                values = PAT_SPLIT.split(value)
                if len(values) != 3:
                    raise ValueError("Alias %s does not have 3 values, but %i"
                                     % (marker, len(values)))
                if PAT_SEQ_RAW.match(values[1]) is None:
                    raise ValueError(
                        "Alias sequence '%s' of alias %s is invalid" %
                        (values[1], marker))
                if PAT_ALIAS.match(values[2]) is None:
                    raise ValueError(
                        "Allele name '%s' of alias %s is invalid" %
                        (values[2], marker))
                library["aliases"][marker] = {
                    "marker": values[0],
                    "sequence": values[1],
                    "name": values[2]
                }
                markers.add(marker)
            elif section_low == "repeat":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "Marker %s was encountered in both [repeat] and "
                        "[no_repeat] sections" % marker)
                if PAT_STR_DEF.match(value) is None:
                    raise ValueError(
                        "STR definition '%s' of marker %s is invalid" %
                        (value, marker))
                library["regex"][marker] = value
                markers.add(marker)
            elif section_low == "no_repeat":
                if marker in library["regex"]:
                    raise ValueError(
                        "Marker %s was encountered in both [repeat] and "
                        "[no_repeat] sections" % marker)
                if marker in library["prefix"] or marker in library["suffix"]:
                    raise ValueError(
                        "A prefix or suffix was defined for non-STR marker %s"
                        % marker)
                if PAT_SEQ_RAW.match(value) is None:
                    raise ValueError(
                        "Reference sequence '%s' of marker %s is invalid" %
                        (value, marker))
                library["nostr_reference"][marker] = value
                markers.add(marker)
            elif section_low == "expected_allele_length":
                values = PAT_SPLIT.split(value)
                try:
                    min_length = int(values[0])
                except:
                    raise ValueError(
                        "Minimum expected allele length '%s' of marker %s "
                        "is not a valid integer" % (values[0], marker))
                if len(values) > 1:
                    try:
                        max_length = int(values[1])
                    except:
                        raise ValueError(
                            "Maximum expected allele length '%s' of marker %s "
                            "is not a valid integer" % (values[1], marker))
                else:
                    max_length = sys.maxint
                library["expected_length"][marker] = (min_length, max_length)
                markers.add(marker)

    # Sanity check: prohibit prefix/suffix for aliases of non-STRs.
    for alias in library["aliases"]:
        if library["aliases"][alias]["marker"] in library["nostr_reference"] \
                and (alias in library["prefix"] or alias in library["suffix"]):
            raise ValueError(
                "A prefix or suffix was defined for alias %s of non-STR "
                "marker %s" % (alias, library["aliases"][alias]["marker"]))

    # Sanity check: end position of marker should reflect ref length.
    for marker in library["genome_position"]:
        if marker not in library["nostr_reference"]:
            continue
        pos = library["genome_position"][marker]
        reflength = len(library["nostr_reference"][marker])
        length = 0
        for i in range(2, len(pos), 2):
            length += pos[i] - pos[i-1] + 1
        if reflength < length or (len(pos) % 2 and reflength != length):
            raise ValueError(
                "Length of reference sequence of marker %s is %i bases, but "
                "genome positions add up to %i bases" %
                (marker, reflength, length))

    # Compile regular expressions.
    for marker in markers:
        parts = []
        blocksm = []
        if marker in library["prefix"]:
            parts += ("(%s){0,1}" % iupac_pattern_ambiguous(x)
                for x in library["prefix"][marker])
        elif (marker in library["aliases"] and
                library["aliases"][marker]["marker"] in library["prefix"]):
            parts += ("(%s){0,1}" % iupac_pattern_ambiguous(x)
                for x in library["prefix"][
                    library["aliases"][marker]["marker"]])
        if marker in library["aliases"]:
            blocksm.append((library["aliases"][marker]["sequence"], 0, 1))
        elif marker in library["regex"]:
            blocksm += [(block[0], int(block[1]), int(block[2])) for block in
                        PAT_STR_DEF_BLOCK.findall(library["regex"][marker])]
        parts += ("(%s){%s,%s}" % x for x in blocksm)
        if marker in library["suffix"]:
            parts += ("(%s){0,1}" % iupac_pattern_ambiguous(x)
                for x in library["suffix"][marker])
        elif (marker in library["aliases"] and
                library["aliases"][marker]["marker"] in library["suffix"]):
            parts += ("(%s){0,1}" % iupac_pattern_ambiguous(x)
                for x in library["suffix"][
                    library["aliases"][marker]["marker"]])
        if parts:
            library["regex"][marker] = re.compile(
                "".join(["^"] + parts + ["$"]))
        if blocksm:
            library["blocks_middle"][marker] = blocksm
    return library
#parse_library_ini


def load_profiles(profilefile, library=None):
    # TODO: To be replaced with load_profiles_new (less memory).
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    (colid_marker, colid_allele, colid_sequence, colid_fmean, colid_rmean,
     colid_tool) = get_column_ids(column_names, "marker", "allele", "sequence",
        "fmean", "rmean", "tool")

    profiles = {}
    for line in profilefile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        marker = line[colid_marker]
        if marker not in profiles:
            profiles[marker] = {
                "m": set(),  # To be replaced by its length below.
                "n": set(),  # To be replaced by its length below.
                "seqs": [],
                "forward": {},  # To be replaced by a list below.
                "reverse": {},  # To be replaced by a list below.
                "tool": {}  # To be replaced by a list below.
                }
        allele = ensure_sequence_format(line[colid_allele], "raw",
            library=library, marker=marker)
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
            library=library, marker=marker)
        if (allele, sequence) in profiles[marker]["forward"]:
            raise ValueError(
                "Invalid background noise profiles file: encountered "
                "multiple values for marker '%s' allele '%s' sequence '%s'" %
                (marker, allele, sequence))
        profiles[marker]["forward"][allele,sequence] = float(line[colid_fmean])
        profiles[marker]["reverse"][allele,sequence] = float(line[colid_rmean])
        profiles[marker]["tool"][allele, sequence] = line[colid_tool]
        profiles[marker]["m"].update((allele, sequence))
        profiles[marker]["n"].add(allele)

    # Check completeness and reorder true alleles.
    for marker in profiles:
        profiles[marker]["seqs"] = list(profiles[marker]["n"]) + \
            list(profiles[marker]["m"]-profiles[marker]["n"])
        profiles[marker]["n"] = len(profiles[marker]["n"])
        profiles[marker]["m"] = len(profiles[marker]["m"])
        newprofiles = {"forward": [], "reverse": []}
        tools = []
        for i in range(profiles[marker]["n"]):
            allele = profiles[marker]["seqs"][i]
            for direction in newprofiles:
                newprofiles[direction].append([0] * profiles[marker]["m"])
            tools.append([""] * profiles[marker]["m"])
            for j in range(profiles[marker]["m"]):
                sequence = profiles[marker]["seqs"][j]
                if (allele, sequence) in profiles[marker]["forward"]:
                    for direction in newprofiles:
                        newprofiles[direction][i][j] = \
                            profiles[marker][direction][allele, sequence]
                    tools[i][j] = profiles[marker]["tool"][allele, sequence]
        profiles[marker]["forward"] = newprofiles["forward"]
        profiles[marker]["reverse"] = newprofiles["reverse"]
        profiles[marker]["tool"] = tools

    return profiles
#load_profiles


def load_profiles_new(profilefile, library=None):
    # TODO, rename this to load_profiles to complete transition.
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    (colid_marker, colid_allele, colid_sequence, colid_fmean, colid_rmean,
     colid_tool) = get_column_ids(column_names, "marker", "allele", "sequence",
        "fmean", "rmean", "tool")

    profiles = {}
    for line in profilefile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        marker = line[colid_marker]
        if marker not in profiles:
            profiles[marker] = {}
        allele = ensure_sequence_format(line[colid_allele], "raw",
            library=library, marker=marker)
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
            library=library, marker=marker)
        if allele not in profiles[marker]:
            profiles[marker][allele] = {}
        elif sequence in profiles[marker][allele]:
            raise ValueError(
                "Invalid background noise profiles file: encountered "
                "multiple values for marker '%s' allele '%s' sequence '%s'" %
                (marker, allele, sequence))
        profiles[marker][allele][sequence] = {
            "forward": float(line[colid_fmean]),
            "reverse": float(line[colid_rmean]),
            "tool": line[colid_tool]}

    return profiles
#load_profiles_new


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
        if match is None or m.end()-m.start() > match.end()-match.start():
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
            if m.end()-m.start() > match.end()-match.start():
                match = m
                reverse = True
            pos = m.start() + 1

    # Extract the blocks from the match.
    match = [] if match is None or not match.group() else reduce(
        lambda x, i: (
            x[0] + [match.group(i)]*((match.end(i)-x[1])/len(match.group(i))),
            match.end(i)) if match.group(i) else x,
        range(1, match.lastindex+1), ([], match.start()))[0]

    # Return the match in the same sequence orientation as the input.
    return map(reverse_complement, reversed(match)) if reverse else match
#pattern_longest_match


def pattern_longest_match_veryslow(pattern, subject):
    """Return the longest match of the pattern in the subject string."""
    longest = 0
    the_match = []
    # Generate all possible matching sequences for this pattern.
    #print("Finding match of pattern %r to sequence %s" % (pattern, subject))
    for matching_blocks in itertools.product(*(
            [[block[0]]*i for i in range(block[1], block[2]+1)]
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
        return 'raw'
    if PAT_SEQ_TSSV.match(seq):
        return 'tssv'
    if PAT_SEQ_ALLELENAME_STR.match(seq) or PAT_SEQ_ALLELENAME_MT.match(seq) \
            or PAT_SEQ_ALLELENAME_SNP.match(seq):
        return 'allelename'
    raise ValueError("Unrecognised sequence format")
#detect_sequence_format


def convert_sequence_tssv_raw(seq):
    return "".join(block[0] * int(block[1])
                   for block in PAT_TSSV_BLOCK.findall(seq))
#convert_sequence_tssv_raw


def convert_sequence_raw_tssv(seq, library, marker, return_alias=False):
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
            match.lastindex is None else match.lastindex+1) if match.group(i))

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
        middle = [(seq, len(pre_suf[0])+len(seq))] if seq else []
        if middle and marker in library["blocks_middle"]:
            match = pattern_longest_match(library["blocks_middle"][marker],seq)
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
                if (not pre_suf[0] and "prefix" in library
                        and marker in library["prefix"]):
                    ref = library["prefix"][marker][0]
                    i = min(len(ref), len(matched))
                    while i > 0:
                        if ref.endswith(matched[:i]):
                            start += i
                            matched = matched[i:]
                            modified = True
                            break
                        i -= 1
                if (not pre_suf[1] and "suffix" in library
                        and marker in library["suffix"]):
                    ref = library["suffix"][marker][0]
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
                if middle:
                    middle = reduce(
                        lambda x, y: (x[:-1] if x[-1][0] == y else x) +
                            [(y, x[-1][1]+len(y))], middle[1:],
                            [(middle[0],
                              start+len(middle[0])+len(pre_suf[0]))])

                pre_suf[0] += seq[:start]
                pre_suf[1] = seq[end:] + pre_suf[1]
                seq = matched

        # Now construct parts.
        parts = []
        if pre_suf[0]:
            parts.append((pre_suf[0], len(pre_suf[0])))
        parts += middle
        if pre_suf[1]:
            parts.append((pre_suf[1], sum(map(len,pre_suf))+len(seq)))

    seq = reduce(
        lambda a, b: (a[0] + "%s(%i)" % (b[0], (b[1]-a[1])/len(b[0])), b[1]),
        reduce(
            lambda x, y: x[:-1] + [y] if x[-1][0] == y[0] else x + [y],
            parts,
            [("", 0)]))[0]
    return (seq, marker) if return_alias else seq
#convert_sequence_raw_tssv


def convert_sequence_allelename_tssv(seq, library, marker):
    # Check whether there is an alias for this sequence.
    alias_of = None
    if "aliases" in library:
        for alias in library["aliases"]:
            if library["aliases"][alias]["marker"] == marker and (
                    seq == library["aliases"][alias]["name"] or
                    seq.startswith(library["aliases"][alias]["name"] + "_")):
                alias_of = marker
                marker = alias
                seq = "".join([
                    "0_",
                    library["aliases"][alias]["sequence"] + "[1]",
                    seq[len(library["aliases"][alias]["name"]):]])
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
                "no reference sequence was found in the library" %
                (nameformat, seq, marker))
        if seq == "REF":
            return reference + "(1)"
        return mutate_sequence(reference, seq.split(),
            library["genome_position"].get(marker,
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
    # We actually convert raw->allelename via TSSV format.
    seq, alias = convert_sequence_raw_tssv(seq, library, marker, True)
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
                library["genome_position"].get(marker, ("?", 1))))

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
        variants += call_variants(prefix, this_prefix, "prefix")
        length += len(this_prefix) - len(prefix)
    if suffix != this_suffix:
        variants += call_variants(suffix, this_suffix, "suffix")
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

    allelename = "CE" + str(length / blocksize)
    if length % blocksize:
        allelename += "." + str(length % blocksize)
    return "_".join([allelename, "".join(blocknames)] + variants)
#convert_sequence_raw_allelename


def ensure_sequence_format(seq, to_format, from_format=None, library=None,
                           marker=None):
    """Convert seq to 'raw', 'tssv', or 'allelename' format."""
    if seq in SEQ_SPECIAL_VALUES:
        # Special case.
        return seq
    known_formats = ("raw", "tssv", "allelename")
    if to_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" %
                         (to_format, known_formats))
    if from_format is None:
        from_format = detect_sequence_format(seq) if seq else "raw"
    elif from_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" %
                         (from_format, known_formats))

    # No conversion needed?
    if to_format == from_format:
        return seq

    # From TSSV to raw sequence is easy.
    # We'll need a library and marker name for anything else.
    if (library is None or marker is None) and (from_format != "tssv" or
            to_format != "raw"):
        raise ValueError("Sequence needs to be converted from %s to %s, this "
                          "conversion requires a library file" %
                          (from_format, to_format))

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


def reverse_complement(sequence):
    """Return the reverse complement of the given DNA sequence."""
    return "".join(COMPL[x] if x in COMPL else x for x in reversed(sequence))
#reverse_complement


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
        B = np.matrix(np.zeros([A.shape[1], C.shape[1]]))
    E = A.T * A
    F = A.T * C
    prev_score = cur_score = sys.float_info.max
    for i in range(max_iter):
        for n in range(B.shape[0]):
            nn = list(range(n)) + list(range(n + 1, B.shape[0]))
            tmp = (F[n, :] - E[n, nn] * B[nn, :]) / E[n, n]
            tmp[np.isnan(tmp)] = 0
            tmp[tmp < 0] = 0
            B[n, :] = tmp
        prev_score = cur_score
        cur_score = np.square(C - A * B).sum()
        score_change = (prev_score-cur_score)/prev_score

        if debug:
            if i:
                print("%4i %15.6f %15.6f %6.2f" % (i, cur_score,
                    prev_score-cur_score, 100*score_change))
            else:
                print("%4i %15.6f" % (i, cur_score))

        if not cur_score or score_change < min_change:
            # We have converged.
            break

    return B
#nnls


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
        return {"n": 1, "min": value, "max": value, "mean": value, "m2": 0.0,
                "variance": 0.0}
    stats["n"] += 1
    delta = value - stats["mean"]
    stats["mean"] += delta / stats["n"]
    stats["m2"] += delta * (value - stats["mean"])
    try:
        stats["variance"] = stats["m2"] / (stats["n"] - 1)
        stats["min"] = min(stats["min"], value)
        stats["max"] = max(stats["max"], value)
    except ZeroDivisionError:
        stats["variance"] = 0
        stats["min"] = value
        stats["max"] = value
    return stats
#adjust_stats


def get_repeat_pattern(seq):
    """Return compiled regular expression that finds repeats of seq."""
    return re.compile("".join(             # For AGAT, one obtains:
        ["(?:" * (len(seq)-1)] +           # (?:(?:(?:
        ["%s)?" % x for x in seq[1:]] +    # G)?A)?T)?
        ["(?:", seq, ")+"] +               # (?AGAT)+
        ["(?:%s" % x for x in seq[:-1]] +  # (?:A(?:G(?:A
        [")?" * (len(seq)-1)]))            # )?)?)?
#get_repeat_pattern


def read_sample_data_file(infile, data, annotation_column=None, seqformat=None,
                          library=None, default_marker=None,
                          drop_special_seq=False, after_correction=False,
                          extra_columns=None):
    """Add data from infile to data dict as [marker, sequence]=reads."""
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_sequence = get_column_ids(column_names, "sequence")
    colid_forward = None
    colid_reverse = None
    numtype = int
    if after_correction:
        colid_forward, colid_reverse = get_column_ids(column_names,
            "forward_corrected", "reverse_corrected",
            optional=(after_correction != "require"))
    if colid_forward is None:
        colid_forward = get_column_ids(column_names, "forward")
    else:
        numtype = float
    if colid_reverse is None:
        colid_reverse = get_column_ids(column_names, "reverse")
    else:
        numtype = float
    if extra_columns is not None:
        extra_colids = {c: i for c, i in
            ((c, get_column_ids(column_names, c, optional=extra_columns[c]))
                for c in extra_columns)
            if i is not None}

    # Get marker name column if it exists.
    colid_marker = get_column_ids(column_names, "marker", optional=True)

    # Also try to get annotation column if we have one.
    if annotation_column is not None:
        try:
            colid_annotation = get_column_ids(column_names, annotation_column)
        except:
            annotation_column = None

    found_alleles = []
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if drop_special_seq and line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        marker = line[colid_marker] if colid_marker is not None \
            else default_marker
        sequence = line[colid_sequence] if seqformat is None \
            else ensure_sequence_format(line[colid_sequence], seqformat,
                                        library=library, marker=marker)
        if (annotation_column is not None and
                line[colid_annotation].startswith("ALLELE")):
            found_alleles.append((marker, sequence))
        data[marker, sequence] = map(numtype,
            (line[colid_forward], line[colid_reverse]))
        if extra_columns is not None:
            data[marker, sequence].append(
                {c: line[extra_colids[c]] for c in extra_colids})

    return found_alleles
#read_sample_data_file


def reduce_read_counts(data, limit_reads):
    sum_reads = 0
    for markerallele in data:
        sum_reads += sum(data[markerallele])
    if sum_reads <= limit_reads:
        return

    remove = sorted(random.sample(xrange(sum_reads), sum_reads - limit_reads))
    removed = 0
    seen = 0
    for markerallele in data:
        for direction in (0, 1):
            seen += data[markerallele][direction]
            while removed < len(remove) and seen > remove[removed]:
                data[markerallele][direction] -= 1
                removed += 1
#reduce_read_counts


def get_sample_data(tags_to_files, callback, allelelist=None,
                    annotation_column=None, seqformat=None, library=None,
                    marker=None, homozygotes=False, limit_reads=sys.maxint,
                    drop_samples=0, drop_special_seq=False,
                    after_correction=False, extra_columns=None):
    if drop_samples:
        sample_tags = tags_to_files.keys()
        for tag in random.sample(xrange(len(sample_tags)),
                                 int(len(sample_tags) * drop_samples)):
            del tags_to_files[sample_tags[tag]]

    for tag in tags_to_files:
        data = {}
        alleles = set()
        for infile in tags_to_files[tag]:
            infile = sys.stdin if infile == "-" else open(infile, "r")
            alleles.update(read_sample_data_file(
                infile, data, annotation_column, seqformat, library, marker,
                drop_special_seq, after_correction, extra_columns))
            if infile != sys.stdin:
                infile.close()
        if limit_reads < sys.maxint:
            reduce_read_counts(data, limit_reads)
        if allelelist is not None:
            if tag not in allelelist:
                allelelist[tag] = {}
            for markerx, allele in alleles:
                if markerx not in allelelist[tag]:
                    allelelist[tag][markerx] = set()
                allelelist[tag][markerx].add(allele)
            if marker:
                if marker in allelelist[tag]:
                    allelelist[tag] = {marker: allelelist[tag][marker]}
                else:
                    allelelist[tag] = {}
            if homozygotes:
                for markerx in allelelist[tag].keys():
                    if len(allelelist[tag][markerx]) > 1:
                        del allelelist[tag][markerx]
        callback(tag, data)
#get_sample_data


def get_column_ids(column_names, *names, **optional):
    """Find all names in column_names and return their indices."""
    result = []
    for name in names:
        try:
            result.append(column_names.index(name))
        except ValueError:
            if "optional" in optional and optional["optional"]:
                result.append(None)
            else:
                raise ValueError("Column not found in input file: %s" % name)
    if len(result) == 1:
        return result[0]
    return tuple(result)
#get_column_ids


def parse_allelelist(allelelist, convert=None, library=None):
    """Read allele list from open file handle."""
    column_names = allelelist.readline().rstrip("\r\n").split("\t")
    colid_sample, colid_marker, colid_allele = get_column_ids(column_names,
        "sample", "marker", "allele")
    alleles = {}
    for line in allelelist:
        line = line.rstrip("\r\n").split("\t")
        sample = line[colid_sample]
        marker = line[colid_marker]
        allele = line[colid_allele]
        if convert is not None:
            allele = ensure_sequence_format(allele, convert, library=library,
                                            marker=marker)
        if sample not in alleles:
            alleles[sample] = {}
        if marker not in alleles[sample]:
            alleles[sample][marker] = set()
        alleles[sample][marker].add(allele)
    return alleles
#parse_allelelist


def pos_int_arg(value):
    """Convert str to int, raise ArgumentTypeError if not positive."""
    if not value.isdigit() or not int(value):
        raise argparse.ArgumentTypeError(
            "invalid positive int value: '%s'" % value)
    return int(value)
#pos_int_arg


def regex_arg(value):
    """Compile value into a regular expression."""
    try:
        return re.compile(value)
    except re.error as err:
        raise argparse.ArgumentTypeError(err)
#regex_arg


def add_allele_detection_args(parser):
    group = parser.add_argument_group("allele detection options")
    group.add_argument('-a', '--allelelist', metavar="ALLELEFILE",
        type=argparse.FileType('r'),
        help="file containing a list of the true alleles of each sample "
             "(e.g., obtained from allelefinder)")
    group.add_argument('-c', '--annotation-column', metavar="COLNAME",
        help="name of a column in the sample files, which contains a value "
             "beginning with 'ALLELE' for the true alleles of the sample")
#add_allele_detection_args


def add_random_subsampling_args(parser):
    group = parser.add_argument_group("random subsampling options (advanced)")
    group.add_argument('-Q', '--limit-reads', metavar="N", type=pos_int_arg,
        default=sys.maxint,
        help="simulate lower sequencing depth by randomly dropping reads down "
             "to this maximum total number of reads for each sample")
    group.add_argument('-x', '--drop-samples', metavar="N", type=float,
        default=0, help="randomly drop this fraction of input samples")
#add_random_subsampling_args


def add_sequence_format_args(parser, default_format=None, force=False,
                             require_library=False):
    group = parser.add_argument_group("sequence format options")
    if force:
        group.set_defaults(sequence_format=default_format)
    else:
        group.add_argument('-F', '--sequence-format', metavar="FORMAT",
            choices=("raw", "tssv", "allelename"),
            default=default_format,
            help="convert sequences to the specified format: one of "
                 "%(choices)s (default: " + (
                 "no conversion" if default_format is None else default_format)
                 + ")")
    if require_library:
        parser.add_argument('library', metavar="LIBRARY", type=parse_library,
            help="library file with marker definitions")
    else:
        group.add_argument('-l', '--library', metavar="LIBRARY",
            type=parse_library,
            help="library file for sequence format conversion")
#add_sequence_format_args


def add_input_output_args(parser, single_in=False, batch_support=False,
                          report_out=False):
    """Add arguments for opening sample files to the given parser."""
    # Input file options group.
    if not single_in:
        # Multiple input files: positionals.
        parser.add_argument('infiles', nargs='*', metavar="FILE",
            default=["-"],
            help="the sample data file(s) to process (default: read from "
                 "stdin)")
    elif not batch_support:
        # Single input file and no batches: single positional.
        parser.add_argument('infile', nargs='?', metavar="IN",
            default="-",
            help="the sample data file to process (default: read from stdin)")
    else:
        # Single input file with batch support: single positional and -i
        # option for batches, which are mutually exclusive.
        mutex = parser.add_argument_group(
                    "input file options").add_mutually_exclusive_group()
        mutex.add_argument('infile', nargs='?', metavar="IN",
            default="-",
            help="single sample data file to process (default: read from "
                 "stdin)")
        mutex.add_argument("-i", "--input", dest="infiles", nargs="+",
            metavar="IN",
            help="multiple sample data files to process (use with "
                 "-o/--output)")

    # Output file options group.
    group = parser.add_argument_group("output file options")
    if batch_support and single_in:
        # Single input file with batch support: single positional and -o
        # option for batches, which are mutually exclusive.
        mutex = group.add_mutually_exclusive_group()
        mutex.add_argument('outfile', nargs='?', metavar="OUT",
            default=sys.stdout,
            help="the file to write the output to (default: write to stdout)")
        mutex.add_argument('-o', '--output', dest="outfiles", nargs="+",
            metavar="OUT",
            help="list of names of output files to match with input files "
                 "specified with -i/--input, or a format string to construct "
                 "file names from sample tags; e.g., the default value is "
                 "'\\1-%s.out', which expands to 'sampletag-%s.out'" %
                    ((parser.prog.rsplit(" ", 1)[-1],)*2))
    elif single_in:
        # Single input file and no batch support: single positional.
        parser.add_argument('outfile', nargs='?', metavar="OUT",
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="the file to write the output to (default: write to stdout)")
    elif batch_support:
        # Multiple input files and batch support: use -o option.
        # (This is multi-in, multi-out).
        group.add_argument('-o', '--output', dest="outfiles", nargs="+",
            metavar="OUT",
            default=[sys.stdout],
            help="a single file name to write all output to (default: write "
                 "to stdout) OR a list of names of output files to match with "
                 "input files OR a format string to construct file names from "
                 "sample tags; e.g., the value '\\1-%s.out' expands to "
                 "'sampletag-%s.out'" % ((parser.prog.rsplit(" ", 1)[-1],)*2))
    else:
        # Multiple input files and no batch support: use -o option.
        group.add_argument('-o', '--output', dest="outfile", metavar="FILE",
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="file to write output to (default: write to stdout)")
    if report_out:
        group.add_argument('-R', '--report', metavar="FILE",
            type=argparse.FileType('w'),
            default=sys.stderr,
            help="file to write a report to (default: write to stderr)")

    # Sample tag parsing options group.
    if not single_in or batch_support:
        group = parser.add_argument_group("sample tag parsing options",
            "for details about REGEX syntax and capturing groups, check "
            "https://docs.python.org/howto/regex")
        group.add_argument('-e', '--tag-expr', metavar="REGEX", type=regex_arg,
            default=DEF_TAG_EXPR,
            help="regular expression that captures (using one or more "
                 "capturing groups) the sample tags from the file names; by "
                 "default, the entire file name except for its extension (if "
                 "any) is captured")
        group.add_argument('-f', '--tag-format', metavar="EXPR",
            default=DEF_TAG_FORMAT,
            help="format of the sample tags produced; a capturing group "
                 "reference like '\\n' refers to the n-th capturing group in "
                 "the regular expression specified with -e/--tag-expr (the "
                 "default of '\\1' simply uses the first capturing group); "
                 "with a single sample, you can enter the sample tag here "
                 "explicitly")
#add_input_output_args


def get_tag(filename, tag_expr, tag_format):
    """Return formatted sample tag from filename using regex."""
    try:
        return tag_expr.search(filename).expand(tag_format)
    except:
        return filename
#get_tag


def glob_path(pathname):
    """Yield filenames matching pathname, or pathname if none match."""
    success = False
    for file in glob.iglob(pathname):
        success = True
        yield file
    if not success:
        yield pathname
#glob_path


def get_input_output_files(args, single=False, batch_support=False):
    if single and not batch_support:
        # One infile, one outfile.  Return 2-tuple (infile, outfile).
        if args.infile == "-" and sys.stdin.isatty():
            return False  # No input specified.
        return args.infile, args.outfile


    if not single and not batch_support:
        # N infiles, one outfile.  Return 2-tuple ({tag: infiles}, out).
        if args.infiles == ["-"] and sys.stdin.isatty():
            return False  # No input specified.

        # Glob args.infiles in case the shell didn't (e.g, on Windows).
        tags_to_files = {}
        for infile in (x for x in args.infiles for x in glob_path(x)):
            tag = get_tag(infile, args.tag_expr, args.tag_format)
            try:
                tags_to_files[tag].append(infile)
            except KeyError:
                tags_to_files[tag] = [infile]
        return tags_to_files, args.outfile


    if single and batch_support:
        # N infiles, N outfiles.  Return generator of (tag, [ins], out).
        # Each yielded tuple should cause a separate run of the tool.

        # Glob args.infiles in case the shell didn't (e.g, on Windows).
        infiles = [x for x in args.infiles for x in glob_path(x)] if "infiles"\
                  in args and args.infiles is not None else [args.infile]
        if infiles == ["-"] and sys.stdin.isatty():
            return False  # No input specified.

        outfiles = args.outfiles if "outfiles" in args \
                   and args.outfiles is not None else [args.outfile]
        if len(outfiles) > 1 and len(outfiles) != len(infiles):
            raise ValueError(
                "Number of input files (%i) is not equal to number of output "
                "files (%i)." % (len(infiles), len(outfiles)))

        tags = [get_tag(infile, args.tag_expr, args.tag_format)
                for infile in infiles]

        if len(outfiles) == 1:
            outfile = sys.stdout if outfiles[0] == "-" else outfiles[0]

            if outfile == sys.stdout and len(set(tags)) == 1:
                # Write output of single sample to stdout.
                return ((tag, infiles, outfile) for tag in set(tags))

            # Write output of each sample to its own outfile.
            if outfile == sys.stdout:
                outfile = DEF_OUTFILE_FORMAT
            return ((tag,
                    [infiles[i] for i in range(len(tags)) if tags[i]==tag],
                    open(outfile.replace("\\1", tag).replace("\\2",
                         args.tool), "w")) for tag in set(tags))

        # Link each output file to each input file.
        # Treating files with the same sample tag as separate samples.
        return ((tags[i], [infiles[i]], open(outfiles[i], 'w'))
                for i in range(len(tags)))

    if not single and batch_support:
        # N infiles, one or N outfiles.
        # If one outfile, return ({tag: [infiles]}, outfile).
        # If N outfiles, return generator of (tag, [infiles], outfile).
        raise NotImplementedError(
            "Multi-input with optional multi-output not supported yet.")
#get_input_output_files


def split_quoted_string(text):
    return reduce(
        lambda x, y: x + ["".join([
            y[0].replace("\\\"", "\""),
            y[1].replace("\\'", "'"),
            y[2]])],
        PAT_SPLIT_QUOTED.findall(text), [])
#split_quoted_string


def print_db(text, debug):
    """Print text if debug is True."""
    if debug:
        print(text)
#print_db
