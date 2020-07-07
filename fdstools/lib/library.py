#!/usr/bin/env python3

#
# Copyright (C) 2020 Jerry Hoogenboom
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

import io
import re
import sys
import textwrap

from configparser import RawConfigParser, MissingSectionHeaderError

from .seq import PAT_SEQ_RAW, PAT_SEQ_IUPAC, IUPAC

# Patterns that match (parts of) an STR definition.
PAT_STR_DEF = re.compile("^(?:(?:(?<=^)|(?<!^)\s+)[ACGT]+\s+\d+\s+\d+)*$")
PAT_STR_DEF_BLOCK = re.compile("([ACGT]+)\s+(\d+)\s+(\d+)")

# Pattern that matches a valid sequence alias name.
PAT_ALIAS = re.compile("^(?!_).+$")

# Pattern that matches a chromosome name/number.
PAT_CHROMOSOME = re.compile("^(?:[Cc][Hh][Rr](?:[Oo][Mm])?)?([1-9XYM]|1\d|2[0-2])$")

# Pattern to split a comma-, semicolon-, or space-separated list.
PAT_SPLIT = re.compile("\s*[,; \t]\s*")

# TextWrapper object for formatting help texts in generated INI files.
INI_COMMENT = textwrap.TextWrapper(width=79, initial_indent="; ",
    subsequent_indent="; ", break_on_hyphens=False)


def iupac_pattern_ambiguous(seq):
    """Return regex pattern that matches the ambiguous seq."""
    return "".join(
        (("[%s]" if len(IUPAC[x.upper()]) > 1 else "%s") +
        ("?" if x.islower() else "")) % IUPAC[x.upper()] for x in seq)
#iupac_pattern_ambiguous


def parse_library(libfile, *, stream=False):
    """Parse the given libfile or stream and return a library."""
    if not stream:
        libfile = sys.stdin if libfile == "-" else open(libfile, "tr")
    if libfile == sys.stdin or not libfile.seekable():
        # Can't seek on pipes, so read it into a buffer first.
        # On Windows, sys.stdin may falsely report seekable=True.
        contents = libfile.read()
        if not stream and libfile != sys.stdin:
            libfile.close()
        libfile = io.StringIO(contents)
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
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" % (line[1], marker))
        if PAT_SEQ_RAW.match(line[2]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" % (line[2], marker))
        if PAT_STR_DEF.match(line[3]) is None:
            raise ValueError("STR definition '%s' of marker %s is invalid" % (line[3], marker))
        library["flanks"][marker] = line[1 : 3]
        library["blocks_middle"][marker] = [
            (block[0], int(block[1]), int(block[2])) for block in
                PAT_STR_DEF_BLOCK.findall(line[3])]
        # NOTE: Libconvert tool expects "(seq){num,num}" blocks ONLY!
        library["regex"][marker] = re.compile("".join((
            "^", "".join("(%s){%s,%s}" % x for x in PAT_STR_DEF_BLOCK.findall(line[3])), "$")))
    return library
#parse_library_tsv


def parse_library_ini(handle):
    """Parse an FDSTools library file as made with the library tool."""
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
                        "For marker %s, %i flanking sequences were given, "
                        "need exactly 2" % (marker, len(values)))
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Flanking sequence '%s' of marker %s is invalid" % (value, marker))
                library["flanks"][marker] = values
                markers.add(marker)
            elif section_low == "prefix":
                if marker in library["nostr_reference"]:
                    raise ValueError("A prefix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for i in range(len(values)):
                    if (not i and PAT_SEQ_RAW.match(values[i]) is None or
                            i and PAT_SEQ_IUPAC.match(values[i]) is None):
                        raise ValueError(
                            "Prefix sequence '%s' of marker %s is invalid" % (values[i], marker))
                library["prefix"][marker] = values
                markers.add(marker)
            elif section_low == "suffix":
                if marker in library["nostr_reference"]:
                    raise ValueError("A suffix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for i in range(len(values)):
                    if (not i and PAT_SEQ_RAW.match(values[i]) is None or
                            i and PAT_SEQ_IUPAC.match(values[i]) is None):
                        raise ValueError(
                            "Suffix sequence '%s' of marker %s is invalid" % (values[i], marker))
                library["suffix"][marker] = values
                markers.add(marker)
            elif section_low == "genome_position":
                values = PAT_SPLIT.split(value)
                chromosome = PAT_CHROMOSOME.match(values[0])
                if chromosome is None:
                    raise ValueError("Invalid chromosome '%s' for marker %s" % (values[0], marker))
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
                            "corresponding start position %i" % (pos[-1], marker, pos[-2]))
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
                        "Block length '%s' of marker %s is not a valid integer" % (value, marker))
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
                        "Alias sequence '%s' of alias %s is invalid" % (values[1], marker))
                if PAT_ALIAS.match(values[2]) is None:
                    raise ValueError(
                        "Allele name '%s' of alias %s is invalid" % (values[2], marker))
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
                        "STR definition '%s' of marker %s is invalid" % (value, marker))
                library["regex"][marker] = value
                markers.add(marker)
            elif section_low == "no_repeat":
                if marker in library["regex"]:
                    raise ValueError(
                        "Marker %s was encountered in both [repeat] and "
                        "[no_repeat] sections" % marker)
                if marker in library["prefix"] or marker in library["suffix"]:
                    raise ValueError(
                        "A prefix or suffix was defined for non-STR marker %s" % marker)
                if PAT_SEQ_RAW.match(value) is None:
                    raise ValueError(
                        "Reference sequence '%s' of marker %s is invalid" % (value, marker))
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
                    max_length = sys.maxsize
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
            length += pos[i] - pos[i - 1] + 1
        if reflength < length or (len(pos) % 2 and reflength != length):
            raise ValueError(
                "Length of reference sequence of marker %s is %i bases, but "
                "genome positions add up to %i bases" % (marker, reflength, length))

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
            library["regex"][marker] = re.compile("".join(["^"] + parts + ["$"]))
        if blocksm:
            library["blocks_middle"][marker] = blocksm
    return library
#parse_library_ini
