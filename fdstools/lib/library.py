#!/usr/bin/env python3

#
# Copyright (C) 2023 Jerry Hoogenboom
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
import sys
import textwrap

from configparser import RawConfigParser, MissingSectionHeaderError
from pathlib import Path
from strnaming import classes, libstrnaming, libsequence

from .seq import PAT_SEQ_RAW, PAT_SEQ_IUPAC

# Patterns that match (parts of) an STR definition.
PAT_STR_DEF = re.compile("^(?:(?:(?<=^)|(?<!^)\s+)[ACGT]+\s+\d+\s+\d+)*$")
PAT_STR_DEF_BLOCK = re.compile("([ACGT]+)\s+(\d+)\s+(\d+)")

# Pattern that matches a chromosome name/number.
PAT_CHROMOSOME = re.compile("^(?:[Cc][Hh][Rr](?:[Oo][Mm])?)?([1-9XYM]|1\d|2[0-2])$")

# Pattern to split a comma-, tab- or space-separated list.
PAT_SPLIT = re.compile("\s*[, \t]\s*")

# TextWrapper object for formatting help texts in generated INI files.
INI_COMMENT = textwrap.TextWrapper(width=79, initial_indent="; ",
    subsequent_indent="; ", break_on_hyphens=False)

# Dict of built-in library files.
BUILTIN_LIBS = {_.stem: _
    for _ in (Path(__file__).parent.parent / "data" / "libraries").glob("*.ini")}
BUILTIN_NAMES = tuple(sorted(BUILTIN_LIBS))


def add_legacy_range(reported_range_store, marker, prefix, suffix, blocks, options, genome_position=None):
    """
    Add a new marker to the provided STRnaming ReportedRangeStore,
    defined only by the provided TSSV block structure. The overlong gap
    (if any) is automatically detected in the block structure.

    A fake refence sequence will be generated. Not necessarily one that
    conforms to the input block pattern, but it will include the flanks,
    prefix, suffix, stretches of repeat units and the overlong gap if
    detected. The flanks are read from options["flanks"].

    The genome_position should be a two-tuple (chromosome, position) and
    defaults to (marker, len(flanks[0]) + 1). If no blocks are specified,
    the genome_position can consist of multiple ranges by specifying
    (chromosome, start1, end1, start2, end2, ...)
    """
    # Make sure the following are met to make STRNaming work:
    # # Construct reference sequence that starts with flank+prefix and ends with suffix+flank.
    # # Repeats must cover minimum structure length.
    # # Each given unit occurs once; stretch length (in bp) reflects relative number of repeats in the input.
    # # Overlong gap (if any) is inserted between two stretches.
    # TODO:  Strictly speaking, a gap of >20nt should result in multiple structures.
    if genome_position is None:
        # NOTE: For legacy ranges with no genome_position, the reported
        # range starts at position 1. This necessarily means the left
        # flank, if any, goes into nonpositive positions, which is OK.
        genome_position = (marker, 1)
    if blocks and len(genome_position) > 3:
        raise ValueError(
            "Gapped or combined genomic range is not supported for STR marker %s" % marker)

    # Find units and overlong_gap.
    overlong_gap = ""
    units = {}
    for unit, min_repeats, max_repeats in blocks:
        if max_repeats == 1:
            if len(unit) > max(libstrnaming.NAMING_OPTIONS["max_gap"], len(overlong_gap)):
                overlong_gap = unit
        else:
            units[unit] = units.get(unit, 0) + max_repeats

    # Start building refseq and stretches.
    flanks = options["flanks"]
    refseq = [flanks[0], prefix]
    stretches = []
    start = genome_position[1]
    end = start + len(prefix)
    min_required_length = 1
    units = sorted((total_count, len(unit), unit) for unit, total_count in units.items())
    for total_count, len_unit, unit in units:
        if stretches and overlong_gap and not overlong_gap.startswith(stretches[-1][2]) and not overlong_gap.endswith(unit):
            # This would be an ideal place to insert the overlong gap.
            refseq.append(overlong_gap)
            end += len(overlong_gap)
            overlong_gap = ""  # To avoid adding it again.
        repeats = (min_required_length // len_unit or 1) + 1
        stretch_seq = unit * repeats
        stretch_length = len(stretch_seq)
        refseq.append(stretch_seq)
        stretches.append([end, end + stretch_length, unit])
        end += stretch_length
        min_required_length = stretch_length + 1

    # Make sure the overlong gap is inserted into the structure.
    if overlong_gap:
        refseq.insert(-1, overlong_gap)
        overlong_length = len(overlong_gap)
        end += overlong_length
        if len(stretches) == 1:
            # Duplicate the stretch if there is just one (otherwise it's not a gap).
            refseq.insert(-2, refseq[-1])
            stretch_length = stretches[0][1] - stretches[0][0]
            stretches.append([end, end + stretch_length, stretches[0][2]])
            end += stretch_length
        else:
            # Move the positions of the last stretch.
            stretches[-1][0] += overlong_length
            stretches[-1][1] += overlong_length

    if stretches:
        # Make sure the structure is long enough to be accepted by STRNaming.
        while end - stretches[0][0] < libstrnaming.NAMING_OPTIONS["min_structure_length"]:
            stretches[-1][1] += len(stretches[-1][2])
            end += len(stretches[-1][2])
            refseq[-1] += stretches[-1][2]

    # Add suffix.
    end += len(suffix)
    refseq.append(suffix)
    refseq.append(flanks[1])
    refseq = "".join(refseq)

    # Add data to the reported range store.
    struct_store = reported_range_store.get_structure_store()
    refseq_store = struct_store.get_refseq_store()
    pos = 0
    for i in range(1, len(genome_position), 2):
        location = genome_position[i]
        length = len(refseq) if i + 2 >= len(genome_position) else genome_position[i+1] - genome_position[i] + 1
        if i == 1:
            location -= len(flanks[0])
            length += len(flanks[0])
        refseq_store.add_refseq(genome_position[0], location, refseq[pos:pos+length])
        pos += length
    if stretches:
        struct_store.add_structure(genome_position[0],
            [[start, end, len(unit)] for start, end, unit in stretches])
    if len(genome_position) > 3:
        return reported_range_store.add_complex_range(marker, genome_position, options=options)
    return reported_range_store.add_range(marker, genome_position[0], start, end, options=options)
#add_legacy_range


def parse_library(handle):
    """Parse an FDSTools library file as made with the library tool."""
    markers = {}
    ini = RawConfigParser()
    ini.optionxform = str
    ini.readfp(handle)
    for section in ini.sections():
        for marker in ini.options(section):
            value = ini.get(section, marker).split(";", 1)[0].strip()
            section_low = section.lower()
            if section_low == "flanks":
                value = PAT_SPLIT.split(value)
                if len(value) != 2:
                    raise ValueError(
                        "For marker %s, %i flanking sequences were given, "
                        "need exactly 2" % (marker, len(value)))
                for i, val in enumerate(value):
                    if PAT_SEQ_IUPAC.match(val) is None:
                        try:
                            value[i] = int(val)
                            if value[i] < 1:
                                raise ValueError
                        except:
                            raise ValueError(
                                "Flanking sequence '%s' of marker %s is invalid" % (val, marker))
            elif section_low in ("prefix", "suffix"):
                if not PAT_SEQ_RAW.match(value):
                    raise ValueError(
                        "The %s sequence '%s' of marker %s is invalid (note: only one sequence "
                        "needs to be specified for reference)" % (section_low, value, marker))
            elif section_low == "genome_position":
                value = PAT_SPLIT.split(value)
                chromosome = PAT_CHROMOSOME.match(value[0])
                if chromosome is None:
                    raise ValueError("Invalid chromosome '%s' for marker %s" % (value[0], marker))
                pos = [chromosome.group(1)]
                for i in range(1, len(value)):
                    try:
                        pos.append(int(value[i]))
                    except:
                        raise ValueError(
                            "Position '%s' of marker %s is not a valid integer"
                            % (value[i], marker))
                    if not i % 2 and pos[-2] >= pos[-1]:
                        raise ValueError(
                            "End position %i of marker %s must be higher than "
                            "corresponding start position %i" % (pos[-1], marker, pos[-2]))
                if len(value) == 1:
                    pos.append(1)
                value = tuple(pos)
            elif section_low == "length_adjust":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Length adjustment '%s' of marker %s is not a valid "
                        "integer" % (value, marker))
            elif section_low == "block_length":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Block length '%s' of marker %s is not a valid integer" % (value, marker))
            elif section_low == "max_expected_copies":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Maximum number of expected copies '%s' of marker %s "
                        "is not a valid integer" % (value, marker))
            elif section_low == "repeat":
                if PAT_STR_DEF.match(value) is None:
                    raise ValueError(
                        "STR definition '%s' of marker %s is invalid" % (value, marker))
            elif section_low == "no_repeat":
                if PAT_SEQ_RAW.match(value) is None:
                    raise ValueError(
                        "Reference sequence '%s' of marker %s is invalid" % (value, marker))
            elif section_low == "microhaplotype_positions":
                value = PAT_SPLIT.split(value)
                for i in range(len(value)):
                    try:
                        value[i] = int(value[i])
                    except:
                        raise ValueError(
                            "Invalid position number '%s' for microhaplotype marker %s" %
                            (value[i], marker))
            elif section_low == "expected_allele_length":
                value = PAT_SPLIT.split(value)
                try:
                    min_length = int(value[0])
                except:
                    raise ValueError(
                        "Minimum expected allele length '%s' of marker %s "
                        "is not a valid integer" % (value[0], marker))
                if len(value) == 2:
                    try:
                        max_length = int(value[1])
                    except:
                        raise ValueError(
                            "Maximum expected allele length '%s' of marker %s "
                            "is not a valid integer" % (value[1], marker))
                elif len(value) > 2:
                    raise ValueError(
                        "%i values specified for expected_allele_length of marker %s; specify "
                        "only a minimum and optionally a maximum length" % (len(value), marker))
                else:
                    max_length = sys.maxsize
                value = (min_length, max_length)

            # Store the validated value.
            if marker not in markers:
                markers[marker] = {}
            markers[marker][section_low] = value

    # Create a ReportedRangeStore to store data about each marker.
    reported_range_store = classes.ReportedRangeStore()
    MUTEX_GROUPS = {
        "explicit STR": ("prefix", "suffix", "repeat", "length_adjust", "block_length"),
        "explicit non-STR": ("no_repeat",),
        "STRNaming-specific": tuple()}  # NOTE: No STRNaming-specific settings now...

    # Get a list of ranges, so that we can load the STR structures in one call.
    ranges = {}
    for marker, settings in markers.items():
        groups = [mutex_group for mutex_group, sections in MUTEX_GROUPS.items()
            if any(section in settings for section in sections)]
        if len(groups) > 1:
            raise ValueError("The definition of marker %s is ambiguous, because it appears in %s" %
                (marker, " and ".join(
                    "%s sections (%s)" % (group, ", ".join(
                        section for section in MUTEX_GROUPS[group] if section in settings))
                    for group in groups)))
        if "explicit STR" not in groups and "explicit non-STR" not in groups:
            try:
                genome_position = settings["genome_position"]
            except KeyError:
                raise ValueError("No genome_position or explicit repeat or no_repeat "
                                 "configuration provided for marker %s" % marker)
            if not len(genome_position) % 2:
                raise ValueError(
                    "Invalid genomic position given for marker %s: need an odd number of values "
                    "(chromosome, start position, end position[, start2, end2, ...])" % marker)
            if len(genome_position) == 3:
                chromosome, start, end = genome_position
                if chromosome not in ranges:
                    ranges[chromosome] = [(start, end + 1)]
                else:
                    ranges[chromosome].append((start, end + 1))

    # Load the STR structures on the given ranges.
    structure_store = reported_range_store.get_structure_store()
    for chromosome, ranges_here in ranges.items():
        structure_store.load_within_ranges(chromosome, sorted(ranges_here))

    # Now add all markers to the ReportedRangeStore.
    for marker, settings in markers.items():
        groups = [mutex_group for mutex_group, sections in MUTEX_GROUPS.items()
            if any(section in settings for section in sections)]
        options = {option: settings[option] for option in {"flanks", "max_expected_copies",
                "expected_allele_length", "microhaplotype_positions"} & settings.keys()}
        if "explicit STR" in groups:
            # Legacy FDSTools-style definition of an STR marker.
            # TODO: Alias of STR markers was defined as excluding the prefix/suffix!
            if "flanks" not in options:
                options["flanks"] = ("", "")
            elif any(isinstance(flank, int) for flank in options["flanks"]):
                raise ValueError(
                    "Please specify an explit flanking sequence, not just a length, for marker %s"
                        % marker)
            reported_range = add_legacy_range(reported_range_store, marker,
                settings.get("prefix", ""),
                settings.get("suffix", ""),
                [(unit, int(min_repeats), int(max_repeats)) for unit, min_repeats, max_repeats in
                    PAT_STR_DEF_BLOCK.findall(settings.get("repeat", ""))],
                options,
                settings.get("genome_position", None))
            if "length_adjust" in settings:
                reported_range.length_adjust -= settings["length_adjust"]
            if "block_length" in settings:
                reported_range.block_length = settings["block_length"]
        elif "explicit non-STR" in groups:
            # Legacy FDSTools-style definition of a non-STR marker.
            if "flanks" not in options:
                options["flanks"] = ("", "")
            elif any(isinstance(flank, int) for flank in options["flanks"]):
                raise ValueError(
                    "Please specify an explit flanking sequence, not just a length, for marker %s"
                        % marker)
            refseq = settings["no_repeat"]
            pos = None
            if "genome_position" in settings:
                pos = settings["genome_position"]

                # Sanity check: end position should reflect ref length.
                length = sum(pos[i] - pos[i - 1] + 1 for i in range(2, len(pos), 2))
                if len(refseq) < length or (len(pos) % 2 and len(refseq) != length):
                    raise ValueError(
                        "Length of reference sequence of marker %s is %i bases, but "
                        "genome positions add up to %i bases" % (marker, len(refseq), length))
            add_legacy_range(reported_range_store, marker, refseq, "", [], options, pos)
        else:
            # Use STRNaming for this marker.
            # TODO: Alias of STR markers was defined as excluding the prefix/suffix!
            genome_position = settings["genome_position"]
            if len(genome_position) == 3:
                chromosome, start, end = genome_position
                reported_range_store.add_range(
                    marker, chromosome, start, end + 1, options=options)
            else:
                reported_range_store.add_complex_range(marker, genome_position, options=options)

        if "microhaplotype_positions" in settings:
            # Put Ns in reporting range refseq for microhaplotype markers.
            reported_range = reported_range_store.get_range(marker)
            if reported_range.library:
                raise ValueError(
                    "Cannot define microhaplotype positions for STR marker %s" % marker)
            refseq = list(reported_range.refseq)
            location = reported_range.location
            for position in settings["microhaplotype_positions"]:
                refseq[libsequence.get_genome_pos(location, position, invert=True)] = "N"
            reported_range.refseq = "".join(refseq)
    return reported_range_store
#parse_library


def get_builtin_library(name):
    try:
        with BUILTIN_LIBS[name].open("rt", encoding="UTF-8") as handle:
            return parse_library(handle)
    except KeyError:
        return None
#get_builtin_library


def get_max_expected_alleles(max_alleles, marker, library):
    if max_alleles is not None:
        return max_alleles
    if library is not None:
        r = library.get_range(marker)
        return r.get_option("max_expected_copies", 1 if r.location[0] in ("MY") else 2)
    return 2
#get_max_expected_alleles
