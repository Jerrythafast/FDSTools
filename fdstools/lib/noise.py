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

from itertools import islice
from .io import get_column_ids
from .seq import ensure_sequence_format, mutate_sequence


def load_profiles(profilefile, library=None):
    # TODO: Docs.
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return {}  # Empty file.
    colid_marker, colid_allele, colid_sequence = get_column_ids(
        column_names, "marker", "allele", "sequence")

    colid_tool, colid_tools, colid_fmean, colid_rmean, colid_tmean = get_column_ids(
        column_names, "tool", "tools", "fmean", "rmean", "tmean", optional=True)

    if colid_fmean is None and colid_rmean is None and colid_tmean is None:
        raise ValueError("Invalid background noise profiles file: the columns 'fmean', 'rmean', "
            "and 'tmean' were not found; at least one of these must be present")

    profiles = {}
    for line in profilefile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        marker = line[colid_marker]
        if marker not in profiles:
            profiles[marker] = {}
        allele = ensure_sequence_format(line[colid_allele], "raw", library=library, marker=marker)
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
            library=library, marker=marker)
        if allele not in profiles[marker]:
            profiles[marker][allele] = {}
        elif sequence in profiles[marker][allele]:
            raise ValueError(
                "Invalid background noise profiles file: encountered multiple values for marker "
                "'%s' allele '%s' sequence '%s'" % (marker, allele, sequence))
        tools = set([] if colid_tools is None else map(str.strip, line[colid_tools].split(",")))
        if colid_tool is not None:
            tools.add(line[colid_tool].strip())
        profiles[marker][allele][sequence] = {
            "tools": tools,
            "forward": float(line[colid_fmean]) if colid_fmean is not None else 0.,
            "reverse": float(line[colid_rmean]) if colid_rmean is not None else 0.,
            "total": float(line[colid_tmean]) if colid_tmean is not None else 0.}

    return profiles
#load_profiles


def iter_overlapping_repeats(allele, repeat_pattern, unit_length):
    """Similar to repeat_pattern.finditer(allele) but allowing overlap."""
    pos = 0
    while True:
        match = repeat_pattern.search(allele, pos)
        if not match:
            return
        yield match
        pos = match.end() - unit_length + 1
#iter_overlapping_repeats


def get_all_combinations(stutters, minlen, maxlen, indexes_out=False, appendto=None, start=0):
    """
    Generate all compatible combinations of stutters.
    The input stutters must be sorted by increasing (start, end).
    """
    if appendto is None:
        appendto = []
    for i, stutter_i in enumerate(islice(stutters, start, None), start):
        if appendto:
            prev_stutter = stutters[appendto[-1]] if indexes_out else appendto[-1]
            if stutter_i["end"] <= prev_stutter["end"]:
                continue  # Engulfed in or same repeat stretch as previous.
            if prev_stutter["fold"] < 0 and -prev_stutter["fold"] * len(prev_stutter["seq"]) > \
                    stutter_i["start"] - prev_stutter["start"]:
                continue  # Minus stutter of prev overlaps next stretch.
            if stutter_i["fold"] < 0 and -stutter_i["fold"] * len(stutter_i["seq"]) > \
                    stutter_i["end"] - prev_stutter["end"]:
                continue  # Minus stutter overlaps prev stretch.
        withnewelement = [x for x in appendto]
        withnewelement.append(i if indexes_out else stutter_i)
        if len(withnewelement) >= minlen:
            yield tuple(withnewelement)
        if len(withnewelement) < maxlen:
            yield from get_all_combinations(stutters, minlen, maxlen, indexes_out, withnewelement, i + 1)
#get_all_combinations


def get_stutter_sequence(allele, stutter_combination):
    return mutate_sequence(allele, [
        "%+i.1->%s" % (s["end"], allele[s["end"] - s["stutlen"] : s["end"]])
        if s["fold"] > 0 else "%+i%s>-" % (
            s["end"] - s["stutlen"] + 1,
            allele[s["end"] - s["stutlen"] : s["end"]])
        for s in stutter_combination])
#get_stutter_sequence