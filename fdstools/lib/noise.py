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

from .io import get_column_ids
from .seq import ensure_sequence_format


def load_profiles(profilefile, library=None):
    # TODO: To be replaced with load_profiles_new (less memory).
    # TODO: Docs.
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return {}  # Empty file.
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
        allele = ensure_sequence_format(line[colid_allele], "raw", library=library, marker=marker)
        sequence = ensure_sequence_format(
            line[colid_sequence], "raw", library=library, marker=marker)
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
            list(profiles[marker]["m"] - profiles[marker]["n"])
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
    # TODO: Docs.
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return {}  # Empty file.
    colid_marker, colid_allele, colid_sequence, colid_tool = get_column_ids(
        column_names, "marker", "allele", "sequence", "tool")

    colid_fmean, colid_rmean, colid_tmean = get_column_ids(
        column_names, "fmean", "rmean", "tmean", optional=True)

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
        profiles[marker][allele][sequence] = {
            "tools": set(map(str.strip, line[colid_tool].split(","))),
            "forward": float(line[colid_fmean]) if colid_fmean is not None else 0.,
            "reverse": float(line[colid_rmean]) if colid_rmean is not None else 0.,
            "total": float(line[colid_tmean]) if colid_tmean is not None else 0.}

    return profiles
#load_profiles_new
