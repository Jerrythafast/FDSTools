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

"""
Convert between TSSV (tab-separated) and FDSTools (ini-style) library
formats.

This is a convenience tool for users migrating from the standalone
'TSSV' programme.  Use the 'library' tool if you wish to create a new,
empty FDSTools library file to start with.

Both FDSTools and the standalone 'TSSV' programme use a library file to
store the names, flanking (primer) sequences, and STR repeat structure
of forensic STR markers.  However, the TSSV library file format is not
suitable for non-STR markers and automatic generation of allele names.
FDSTools therefore employs a different (ini-style) library file format
that can store more details about the markers used.  The libconvert tool
can be used to convert between the two formats.

Please refer to the help of the 'library' tool for more information
about FDSTools library files.
"""
import argparse
import sys
import re
import os.path

from ..lib import parse_library, INI_COMMENT
from library import make_empty_library_ini
from ConfigParser import RawConfigParser

__version__ = "1.1.0"


def convert_library(infile, outfile, aliases=False):
    pattern_reverse = re.compile("\(([ACGT]+)\)\{(\d+),(\d+)\}")
    if infile is not None:
        library = parse_library(infile, stream=True)
    else:
        library = {"flanks": {}, "regex": {}}
    if "aliases" in library:
        # FDSTools -> TSSV
        markers = set()
        for marker in library["flanks"]:
            markers.add(marker)
        for marker in library["prefix"]:
            markers.add(marker)
        for marker in library["suffix"]:
            markers.add(marker)
        for marker in library["regex"]:
            markers.add(marker)
        for marker in library["nostr_reference"]:
            markers.add(marker)

        marker_aliases = {}
        for alias in library["aliases"]:
            marker = library["aliases"][alias]["marker"]
            markers.add(marker)
            try:
                marker_aliases[marker].append(alias)
            except KeyError:
                marker_aliases[marker] = [alias]

        for marker in sorted(markers):
            pattern = []
            if marker in library["aliases"] and not aliases:
                # Ignore this alias, it will be merged into its marker.
                continue
            if marker in library["aliases"] and aliases:
                # Output this alias as a separate marker.
                if marker in library["flanks"]:
                    flanks = library["flanks"][marker]
                elif library["aliases"][marker]["marker"] in library["flanks"]:
                    flanks = library["flanks"][
                        library["aliases"][marker]["marker"]]
                else:
                    continue  # Worthless, no flanks.
                if marker in library["regex"]:
                    pattern = pattern_reverse.findall(
                        library["regex"][marker].pattern)
            elif aliases or marker not in marker_aliases:
                # Normal marker, or separately from its aliases.
                if marker not in library["flanks"]:
                    continue  # Worthless, no flanks.
                flanks = library["flanks"][marker]
                if marker in library["regex"]:
                    pattern = pattern_reverse.findall(
                        library["regex"][marker].pattern)
                elif marker in library["nostr_reference"]:
                    pattern = [(library["nostr_reference"][marker], "1", "1")]
            else:
                # Merge marker with its aliases.
                flanks = False
                if marker in library["flanks"]:
                    flanks = library["flanks"][marker]
                else:
                    for alias in marker_aliases[marker]:
                        if alias in library["flanks"]:
                            flanks = library["flanks"][alias]
                            break
                if not flanks:
                    continue  # Worthless, no flanks.

                prefixes = []
                suffixes = []
                if marker in library["prefix"]:
                    prefixes += library["prefix"][marker]
                if marker in library["suffix"]:
                    suffixes += library["suffix"][marker]
                middle = []

                if marker in library["regex"]:
                    # This marker has a regex next to its aliases.
                    # Check if the aliases fit the regex without change.
                    unmatched = []
                    for alias in marker_aliases[marker]:
                        allele = []
                        if marker in library["prefix"]:
                            allele.append(library["prefix"][marker][0])
                        allele.append(library["aliases"][alias]["sequence"])
                        if marker in library["suffix"]:
                            allele.append(library["suffix"][marker][0])
                        allele = "".join(allele)
                        if library["regex"][marker].match(allele) is None:
                            unmatched.append(
                                library["aliases"][alias]["sequence"])

                    middle = pattern_reverse.findall(
                        library["regex"][marker].pattern)[len(prefixes):]
                    if len(suffixes):
                        middle = middle[:-len(suffixes)]
                    if unmatched:
                        middle = [(x, "0", "1") for x in unmatched] + \
                                 [(x[0], "0", x[2]) for x in middle]

                elif marker in library["nostr_reference"]:
                    middle = [(library["nostr_reference"][marker],
                        "0" if marker in marker_aliases else "1", "1")]

                # Add prefixes and suffixes of aliases.
                if marker in marker_aliases:
                    for alias in marker_aliases[marker]:
                        if alias in library["prefix"]:
                            prefixes += library["prefix"][alias]
                        if alias in library["suffix"]:
                            suffixes += library["suffix"][alias]
                        if marker not in library["regex"] and (
                                marker not in library["nostr_reference"] or
                                library["nostr_reference"][marker] !=
                                    library["aliases"][alias]["sequence"]):
                            middle.append((
                                library["aliases"][alias]["sequence"],
                                "0", "1"))

                # Final regex is prefixes + middle + suffixes.
                pattern = []
                for i in range(len(prefixes)):
                    if i == prefixes.index(prefixes[i]):
                        pattern.append((prefixes[i], "0", "1"))
                pattern += middle
                for i in range(len(suffixes)):
                    if i == suffixes.index(suffixes[i]):
                        pattern.append((suffixes[i], "0", "1"))

            outfile.write("%s\t%s\t%s\t%s\n" % (
                marker, flanks[0], flanks[1],
                " ".join("%s %s %s" % x for x in pattern)))

    else:
        # TSSV -> FDSTools
        outfile.write(INI_COMMENT.fill("Lines beginning with a semicolon (;) "
            "are ignored by FDSTools.") + "\n\n")
        ini = make_empty_library_ini("full", aliases)

        # Enter flanking sequences and STR definitions.
        fmt = "%%-%is" % reduce(max, map(len,
            set(library["flanks"].keys() + library["regex"].keys())), 0)
        for marker in sorted(library["flanks"]):
            ini.set("flanks", fmt%marker, ", ".join(library["flanks"][marker]))
        for marker in sorted(library["regex"]):
            blocks = pattern_reverse.findall(library["regex"][marker].pattern)
            ini.set("repeat", fmt % marker,
                    " ".join("%s %s %s" % x for x in blocks))

            # Try to infer block length from the regular expression.
            length_counts = {0: 0}
            for block in blocks:
                amount = (int(block[1]) + int(block[2])) / 2.
                if len(block[0]) not in length_counts:
                    length_counts[len(block[0])] = amount
                else:
                    length_counts[len(block[0])] += amount
            block_length = sorted(
                length_counts, key=lambda x: -length_counts[x])[0]
            if block_length != 0 and block_length < 10:
                ini.set("block_length", fmt % marker, block_length)

            # Write max_expected_copies=2 for all markers explicitly.
            ini.set("max_expected_copies", fmt % marker, 2)

            # TODO: I could also do some fiddling for prefix/suffix...

        # Write INI file.
        ini.write(outfile)
#convert_library


def add_arguments(parser):
    parser.add_argument('infile', nargs='?', metavar="IN", default=sys.stdin,
        help="input library file, the format is automatically detected "
             "(default: read from stdin)")
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-a', '--aliases', action="store_true",
        help="when converting to TSSV format, aliases in FDSTools libraries "
             "are converted to separate markers in the output library when "
             "this option is specified; otherwise, they are merged into their "
             "respective markers; when converting to FDSTools format, the "
             "[aliases] section is included in the output if this option is "
             "specified")
#add_arguments


def run(args):
    if args.infile == "-":
        args.infile = sys.stdin
    if (args.infile != sys.stdin and args.outfile == sys.stdout
            and not os.path.exists(args.infile)):
        # One filename given, and it does not exist.  Assume outfile.
        args.outfile = open(args.infile, 'w')
        args.infile = sys.stdin

    if args.infile != sys.stdin:
        # Open the specified input file.
        args.infile = open(args.infile, 'r')
    elif args.infile.isatty():
        # No input given.  Produce a default FDSTools library.
        args.infile = None

    convert_library(args.infile, args.outfile, args.aliases)
#run
