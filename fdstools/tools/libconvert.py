#!/usr/bin/env python
"""
Convert between TSSV (tab-separated) and FDSTools (ini-style) library
formats.  When no input is given, an empty FDSTools library is produced.

FDSTools uses library files to convert sequences to allele names and
vice versa.  Because FDSTools was made to work with the output of TSSV,
it has been made compatible with TSSV's library files as well.  The TSSV
library format is much less suitable for the generation of allele names,
which makes FDSTools libraries the recommended choice.  The libconvert
tool can be used to create a compatible TSSV library.

In FDSTools, sequences of STR alleles are split up into three parts: a
prefix, the STR, and a suffix.  The prefix and suffix are optional and
are meant to fill the gap between the STR and the primer binding sites.
The primer binding sites are called 'flanks' in the library file.

Allele names typically consist of an allele number compatible with those
obtained from Capillary Electrophoresis (CE), followed by the STR
sequence in a shortened form and any substitutions or other variants
that occur in the prefix and suffix.  The first prefix/suffix in the
library file is used as the reference sequence for calling variants.

Special alleles, such as the 'X' and 'Y' allele from the Amelogenin
gender test, may be given an explicit allele name by specifying an Alias
in the FDSTools library file.

Run libconvert without any arguments to obtain a default FDSTools
library to start with.  The default library contains commentary lines
that explain the use of each section in more detail.
"""
import argparse
import sys
import re
import os.path

from ..lib import parse_library
from ConfigParser import RawConfigParser
from StringIO import StringIO

__version__ = "0.1dev"


# If no input is given, convert the following to FDSTools format.
_DEFAULT_LIBRARY = "\t".join([
    "MyMarker",
    "ACTAGCTAGCGCTA",
    "GCTCGATCGATCGA",
    "TGAT 0 2 AGAT 3 20 ACCT 0 5"])


def convert_library(infile, outfile, aliases=False):
    pattern_reverse = re.compile("\(([ACGT]+)\)\{(\d+),(\d+)\}")
    library = parse_library(infile)
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
                # Normal marker, or separtely from its aliases.
                if marker not in library["flanks"]:
                    continue  # Worthless, no flanks.
                flanks = library["flanks"][marker]
                if marker in library["regex"]:
                    pattern = pattern_reverse.findall(
                        library["regex"][marker].pattern)
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
                        middle = [(x[0], "0", x[2]) for x in middle] + \
                                 [(x, "0", "1") for x in unmatched]

                # Add prefixes and suffixes of aliases.
                if marker in marker_aliases:
                    for alias in marker_aliases[marker]:
                        if alias in library["prefix"]:
                            prefixes += library["prefix"][alias]
                        if alias in library["suffix"]:
                            suffixes += library["suffix"][alias]
                        if marker not in library["regex"]:
                            middle.append((
                                library["aliases"][alias]["sequence"],
                                "0", "1"))

                # Final regex is prefixes + middle + suffixes.
                pattern = []
                for prefix in prefixes:
                    pattern.append((prefix, "0", "1"))
                pattern += middle
                for suffix in suffixes:
                    pattern.append((suffix, "0", "1"))

            outfile.write("%s\t%s\t%s\t%s\n" % (
                marker, flanks[0], flanks[1],
                " ".join("%s %s %s" % x for x in pattern)))

    else:
        # TSSV -> FDSTools
        ini = RawConfigParser(allow_no_value=True)
        ini.optionxform = str

        # Create sections.  Most of them will be empty but we will put
        # comments in them to explain how to use them.
        ini.add_section("aliases")
        ini.set("aliases", "; Specify three comma-separated values: marker "
                            "name, sequence, and allele name.")
        ini.set("aliases", "; You may use the alias name to specify flanks, "
                           "prefix, and suffix for this")
        ini.set("aliases", "; allele specifically. You cannot specify a "
                           "repeat structure for an alias.")
        ini.set("aliases", ";MyAlias = MyMarker, AGCTAGC, MySpecialAlleleName")
        ini.add_section("flanks")
        ini.set("flanks", "; Specify two comma-separated values: left flank "
                          "and right flank.")
        ini.add_section("prefix")
        ini.set("prefix", "; Specify all known prefix sequences separated "
                          "by commas. The first sequence")
        ini.set("prefix", "; listed is used as the reference sequence when "
                          "generating allele names. The")
        ini.set("prefix", "; prefix is the sequence between the left flank "
                          "and the repeat and is omitted")
        ini.set("prefix", "; from allele names. Deviations from the reference "
                          "are expressed as variants.")
        ini.add_section("suffix")
        ini.set("suffix", "; Specify all known suffix sequences separated "
                          "by commas. The first sequence")
        ini.set("suffix", "; listed is used as the reference sequence when "
                          "generating allele names. The")
        ini.set("suffix", "; suffix is the sequence between the repeat and "
                          "the right flank.")
        ini.add_section("repeat")
        ini.set("repeat", "; Specify the STR repeat structure in "
                          "space-separated triples of sequence,")
        ini.set("repeat", "; minimum number of repeats, and maximum number of "
                          "repeats.")
        ini.add_section("length_adjust")
        ini.set("length_adjust", "; When generating allele names, the CE "
                                 "allele number is based on the length")
        ini.set("length_adjust", "; of the sequence (prefix+repeat+suffix) "
                                 "minus the adjustment specified here.")
        ini.add_section("block_length")
        ini.set("block_length", "; Specify the core repeat unit length of "
                                "each marker. The default length is 4.")

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
                amount = (int(block[1])+int(block[2]))/2.
                if len(block[0]) not in length_counts:
                    length_counts[len(block[0])] = amount
                else:
                    length_counts[len(block[0])] += amount
            block_length = sorted(
                length_counts, key=lambda x: -length_counts[x])[0]
            if block_length != 0 and block_length < 10:
                ini.set("block_length", fmt%marker, block_length)

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
        help="if specified, aliases in FDSTools libraries are converted to "
             "separate markers in the output library; otherwise, they are "
             "merged into their respective markers")
#add_arguments


def run(args):
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
        args.infile = StringIO(_DEFAULT_LIBRARY)

    convert_library(args.infile, args.outfile, args.aliases)
#run


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    try:
        add_arguments(parser)
        run(parser.parse_args())
    except OSError as error:
        parser.error(error)
#main


if __name__ == "__main__":
    main()
