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

"""
Create an empty FDSTools library file.

An FDSTools library file contains various details about the forensic
markers used in the analysis, such as the genomic location, expected
number of alleles, expected length range of alleles, etc.  FDSTools
primarily uses library files for configuring STRNaming, which is
responsible for converting sequences to allele names and vice versa.
This is true even for non-STR markers and fragments on the
mitochondrial genome.

In its simplest form, the library file only contains the positions (on
the human genome reference sequence, GRCh38) of the reported range of
each marker.  This is referred to as a 'smart' library file.
Alternatively, markers can be explicitly configured, which was the
default prior to FDSTools version 2.0.  Explicit configuration is
currently required when the analysed markers are non-human.

Users migrating from the standalone 'TSSV' programme may use the
libconvert tool to convert their TSSV library file to FDSTools format.
"""
import argparse
import sys

from configparser import RawConfigParser
from errno import EPIPE
from types import MethodType

from ..lib.library import BUILTIN_LIBS, BUILTIN_NAMES, INI_COMMENT

__version__ = "1.1.1"


def ini_add_comment(ini, section, comment):
    for line in INI_COMMENT.wrap(comment):
        ini.set(section, line)
#ini_add_comment


def make_empty_library_ini(type, microhaplotypes=False):
    ini = RawConfigParser(allow_no_value=True)
    ini.optionxform = str
    ini.add_comment = MethodType(ini_add_comment, ini)

    # Create sections and add comments to explain how to use them.
    ini.add_section("genome_position")
    ini.add_comment("genome_position", #smart, str, non-str, full
        "Specify the chromosome number and positions of the first and last reported nucleotide of "
        "each marker (both inclusive, using human genome build GRCh38%s). This range should not "
        "include the primer binding sites.%s" % (
            " and rCRS for human mtDNA" if type != "str" else "",
            " This section is required for automatic configuration of markers; it is optional "
            "for markers that are explictily configured in this library file."
                if type != "smart" else ""))
    if type != "str":
        ini.add_comment("genome_position",
            "Specify 'M' as the chromosome name for markers on mitochondrial "
            "DNA. Allele names generated for these markers will follow mtDNA "
            "nomenclature guidelines (Parson et al., 2014). If one of your "
            "mtDNA fragments starts near the end of the reference sequence "
            "and continues at the beginning, you can obtain correct base "
            "numbering by specifying the fragment's genome position as \"M, "
            "(starting position), 16569, 1, (ending position)\". This tells "
            "FDSTools that the marker is a concatenation of two fragments, "
            "where the first fragment ends at position 16569 and the second "
            "fragment starts at position 1. Similarly, for a fragment that "
            "spans position 3107 in the rCRS (which is nonexistent), you may "
            "specify \"M, (starting position), 3106, 3108, (ending "
            "position)\".")
    if microhaplotypes or type == "full":
        ini.add_section("microhaplotype_positions")
        ini.add_comment("microhaplotype_positions",
            "For each microhaplotype marker, specify one or more positions of SNPs that should "
            "be reported as part of the microhaplotype.%s" % (
                " If the [genome_position] of the marker is given, positions must be within the "
                "given range. Otherwise, the reference sequence must be explicitly provided in "
                "the [no_repeat] section position 1 refers to the first base in the reference "
                "sequence."
                    if type in ("non-str", "full") else ""))
    ini.add_section("flanks")
    ini.add_comment("flanks",
        "The TSSV tool will use a pair of short anchor sequences just outside the reported range "
        "of each marker (e.g., primer sequences) to identify which input reads correspond to "
        "which marker. %s The sequence may contain IUPAC codes for ambiguous positions to account "
        "for degenerate bases in the primers or for bisulfite-converted targets in methylation-"
        "based studies (e.g., Y matches either C or T)." % (
            "Specify two comma-separated values: left flank and right flank sequence, in the same "
            "sequence orientation (strand)."
                if type in ("str", "non-str") else
            ("The default length of the anchor sequences used can be specified as an argument to "
            "the TSSV tool. Individual alternative lengths can be specified here for each marker. "
            "Specify two comma-separated values: one for the left and one for the right flank. "
            "The value can be a number (the length of sequence to use) or an explicit anchor "
            "sequence to use.%s" % (
                " For markers configured explicitly in this library file, the anchor sequences "
                "must be specified explicitly as well." if type == "full" else ""))))
    ini.add_section("max_expected_copies")
    ini.add_comment("max_expected_copies",
        "By default, the Allelefinder tool will report up to 2 alleles per marker, but only a "
        "single allele for markers %son the Y chromosome. If this is incorrect, specify the "
        "maximum expected number of copies (i.e., alleles) for each marker in a "
        "single-contributor reference sample here." % (
            "on mitochondrial DNA or " if type != "str" else ""))
    ini.add_section("expected_allele_length")
    ini.add_comment("expected_allele_length",
        "Specify one or two values for each marker. The first value gives the "
        "expected minimum length (in nucleotides, %sexcluding flanks) of the "
        "alleles and the second value (if given) specifies the maximum allele "
        "length expected for that marker (both inclusive). The TSSV tool will filter "
        "sequences that have a length outside this range." %
        ("including prefix and suffix, " if type in ("str", "full") else ""))

    if type in ("str", "full"):
        ini.add_section("prefix")
        ini.add_comment("prefix",
            "For explicitly-configured STR markers: Specify the prefix sequence of each STR "
            "marker. The prefix is the sequence between the left flank and the repeat and is "
            "omitted from allele names. The sequence is used as the reference sequence for that "
            "marker when generating allele names. Deviations from the reference are expressed as "
            "variants.")
        ini.add_section("suffix")
        ini.add_comment("suffix",
            "For explicitly-configured STR markers: Specify the suffix sequence of each STR "
            "marker. The suffix is the sequence between the repeat and the right flank. The "
            "sequence is used as the reference sequence for that marker when generating allele "
            "names.")
        ini.add_section("repeat")
        ini.add_comment("repeat",
            "For explicitly-configured STR markers: Specify the repeat structure of each STR "
            "marker in space-separated triples of sequence, minimum number of repeats, and "
            "maximum number of repeats.")
        ini.add_section("length_adjust")
        ini.add_comment("length_adjust",
            "For explicitly-configured STR markers: To prevent discrepancies between traditional "
            "CE allele numbers and the CE number in FDSTools allele names, the CE allele number "
            "as calculated by FDSTools is based on the length of the repeat sequence minus the "
            "adjustment specified here.")
        ini.add_section("block_length")
        ini.add_comment("block_length",
            "For explicitly-configured STR markers: Specify the repeat unit length of each STR "
            "marker. By default, the length of the repeat unit of the longest repeat is used.")
    if type in ("non-str", "full"):
        ini.add_section("no_repeat")
        ini.add_comment("no_repeat",
            "For explicitly-configured non-STR markers: Specify the reference sequence for each "
            "non-STR marker.")
    return ini
#make_empty_library_ini


def create_library(outfile, type, microhaplotypes=False, builtin=None):
    outfile.write(INI_COMMENT.fill(
        "Lines beginning with a semicolon (;) are ignored by FDSTools.") + "\n\n")
    library = make_empty_library_ini(type, microhaplotypes)
    if builtin is not None:
        with BUILTIN_LIBS[builtin].open("rt", encoding="UTF-8") as handle:
            library.read_file(handle, source=builtin)
    library.write(outfile)
#create_library


def add_arguments(parser):
    parser.add_argument("outfile", nargs="?", metavar="OUT",
        default=sys.stdout, type=argparse.FileType("tw", encoding="UTF-8"),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument("-t", "--type", metavar="TYPE", default="smart",
        choices=("smart", "full", "str", "non-str"),
        help="the type of markers that this library file will be used for; with 'smart' (the "
             "default), only the genomic positions of the analysed ranges (i.e., the amplicon "
             "excluding the primers) need to be specified and FDSTools will automatically detect "
             "and configure allele naming using STRNaming (currently only supported for markers "
             "in the human genome); 'full' will create a library file with all possible sections; "
             "'str' or 'non-str' will only output sections used to explicitly define STR and "
             "non-STR markers, respectively")
    parser.add_argument("-m", "--microhaplotypes", action="store_true",
        help="if specified, the [microhaplotype_positions] section is included, which can be "
             "used to configure allele calling for microhaplotype targets")
    parser.add_argument("-b", "--builtin", metavar="NAME", choices=BUILTIN_NAMES,
        help="start with a built-in library file, choose from '%s'" % "', '".join(BUILTIN_NAMES))
#add_arguments


def run(args):
    try:
        create_library(args.outfile, args.type, args.microhaplotypes, args.builtin)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#run
