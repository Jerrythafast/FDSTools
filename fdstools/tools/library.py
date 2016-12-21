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
Create an empty FDSTools library file.

An FDSTools library file contains various details about the forensic
markers used in the analysis, such as the flanking (primer) sequences,
general STR structure or reference sequence of non-STR markers, genomic
location, expected number of alleles, expected length range of alleles,
etc.  FDSTools uses library files for tasks such as linking raw sequence
reads to markers and converting sequences to allele names or vice versa.

In FDSTools, sequences of STR alleles are split up into three parts: a
prefix, the STR, and a suffix.  The prefix and suffix are optional and
are meant to fill the gap between the STR and the primer binding sites.
The primer binding sites are called 'flanks' in the library file.  For
non-STR markers, FDSTools library files simply contain the reference
sequence of the region between the flanks.

Allele names typically consist of an allele number compatible with those
obtained from Capillary Electrophoresis (CE), followed by the STR
sequence in a shortened form and any substitutions or other variants
that occur in the prefix and suffix.  The first prefix/suffix in the
library file is used as the reference sequence for calling variants.

Special alleles, such as the 'X' and 'Y' allele from the Amelogenin
gender test, may be given an explicit allele name by specifying an Alias
in the FDSTools library file.

Users migrating from the standalone 'TSSV' programme may use the
libconvert tool to convert their TSSV library file to FDSTools format.
"""
import argparse
import sys

from ..lib import INI_COMMENT
from ConfigParser import RawConfigParser
from types import MethodType

__version__ = "1.0.2"


def ini_add_comment(ini, section, comment):
    for line in INI_COMMENT.wrap(comment):
        ini.set(section, line)
#ini_add_comment


def make_empty_library_ini(type, aliases=False):
    ini = RawConfigParser(allow_no_value=True)
    ini.optionxform = str
    ini.add_comment = MethodType(ini_add_comment, ini)

    # Create sections and add comments to explain how to use them.
    if aliases:
        ini.add_section("aliases")
        ini.add_comment("aliases",
            "Assign an explicit allele name to a specific sequence of a "
            "specific marker. Specify three comma-separated values: "
            "marker name, sequence, and allele name. You may use the "
            "alias name in the %s to specify them for this allele "
            "specifically.%s" % (
                "flanks section" if type == "non-str" else
                    "flanks, prefix, and suffix sections",
                " You cannot specify a repeat structure for an alias."
                    if type != "non-str" else ""))
        ini.set("aliases",
            ";AmelX = Amel, TCAGCTATGAGGTAATTTTTCTCTTTACTAATTTTGACCATTGTTTGCGT"
            "TAACAATGCCCTGGGCTCTGTAAAGAATAGTGTGTTGATTCTTTATCCCAGATGTTTCTCAAGTG"
            "GTCCTGATTTTACAGTTCCTACCACCAGCTTCCCA, X")
        ini.set("aliases",
            ";AmelY = Amel, TCAGCTATGAGGTAATTTTTCTCTTTACTAATTTTGATCACTGTTTGCAT"
            "TAGCAGTCCCCTGGGCTCTGTAAAGAATAGTGGGTGGATTCTTCATCCCAAATAAAGTGGTTTCT"
            "CAAGTGGTCCCAATTTTACAGTTCCTACCATCAGCTTCCCA, Y")
    ini.add_section("flanks")
    ini.add_comment("flanks",
        "The flanking sequences (e.g., primer sequences) of each marker. "
        "Specify two comma-separated values: left flank and right flank, "
        "in the same sequence orientation (strand).")
    if type != "non-str":
        ini.set("flanks", ";CSF1P0 = CCTGTGTCAGACCCTGTT, GTTGGAACACTGCCCTGG")
    if type != "str":
        ini.set("flanks",
            ";MitoFrag = ATTATTTATCGCACCTACGT, TGGCGGTATGCACTTTTAACAG")
    if aliases:
        ini.set("flanks", ";Amel = ACCCTGGTTATATCAACT, GTTTAAGCTCTGATGGTT")
    if type != "non-str":
        ini.add_section("prefix")
        ini.add_comment("prefix",
            "Specify all known prefix sequences of each STR marker, "
            "separated by commas. The prefix is the sequence between the left "
            "flank and the repeat and is omitted from allele names. The first "
            "sequence listed is used as the reference sequence for that "
            "marker when generating allele names. Deviations from the "
            "reference are expressed as variants. IUPAC notation for "
            "ambiguous bases (e.g., 'R' for 'A or G') is supported, except "
            "for the first sequence given. Lowercase letters represent "
            "optional bases.")
        ini.set("prefix", ";CSF1P0 = CTAAGTACTTC")
        ini.add_section("suffix")
        ini.add_comment("suffix",
            "Specify all known suffix sequences of each STR marker, "
            "separated by commas. The suffix is the sequence between the "
            "repeat and the right flank. The first sequence listed is used as "
            "the reference sequence for that marker when generating "
            "allele names.")
        ini.set("suffix",
            ";CSF1P0 = CTAATCTATCTATCTTCTATCTATGAAGGCAGTTACTGTTAATATCTTCATTTTA"
            "CAGGTAGGAAAACTGAGACACAGGGTGGTTAGCAACCTGCTAGTCCTTGGCAGACTCAG, CTAA"
            "TCTATCTATCTTCTATCTATGAAGGCAGTTACTGTTAATATCTTCATTTTACAGGTAGGAAAACT"
            "GAGACACAGGGTGGTTAGAAACCTGCTAGTCCTTGGCAGACTCAG, ATAATCTATCTATCTTCT"
            "ATCTATGAAGGCAGTTACTGTTAATATCTTCATTTTACAGGTAGGAAAACTGAGACACAGGGTGG"
            "TTAGCAACCTGCTAGTCCTTGGCAGACTCAG")
        ini.add_section("repeat")
        ini.add_comment("repeat",
            "Specify the repeat structure of each STR marker in space-"
            "separated triples of sequence, minimum number of repeats, "
            "and maximum number of repeats.")
        ini.set("repeat",
            ";CSF1P0 = CTAT 0 19 CTAC 0 1 TTAT 0 1 CAT 0 1 CTAT 0 19")
    if aliases or type != "str":
        ini.add_section("no_repeat")
        ini.add_comment("no_repeat",
            "Specify the reference sequence for each non-STR marker.")
        if type != "str":
            ini.set("no_repeat", ";MitoFrag = TCAATATTACAGGCGAACATACTTACTAAAGT"
            "GTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCC"
            "ACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTT"
            "AAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATT"
            "TTATCTTT")
        if aliases:
            ini.set("no_repeat", ";Amel = TCAGCTATGAGGTAATTTTTCTCTTTACTAATTTTG"
            "ACCATTGTTTGCGTTAACAATGCCCTGGGCTCTGTAAAGAATAGTGTGTTGATTCTTTATCCCAG"
            "ATGTTTCTCAAGTGGTCCTGATTTTACAGTTCCTACCACCAGCTTCCCA")
    ini.add_section("genome_position")
    ini.add_comment("genome_position",
        "Specify the chromosome number and position of the first base after "
        "the first flank of each marker. Optionally, you may specify the "
        "position of the last base of the fragment as well.%s" % (
            " Specify 'M' as the chromosome name for markers on mitochondrial "
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
            "position)\"." if type != "str" else ""))
    ini.add_comment("genome_position",
        "Using human genome build GRCh38%s." % (
            " and rCRS for human mtDNA" if type != "str" else ""))
    if type != "non-str":
        ini.set("genome_position", ";CSF1P0 = 5, 150076311, 150076487")
    if type != "str":
        ini.set("genome_position", ";MitoFrag = M, 173, 407")
    if aliases:
        ini.set("genome_position", ";Amel = X, 11296816, 11296965")
    if type == "str" or type == "full":
        ini.add_section("length_adjust")
        ini.add_comment("length_adjust",
            "To prevent discrepancies between traditional CE allele numbers "
            "and the CE number in FDSTools allele names, the CE allele number "
            "as calculated by FDSTools is based on the length of the sequence "
            "(prefix+repeat+suffix) minus the adjustment specified here.")
        ini.set("length_adjust", ";CSF1P0 = 0")
        ini.add_section("block_length")
        ini.add_comment("block_length",
            "Specify the repeat unit length of each STR marker. The default "
            "length is 4.")
        ini.set("block_length", ";CSF1P0 = 4")
    ini.add_section("max_expected_copies")
    ini.add_comment("max_expected_copies",
        "Specify the maximum expected number of copies (i.e., alleles) for "
        "each marker in a single reference sample (only used for "
        "allelefinder). The default is 2. Specify 1 here for haploid markers "
        "(i.e., those %son the Y chromosome)." % (
            "on mitochondrial DNA or " if type != "str" else ""))
    if type != "non-str":
        ini.set("max_expected_copies", ";CSF1P0 = 2")
    if type != "str":
        ini.set("max_expected_copies", ";MitoFrag = 1")
    if aliases:
        ini.set("max_expected_copies", ";Amel = 2")
    ini.add_section("expected_allele_length")
    ini.add_comment("expected_allele_length",
        "Specify one or two values for each marker. The first value gives the "
        "expected minimum length (in nucleotides, %sexcluding flanks) of the "
        "alleles and the second value (if given) specifies the maximum allele "
        "length expected for that marker (both inclusive). TSSV will filter "
        "sequences that have a length outside this range." %
        ("including prefix and suffix, " if type != "non-str" else ""))
    if type != "non-str":
        ini.set("expected_allele_length", ";CSF1P0 = 100")
    if type != "str":
        ini.set("expected_allele_length", ";MitoFrag = 150")
    if aliases:
        ini.set("expected_allele_length", ";Amel = 100")
    return ini
#make_empty_library_ini


def create_library(outfile, type, aliases=False):
    outfile.write(INI_COMMENT.fill("Lines beginning with a semicolon (;) are "
        "ignored by FDSTools.") + "\n\n")
    make_empty_library_ini(type, aliases).write(outfile)
#create_library


def add_arguments(parser):
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-t', '--type', metavar="TYPE", default="full",
        choices=("full", "str", "non-str"),
        help="the type of markers that this library file will be used for; "
             "'full' (the default) will create a library file with all "
             "possible sections, whereas 'str' or 'non-str' will only output "
             "the sections applicable to STR and non-STR markers, "
             "respectively")
    parser.add_argument('-a', '--aliases', action="store_true",
        help="if specified, the [aliases] section is included, which can be "
             "used to explicitly assign allele names to specific sequences of "
             "specific markers")
#add_arguments


def run(args):
    create_library(args.outfile, args.type, args.aliases)
#run
