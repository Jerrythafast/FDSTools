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

import argparse
import glob
import re
import sys

from .library import parse_library, get_builtin_library, BUILTIN_NAMES

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


def pos_int_arg(value):
    """Convert str to int, raise ArgumentTypeError if not positive."""
    if not value.isdigit() or not int(value):
        raise argparse.ArgumentTypeError("invalid positive int value: '%s'" % value)
    return int(value)
#pos_int_arg


def regex_arg(value):
    """Compile value into a regular expression."""
    try:
        return re.compile(value)
    except re.error as err:
        raise argparse.ArgumentTypeError(err)
#regex_arg


def library_arg(value):
    """Value is a library file name; parse and return the library."""
    try:
        if value == "-":
            return parse_library(sys.stdin)
        builtin = get_builtin_library(value)
        if builtin is not None:
            return builtin
        with open(value, "rt", encoding="UTF-8") as libfile:
            return parse_library(libfile)
    except Exception as err:
        raise argparse.ArgumentTypeError(err)
#library_arg


def comma_separated_arg(list_type, element_type):
    """
    Return a function that converts a comma-separated string to a
    list_type containing elements of element_type.
    """
    def parse_comma_separated_arg(value):
        try:
            return list_type(map(element_type, value.split(",")))
        except Exception as err:
            raise argparse.ArgumentTypeError(err)
    return parse_comma_separated_arg
#comma_separated_arg


def add_sequence_format_args(parser, *, default_format=None, force=False, require_library=False):
    """Add arguments for sequence formatting to the given parser."""
    group = parser.add_argument_group("sequence format options")
    if force:
        group.set_defaults(sequence_format=default_format)
    else:
        group.add_argument("-F", "--sequence-format", metavar="FORMAT",
            choices=("raw", "tssv", "allelename"), default=default_format,
            help="convert sequences to the specified format: one of %(choices)s (default: " + (
                 "no conversion" if default_format is None else default_format)
                 + ")")
    if require_library:
        parser.add_argument("library", metavar="LIBRARY", type=library_arg,
            help="library file with marker definitions; custom file or built-in: '%s'" %
                "', '".join(BUILTIN_NAMES))
    else:
        group.add_argument("-l", "--library", metavar="LIBRARY", type=library_arg,
            help="library file with marker definitions; custom file or built-in: '%s'" %
                "', '".join(BUILTIN_NAMES))
#add_sequence_format_args


def add_input_output_args(parser, *, single_in=False, batch_support=False, report_out=False):
    """Add arguments for opening sample files to the given parser."""
    # Input file options group.
    if not single_in:
        # Multiple input files: positionals.
        parser.add_argument("infiles", nargs="*", metavar="FILE", default=["-"],
            help="the sample data file(s) to process (default: read from stdin)")
    elif not batch_support:
        # Single input file and no batches: single positional.
        parser.add_argument("infile", nargs="?", metavar="IN", default="-",
            help="the sample data file to process (default: read from stdin)")
    else:
        # Single input file with batch support: single positional and -i
        # option for batches, which are mutually exclusive.
        mutex = parser.add_argument_group("input file options").add_mutually_exclusive_group()
        mutex.add_argument("infile", nargs="?", metavar="IN", default="-",
            help="single sample data file to process (default: read from stdin)")
        mutex.add_argument("-i", "--input", dest="infiles", nargs="+", metavar="IN",
            help="multiple sample data files to process (use with -o/--output)")

    # Output file options group.
    group = parser.add_argument_group("output file options")
    if batch_support and single_in:
        # Single input file with batch support: single positional and -o
        # option for batches, which are mutually exclusive.
        mutex = group.add_mutually_exclusive_group()
        mutex.add_argument("outfile", nargs="?", metavar="OUT", default=sys.stdout,
            help="the file to write the output to (default: write to stdout)")
        mutex.add_argument("-o", "--output", dest="outfiles", nargs="+", metavar="OUT",
            help="list of names of output files to match with input files "
                 "specified with -i/--input, or a format string to construct "
                 "file names from sample tags; e.g., the default value is "
                 "'\\1-%s.out', which expands to 'sampletag-%s.out'" %
                    ((parser.prog.rsplit(" ", 1)[-1],) * 2))
    elif single_in:
        # Single input file and no batch support: single positional.
        parser.add_argument("outfile", nargs="?", metavar="OUT",
            type=argparse.FileType("tw", encoding="UTF-8"),
            default=sys.stdout, help="the file to write the output to (default: write to stdout)")
    elif batch_support:
        # Multiple input files and batch support: use -o option.
        # (This is multi-in, multi-out).
        group.add_argument("-o", "--output", dest="outfiles", nargs="+", metavar="OUT",
            default=[sys.stdout],
            help="a single file name to write all output to (default: write "
                 "to stdout) OR a list of names of output files to match with "
                 "input files OR a format string to construct file names from "
                 "sample tags; e.g., the value '\\1-%s.out' expands to "
                 "'sampletag-%s.out'" % ((parser.prog.rsplit(" ", 1)[-1],) * 2))
    else:
        # Multiple input files and no batch support: use -o option.
        group.add_argument("-o", "--output", dest="outfile", metavar="FILE",
            type=argparse.FileType("tw", encoding="UTF-8"), default=sys.stdout,
            help="file to write output to (default: write to stdout)")
    if report_out:
        group.add_argument("-R", "--report", metavar="FILE",
            type=argparse.FileType("tw", encoding="UTF-8"),
            default=sys.stderr, help="file to write a report to (default: write to stderr)")

    # Sample tag parsing options group.
    if not single_in or batch_support:
        group = parser.add_argument_group("sample tag parsing options",
            "for details about REGEX syntax and capturing groups, check "
            "https://docs.python.org/howto/regex")
        group.add_argument("-e", "--tag-expr", metavar="REGEX", type=regex_arg,
            default=DEF_TAG_EXPR,
            help="regular expression that captures (using one or more "
                 "capturing groups) the sample tags from the file names; by "
                 "default, the entire file name except for its extension (if "
                 "any) is captured")
        group.add_argument("-f", "--tag-format", metavar="EXPR", default=DEF_TAG_FORMAT,
            help="format of the sample tags produced; a capturing group "
                 "reference like '\\n' refers to the n-th capturing group in "
                 "the regular expression specified with -e/--tag-expr (the "
                 "default of '\\1' simply uses the first capturing group); "
                 "with a single sample, you can enter the sample tag here "
                 "explicitly")
#add_input_output_args


def add_allele_detection_args(parser):
    """Add arguments for specifying known alleles of a sample."""
    group = parser.add_argument_group("allele detection options")
    group.add_argument("-a", "--allelelist", metavar="ALLELEFILE",
        type=argparse.FileType("tr", encoding="UTF-8"),
        help="file containing a list of the true alleles of each sample "
             "(e.g., obtained from allelefinder)")
    group.add_argument("-c", "--annotation-column", metavar="COLNAME",
        help="name of a column in the sample files, which contains a value "
             "beginning with 'ALLELE' for the true alleles of the sample")
#add_allele_detection_args


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


def get_input_output_files(args, *, single_in=False, batch_support=False):
    """
    If no input has been specified, return False.
    If single_in and not batch_support, return (infile, outfile).
    If not single_in and not batch_support, return ({tag: infiles}, out).
    If single_in and batch_support, return generator of (tag, [ins], out).
    If not single and batch_support, raise NotImplementedError.

    The input files are returned as strings, the outputs are opened.
    """
    if single_in and not batch_support:
        # One infile, one outfile.  Return 2-tuple (infile, outfile).
        if args.infile == "-" and sys.stdin.isatty():
            return False  # No input specified.
        return args.infile, args.outfile


    if not single_in and not batch_support:
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


    if single_in and batch_support:
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

        tags = [get_tag(infile, args.tag_expr, args.tag_format) for infile in infiles]

        if len(outfiles) == 1:
            outfile = sys.stdout if outfiles[0] == "-" else outfiles[0]

            if outfile == sys.stdout and len(set(tags)) == 1:
                # Write output of single sample to stdout.
                return ((tag, infiles, outfile) for tag in set(tags))

            # Write output of each sample to its own outfile.
            if outfile == sys.stdout:
                outfile = DEF_OUTFILE_FORMAT
            return ((tag,
                    [infiles[i] for i in range(len(tags)) if tags[i] == tag],
                    open(outfile.replace("\\1", tag).replace("\\2", args.tool), "wt",
                        encoding="UTF-8"))
                    for tag in set(tags))

        # Link each output file to each input file.
        # Treating files with the same sample tag as separate samples.
        return ((tag, [infiles[i]], open(outfiles[i], "wt", encoding="UTF-8"))
            for i, tag in enumerate(tags))

    if not single_in and batch_support:
        # N infiles, one or N outfiles.
        # If one outfile, return ({tag: [infiles]}, outfile).
        # If N outfiles, return generator of (tag, [infiles], outfile).
        raise NotImplementedError("Multi-input with optional multi-output not supported yet.")
#get_input_output_files
