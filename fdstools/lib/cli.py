#!/usr/bin/env python3

#
# Copyright (C) 2026 Netherlands Forensic Institute
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

from argparse import ArgumentTypeError
from glob import iglob
from itertools import chain

from .io import std_or_open
from .library import parse_library, get_builtin_library, BUILTIN_NAMES

# Default regular expression to capture sample tags in file names.
# This is the default of the -e command line option.
DEF_TAG_EXPR = r"^(.*?)(?:\.[^.]+)?$"

# Default formatting template to write sample tags.
# This is the default of the -f command line option.
DEF_TAG_FORMAT = r"\1"

# Default formatting template to construct output file names for batch
# processing.  \1 and \2 refer to sample tag and tool name.
# This is the default for the -o command line option with batch support.
DEF_OUTFILE_FORMAT = r"\1-\2.out"


def pos_int_arg(value):
    """Convert str to int, raise ArgumentTypeError if not positive."""
    if not value.isdigit() or not int(value):
        raise ArgumentTypeError("invalid positive int value: '%s'" % value)
    return int(value)
#pos_int_arg


def regex_arg(value):
    """Compile value into a regular expression."""
    try:
        return re.compile(value)
    except re.error as err:
        raise ArgumentTypeError(err)
#regex_arg


def library_arg(value):
    """Value is a library file name; parse and return the library."""
    try:
        builtin = get_builtin_library(value)
        if builtin is not None:
            return builtin
        with std_or_open(value) as libfile:
            return parse_library(libfile)
    except Exception as err:
        raise ArgumentTypeError(err)
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
            raise ArgumentTypeError(err)
    return parse_comma_separated_arg
#comma_separated_arg


def add_sequence_format_args(parser, *, default_format=None, force=False,
                             require_library=False, ce_format=False):
    """Add arguments for sequence formatting to the given parser."""
    group = parser.add_argument_group("sequence format options")
    if force:
        group.set_defaults(sequence_format=default_format)
    else:
        if ce_format:
            choices = ("raw", "tssv", "allelename", "ce")
        else:
            choices = ("raw", "tssv", "allelename")
        group.add_argument("-F", "--sequence-format", metavar="FORMAT",
            choices=choices, default=default_format,
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


def add_sample_tag_parsing_args(parser):
    group = parser.add_argument_group("sample tag parsing options",
        "for details about REGEX syntax and capturing groups, check "
        "https://docs.python.org/howto/regex")
    group.add_argument("-e", "--tag-expr", metavar="REGEX", type=regex_arg, default=DEF_TAG_EXPR,
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
#add_sample_tag_parsing_args


def add_batch_input_output_args(parser):
    """
    Add arguments for opening many sample files to the given parser,
    such that it can perform batch processing of many samples.

    The samples can be iterated over using iterate_samples().
    """
    # Input file options group: single positional and -i
    # option for batches, which are mutually exclusive.
    mutex = parser.add_argument_group("input file options").add_mutually_exclusive_group()
    mutex.add_argument("infile", nargs="?", metavar="IN", default="-",
        help="single sample data file to process (default: read from stdin)")
    mutex.add_argument("-i", "--input", dest="infiles", nargs="+", metavar="IN",
        help="multiple sample data files to process (use with -o/--output)")

    # Output file options group: single positional and -o
    # option for batches, which are mutually exclusive.
    group = parser.add_argument_group("output file options")
    mutex = group.add_mutually_exclusive_group()
    mutex.add_argument("outfile", nargs="?", metavar="OUT",
        help="the file to write the output to (default: write to stdout)")
    mutex.add_argument("-o", "--output", dest="outfiles", nargs="+", metavar="OUT",
        help="list of names of output files to match with input files "
             "specified with -i/--input, or a format string to construct "
             "file names from sample tags; e.g., the default value is "
             "'\\1-%s.out', which expands to 'sampletag-%s.out'" %
                ((parser.prog.rsplit(" ", 1)[-1],) * 2))

    # Sample tag parsing options group.
    add_sample_tag_parsing_args(parser)
#add_batch_input_output_args


def add_bulk_input_output_args(parser, *, output_group=False):
    """
    Add arguments for opening many sample files to the given parser.

    If output_group is True, add the -o/--output option to a newly
    created 'output file options' group, and return the created
    group to allow the caller to add more options to it.

    The samples can be iterated over using iterate_input_samples().
    """
    parser.add_argument("infiles", nargs="*", metavar="FILE", default=["-"],
        help="the sample data file(s) to process (default: read from stdin)")

    group = parser.add_argument_group("output file options") if output_group else parser
    group.add_argument("-o", "--output", dest="outfile", metavar="FILE", default="-",
        help="file to write output to (default: write to stdout)")

    add_sample_tag_parsing_args(parser)
    return group
#add_bulk_input_output_args


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
    for file in iglob(pathname):
        success = True
        yield file
    if not success:
        yield pathname
#glob_path


def iterate_input_samples(infiles, tag_expr, tag_format):
    """
    Yield tuples of (tag, infile), raise ValueError on duplicate tag.
    Raise ValueError when any sample tag is encountered a second time.
    The files will be opened and closed automatically during iteration.
    """
    seen_tags = set()
    for infile in chain.from_iterable(map(glob_path, infiles)):
        tag = get_tag(infile, tag_expr, tag_format)
        if tag in seen_tags:
            raise ValueError(f"multiple input files for sample '{tag}' specified")
        seen_tags.add(tag)
        with std_or_open(infile) as infile:
            yield tag, infile
#iterate_input_samples


def iterate_samples(args, *, open_files=True):
    """
    Return a generator of (tag, infile, outfile).

    If a duplicate tag exists, ValueError is raised prior to
    opening any files.

    Each yielded tuple should cause a separate run of the tool.

    If open_files is True, the files will be opened and closed
    automatically during iteration.
    """
    infiles = list(chain.from_iterable(map(glob_path, args.infiles))) if "infiles" \
              in args and args.infiles is not None else [args.infile]
    outfiles = args.outfiles if "outfiles" in args \
               and args.outfiles is not None else [args.outfile]

    tags = [get_tag(infile, args.tag_expr, args.tag_format) for infile in infiles]
    if len(tags) != len(set(tags)):
        duplicate_tags = "', '".join(tag for tag in set(tags) if tags.count(tag) > 1)
        raise ValueError(f"multiple input files for sample '{duplicate_tags}' specified")

    if len(outfiles) == 1:
        outfile = outfiles[0]
        if outfile is None:
            # No output specified: write output of single sample to stdout.
            outfile = "-" if len(infiles) == 1 else DEF_OUTFILE_FORMAT
        for tag, infile in zip(tags, infiles):
            outfilename = outfile.replace(r"\1", tag).replace(r"\2", args.tool)
            if open_files:
                with std_or_open(infile) as infile, \
                     std_or_open(outfilename, "wt") as outstream:
                    yield tag, infile, outstream
            else:
                yield tag, infile, outfilename
    elif len(outfiles) == len(infiles):
        for tag, infile, outfile in zip(tags, infiles, outfiles):
            if open_files:
                with std_or_open(infile) as infile, \
                     std_or_open(outfile, "wt") as outfile:
                    yield tag, infile, outfile
            else:
                yield tag, infile, outfile
    else:
        raise ValueError(
            f"number of input files ({len(infiles)}) is not equal to number of output "
            f"files ({len(outfiles)})")
#iterate_samples
