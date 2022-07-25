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

"""
Automatically run complete, predefined analysis pipelines.  Recommended
starting point for new users.

This tool runs one of three default analysis pipelines automatically,
given a configuration file with tool options and input/output file
names.  The three available analysis options are 'reference-sample',
analysing a single reference sample with TSSV and Stuttermark;
'reference-database', analysing a collection of reference samples with
BGEstimate and Stuttermodel; and 'case-sample', analysing a single case
sample with TSSV, BGPredict, BGMerge, BGCorrect, and Samplestats.  All
results are visualised in interactive graphical reports for presentation
and further interpretation.

This tool takes a single mandatory argument: the name of an INI
configuration file that contains the analysis settings to use.  An easy
way to obtain such an INI file with default values for all settings, is
to run 'fdstools pipeline your-filename.ini --analysis case-sample'.
This will create the file 'your-filename.ini' and fill it with default
values for the given analysis type (in this example: case-sample
analysis).

All settings in the configuration file correspond to options of various
tools in FDSTools.  Please refer to the tool-specific help for a full
description of each tool.  Type 'fdstools -h TOOL' to get help with the
given TOOL.
"""
import os
import pkgutil
import re
import subprocess
import sys
import tempfile
import threading

import fdstools.tools

from configparser import RawConfigParser, NoSectionError, NoOptionError
from functools import reduce

from ..lib.cli import DEF_TAG_EXPR, DEF_TAG_FORMAT, get_tag, regex_arg, glob_path
from ..lib.io import print_db
from ..lib.library import INI_COMMENT


__version__ = "1.1.1"


# Pattern to split a quoted string.
PAT_SPLIT_QUOTED = re.compile(r""""((?:\\"|[^"])*)"|'((?:\\'|[^'])*)'|(\S+)""")

# Pattern that matches a long argparse argument name.
PAT_ARGNAME = re.compile("^(?:--)?([a-z0-9]+(?:-[a-z0-9]+)*)$")

# Pattern that matches an argument name in documentation ('-a/--long').
PAT_ARGNAME_IN_DOC = re.compile("-[a-zA-Z]/--([a-z0-9]+(?:-[a-z0-9]+)*)")

# Prefix of the tools subpackage of FDSTools.
PACKAGE_PREFIX = fdstools.tools.__name__ + "."

# Name of this tool.
NAME = __name__[len(PACKAGE_PREFIX):]

# Argument names that will not be written in generated INI files.
HIDDEN_ARGS = ("infile", "infile2", "infiles", "outfile", "outfiles", "raw-outfile",
    "library", "library2", "allelelist", "report", "tag-expr", "tag-format", "uncall-alleles",
    "sequence-format", "type", "title", "stuttermodel", "seqs", "profiles", "combine-strands",
    "annotation-column", "marker-column", "allele-column", "output-column", "reverse-complement")

# Tools that are run for each type of analysis pipeline.
ANALYSIS_TOOLS = {
    "reference-sample":
        ("pipeline", "tssv", "stuttermark", "seqconvert", "vis"),
    "reference-database":
        ("pipeline", "allelefinder", "bgestimate", "bghomraw", "stuttermodel",
        "seqconvert", "vis"),
    "case-sample":
        ("pipeline", "tssv", "bgpredict", "bgmerge", "bgcorrect",
        "samplestats", "seqconvert", "vis")
}

# Pipeline arguments for each analysis pipeline.
ANALYSIS_ARGS = {
    "reference-sample":
        ("analysis", "in-library", "in-sample-raw", "tag-expr", "tag-format"),
    "reference-database":
        ("analysis", "in-library", "in-samples", "in-allelelist", "prefix",
        "tag-expr", "tag-format", "combine-strands"),
    "case-sample":
        ("analysis", "in-library", "in-sample-raw", "in-stuttermodel",
        "in-bgprofiles", "store-predictions", "tag-expr", "tag-format", "combine-strands")
}


class ArgumentCollector:
    class SubCollector:
        def __init__(self, collector, tool):
            self.collector = collector
            self.tool = tool
            self.prog = tool
        #__init__

        def add_argument(self, *args, **kwargs):
            self.collector.add_argument(self.tool, *args, **kwargs)
        #add_argument

        def add_argument_group(self, *args, **kwargs):
            return self
        #add_argument_group

        def add_mutually_exclusive_group(self, *args, **kwargs):
            return self
        #add_mutually_exclusive_group

        def set_defaults(self, *args, **kwargs):
            pass
        #set_defaults
    #SubCollector


    def __init__(self):
        self.arglists = {}
    #__init__

    def add_tool(self, tool):
        self.arglists[tool] = []
        return self.SubCollector(self, tool)
    #add_tool

    def add_argument(self, tool, *args, **kwargs):
        # Determine argument name.
        if "dest" in kwargs:
            name = kwargs["dest"].replace("_", "-")
        else:
            for arg in args:
                match = PAT_ARGNAME.match(arg)
                if match is not None:
                    name = match.group(1)
                    break
            else:
                name = args[0]
        if args[0].startswith("-"):
            kwargs["optname"] = args[0]

        # Replace special default values.
        if "optname" in kwargs and "nargs" in kwargs and kwargs["nargs"] == "?":
            # Value can be False (option unspecified),
            # True (option specified),
            # or an actual value (option specified with argument).
            kwargs["default"] = False
        elif "default" in kwargs:
            if kwargs["default"] is None or kwargs["default"] == "":
                del kwargs["default"]
            elif kwargs["default"] in (sys.stdin, sys.stdout):
                kwargs["default"] = "-"
            elif kwargs["default"] is sys.stderr:
                kwargs["default"] = "/dev/stderr"
            elif kwargs["default"] is sys.maxsize:
                kwargs["default"] = "<MAX>"
            elif isinstance(kwargs["default"], list):
                kwargs["default"] = " ".join(kwargs["default"])
        elif "action" in kwargs and kwargs["action"] == "store_true":
            kwargs["default"] = False
        elif "action" in kwargs and kwargs["action"] == "store_false":
            kwargs["default"] = True
        self.arglists[tool].append((name, kwargs))
    #add_argument

    def get_argument_lists(self):
        return self.arglists
    #get_argument_lists
#ArgumentCollector


def split_quoted_string(text):
    return ["".join((
        y[0].replace("\\\"", "\""),
        y[1].replace("\\'", "'"),
        y[2])) for y in PAT_SPLIT_QUOTED.findall(text)]
#split_quoted_string


def format_help(text, format_args={}):
    text = PAT_ARGNAME_IN_DOC.sub("'\\1'",
        text.replace(
            "-i/--input", "'infiles'").replace(
            "-o/--output", "'outfiles'") % format_args)
    return INI_COMMENT.fill(text[:1].upper() + text[1:])
#format_help


def get_tools(hide_pipeline):
    tools = {}
    for importer, name, ispkg in pkgutil.iter_modules(
            path=fdstools.tools.__path__):
        if hide_pipeline and name == NAME:
            continue
        try:
            tools[name] = importer.find_module(PACKAGE_PREFIX + name).load_module(
                PACKAGE_PREFIX + name)
        except Exception as error:
            sys.stderr.write("FDSTools failed to load '%s': %s\n" % (name, error))
            continue
    return tools
#get_tools


def get_arguments(tools):
    collector = ArgumentCollector()
    for name in tools:
        tools[name].add_arguments(collector.add_tool(name))
    return collector.get_argument_lists()
#get_arguments


def read_ini(infile):
    config = RawConfigParser()
    config.optionxform = str
    config.readfp(infile)
    return config
#read_ini


def parse_bool_arg(section, option, value):
    if value.lower() in ("false", "no", "off", "0", ""):
        return False
    if value.lower() in ("true", "yes", "on", "1"):
        return True
    raise ValueError(
            "Invalid boolean value '%s' for option '%s' of %s" %
            (value, option, section))
#parse_bool_arg


def get_argv(toolname, arg_defs, config):
    """Get argument list for the given tool."""
    arglist = (["fdstools", toolname], [])
    if debug:
        arglist[0].append("-d")
    for arg in arg_defs[toolname]:
        try:
            value = config.get(toolname, arg[0]).strip()
            if value == "":
                raise ValueError
        except (NoSectionError, NoOptionError, ValueError):
            if "optname" not in arg[1] and "default" in arg[1]:
                value = arg[1]["default"]
            else:
                continue
        if "action" in arg[1] and arg[1]["action"] == "store_true":
            if parse_bool_arg(toolname, arg[0], value):
                arglist[0].append(arg[1]["optname"])
            continue
        if "action" in arg[1] and arg[1]["action"] == "store_false":
            if not parse_bool_arg(toolname, arg[0], value):
                arglist[0].append(arg[1]["optname"])
            continue
        if "default" in arg[1] and arg[1]["default"] == value:
            if "optname" in arg[1]:
                continue
            if value == "<MAX>":
                value = str(sys.maxsize)
        if "nargs" in arg[1] and arg[1]["nargs"] not in (1, "?"):
            value = split_quoted_string(value)
        elif "optname" in arg[1] and "nargs" in arg[1] and arg[1]["nargs"] == "?":
            # Value can be False (option unspecified),
            # True (option specified),
            # or an actual value (option specified with argument).
            try:
                if value in ("0", "1"):
                    # Check if 0 or 1 is a valid argument value.
                    arg[1]["type"](value)
                    value = [value]  # Specify with given argument.
                else:
                    raise ValueError  # Proceed to next check.
            except:
                try:
                    # Try to parse the value as a boolean (no argument).
                    value = parse_bool_arg(toolname, arg[0], value)
                    if not value:
                        continue  # Value False, don't specify option.
                    value = None  # Value True, specify without arg.
                except ValueError:
                    value = [value]  # Specify with given argument.
        else:
            value = [value]
        if "optname" in arg[1]:
            if value is not None and len(value) == 1 and value[0].startswith("-"):
                # Fix Argparse's "looks like an option" oddity.
                arglist[0].append("%s=%s" % (arg[1]["optname"], value[0]))
            else:
                arglist[0].append(arg[1]["optname"])
                if value is not None:
                    arglist[0].extend(value)
        elif value is not None:
            arglist[1].extend(value)
    print_db(("%r" % (arglist[0]+arglist[1]))[:200], debug)
    return arglist[0] + arglist[1]
#get_argv


def ini_require_option(config, section, option):
    """
    Return the value of option in section of config or raise ValueError.
    """
    try:
        r = config.get(section, option).strip()
        if not r:
            raise ValueError(
                "The value of the '%s' option cannot be empty" % option)
        return r
    except (NoSectionError, NoOptionError):
        raise ValueError(
            "The configuration file should have a section named [%s] with "
            "an option named '%s'" % (section, option))
#ini_require_option


def ini_try_get_option(config, section, option, default=None):
    """
    Return the value of option in section of config or default.
    """
    try:
        return config.get(section, option).strip() or default
    except (NoSectionError, NoOptionError):
        return default
#ini_try_get_option


def tee(instream, *outstreams):
    """
    Copy data from instream to all outstreams until EOF, close streams.
    """
    for line in instream:
        for outstream in outstreams:
            outstream.write(line)
    instream.close()
    for outstream in outstreams:
        outstream.close()
#tee


def write_to_stream(stream, data):
    """
    Write data to the given stream, then close the stream.
    """
    stream.write(data)
    stream.close()
#write_to_stream


def run_ref_sample_analysis(arg_defs, config):
    # Read input/output file configuration.
    in_library = ini_require_option(config, NAME, "in-library")
    if in_library == "-":
        raise ValueError("The library file cannot be named '-'")
    in_sample_raw = ini_require_option(config, NAME, "in-sample-raw")
    if in_sample_raw == "-":
        raise ValueError("The raw sample input file cannot be named '-'")
    tag = get_tag(
        in_sample_raw,
        regex_arg(ini_try_get_option(config, NAME, "tag-expr", DEF_TAG_EXPR)),
        ini_try_get_option(config, NAME, "tag-format", DEF_TAG_FORMAT))

    # Overwrite any explicit I/O configuration.
    if not config.has_section("tssv"):
        config.add_section("tssv")
    config.set("tssv", "library", in_library)
    config.set("tssv", "infile", in_sample_raw)
    config.set("tssv", "outfile", "-")
    config.set("tssv", "report", tag + "-stats.txt")
    if not config.has_section("stuttermark"):
        config.add_section("stuttermark")
    config.set("stuttermark", "infile", "-")
    config.remove_option("stuttermark", "infiles")
    config.set("stuttermark", "outfile", tag + ".csv")
    config.remove_option("stuttermark", "outfiles")
    config.set("stuttermark", "library", in_library)
    if not config.has_section("seqconvert"):
        config.add_section("seqconvert")
    if not config.has_option("seqconvert", "sequence-format"):
        config.set("seqconvert", "sequence-format", "allelename")
    config.set("seqconvert", "infile", "-")
    config.remove_option("seqconvert", "infiles")
    config.set("seqconvert", "outfile", "-")
    config.remove_option("seqconvert", "outfiles")
    config.set("seqconvert", "library", in_library)
    config.set("seqconvert", "marker-column", "marker")
    config.set("seqconvert", "allele-column", "sequence")
    if not config.has_section("vis"):
        config.add_section("vis")
    config.set("vis", "type", "sample")
    config.set("vis", "infile", "-")
    config.set("vis", "outfile", tag + ".html")
    if not config.has_option("vis", "title"):
        config.set("vis", "title", tag)

    # Start with TSSV.
    p_tssv = subprocess.Popen(get_argv("tssv", arg_defs, config),
        stdout=subprocess.PIPE)

    # Run Stuttermark on TSSV output.
    p_stuttermark = subprocess.Popen(
        get_argv("stuttermark", arg_defs, config), stdin=subprocess.PIPE)

    # Convert sequences to allele names with Seqconvert for Vis.
    p_seqconvert = subprocess.Popen(
        get_argv("seqconvert", arg_defs, config),
        stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    # Visualise TSSV output with Vis.
    p_vis = subprocess.Popen(get_argv("vis", arg_defs, config),
        stdin=p_seqconvert.stdout)
    p_seqconvert.stdout.close()

    # Send TSSV output to Samplestats and Seqconvert/Vis.
    tee(p_tssv.stdout, p_stuttermark.stdin, p_seqconvert.stdin)

    # Wait for the tools to finish.
    p_tssv.wait()
    p_stuttermark.wait()
    p_seqconvert.wait()
    p_vis.wait()
#run_ref_sample_analysis


def run_ref_database_analysis(arg_defs, config):
    # Read input/output file configuration.
    in_library = ini_require_option(config, NAME, "in-library")
    if in_library == "-":
        raise ValueError("The library file cannot be named '-'")
    in_samples = [
        "'" + x.replace("'", "\'") + "'"
        for x in split_quoted_string(ini_require_option(
            config, NAME, "in-samples")) for x in glob_path(x)]
    if "'-'" in in_samples:
        raise ValueError("The sample input files cannot be named '-'")
    in_samples = " ".join(x for x in in_samples if x != "''")
    if not in_samples:
        raise ValueError("No sample input files given")
    tag_expr = ini_try_get_option(config, NAME, "tag-expr", DEF_TAG_EXPR)
    tag_format = ini_try_get_option(config, NAME, "tag-format", DEF_TAG_FORMAT)
    prefix = ini_try_get_option(config, NAME, "prefix", "")
    if prefix and prefix[-1] not in "-_./\\":
        prefix += "-"
    allelefile = ini_try_get_option(config, NAME, "in-allelelist")
    run_allelefinder = False
    if allelefile is None:
        allelefile =  prefix + "allelelist.txt"
        run_allelefinder = True
    combine_strands = ini_try_get_option(config, NAME, "combine-strands", "")

    # Overwrite any explicit I/O configuration.
    if run_allelefinder:
        if not config.has_section("allelefinder"):
            config.add_section("allelefinder")
        config.set("allelefinder", "infiles", in_samples)
        config.set("allelefinder", "outfile", allelefile)
        config.set("allelefinder", "report", prefix + "allelereport.txt")
        config.set("allelefinder", "library", in_library)
        config.set("allelefinder", "tag-expr", tag_expr)
        config.set("allelefinder", "tag-format", tag_format)
        if not config.has_option("allelefinder", "stuttermark-column"):
            config.set("allelefinder", "stuttermark-column", "annotation")
    if not config.has_section("bgestimate"):
        config.add_section("bgestimate")
    config.set("bgestimate", "allelelist", allelefile)
    config.set("bgestimate", "infiles", in_samples)
    config.set("bgestimate", "outfile", prefix + "bgprofiles.txt")
    config.set("bgestimate", "library", in_library)
    config.set("bgestimate", "tag-expr", tag_expr)
    config.set("bgestimate", "tag-format", tag_format)
    config.set("bgestimate", "combine-strands", combine_strands)
    if not config.has_section("stuttermodel"):
        config.add_section("stuttermodel")
    config.set("stuttermodel", "allelelist", allelefile)
    config.set("stuttermodel", "infiles", in_samples)
    config.set("stuttermodel", "outfile", prefix + "stuttermodel.txt")
    config.set("stuttermodel", "library", in_library)
    config.set("stuttermodel", "tag-expr", tag_expr)
    config.set("stuttermodel", "tag-format", tag_format)
    config.set("stuttermodel", "raw-outfile", prefix + "stuttermodel-raw.txt")
    config.set("stuttermodel", "combine-strands", combine_strands)
    if not config.has_section("bghomraw"):
        config.add_section("bghomraw")
    config.set("bghomraw", "allelelist", allelefile)
    config.set("bghomraw", "infiles", in_samples)
    config.set("bghomraw", "outfile", prefix + "bgprofiles-raw.txt")
    config.set("bghomraw", "library", in_library)
    config.set("bghomraw", "tag-expr", tag_expr)
    config.set("bghomraw", "tag-format", tag_format)
    config.set("bghomraw", "combine-strands", combine_strands)
    if not config.has_option("bghomraw", "sequence-format"):
        config.set("bghomraw", "sequence-format", "allelename")
    if not config.has_section("vis"):
        config.add_section("vis")
    if not config.has_section("seqconvert"):
        config.add_section("seqconvert")
    if not config.has_option("seqconvert", "sequence-format"):
        config.set("seqconvert", "sequence-format", "allelename")
    config.remove_option("seqconvert", "infiles")
    config.set("seqconvert", "outfile", "-")
    config.remove_option("seqconvert", "outfiles")
    config.set("seqconvert", "library", in_library)
    config.set("seqconvert", "marker-column", "marker")
    config.set("seqconvert", "allele-column", "sequence")

    if run_allelefinder:
        # Start with Allelefinder.
        p_allelefinder = subprocess.Popen(
            get_argv("allelefinder", arg_defs, config))
        p_allelefinder.wait()

    p_bgestimate = subprocess.Popen(get_argv("bgestimate", arg_defs, config))

    p_stuttermodel = subprocess.Popen(
        get_argv("stuttermodel", arg_defs, config))

    p_bghomraw = subprocess.Popen(get_argv("bghomraw", arg_defs, config))

    config.set("vis", "type", "allele")
    config.set("vis", "infile", allelefile)
    config.remove_option("vis", "infile2")
    config.set("vis", "outfile", prefix + "allelegraph.html")
    p_allelevis = subprocess.Popen(get_argv("vis", arg_defs, config))

    # Wait for the tools to finish (in order of expected running time).
    p_allelevis.wait()
    p_bghomraw.wait()
    p_bgestimate.wait()
    p_stuttermodel.wait()

    # Make background profiles visualisation.
    vis_custom_title = config.has_option("vis", "title")
    config.set("seqconvert", "infile", prefix + "bgprofiles.txt")
    p_seqconvert1 = subprocess.Popen(
        get_argv("seqconvert", arg_defs, config), stdout=subprocess.PIPE)
    config.set("seqconvert", "infile", "-")
    config.set("seqconvert", "allele-column", "allele")
    p_seqconvert2 = subprocess.Popen(
        get_argv("seqconvert", arg_defs, config),
        stdin=p_seqconvert1.stdout,
        stdout=subprocess.PIPE)
    p_seqconvert1.stdout.close()
    config.set("vis", "type", "profile")
    config.set("vis", "infile", "-")
    config.set("vis", "infile2", prefix + "bgprofiles-raw.txt")
    config.set("vis", "outfile", prefix + "bgprofiles.html")
    if not vis_custom_title:
        config.set("vis", "title", prefix + "bgprofiles")
    p_profilevis = subprocess.Popen(
        get_argv("vis", arg_defs, config), stdin=p_seqconvert2.stdout)
    p_seqconvert2.stdout.close()

    # Make Stuttermodel visualisation.
    config.set("vis", "type", "stuttermodel")
    config.set("vis", "infile", prefix + "stuttermodel.txt")
    config.set("vis", "infile2", prefix + "stuttermodel-raw.txt")
    config.set("vis", "outfile", prefix + "stuttermodel.html")
    if not vis_custom_title:
        config.remove_option("vis", "title")
    p_stuttermodelvis = subprocess.Popen(get_argv("vis", arg_defs, config))

    # Make raw background noise visualisation.
    config.set("vis", "type", "bgraw")
    config.set("vis", "infile", prefix + "bgprofiles-raw.txt")
    config.remove_option("vis", "infile2")
    config.set("vis", "outfile", prefix + "bgprofiles-raw.html")
    p_bgrawvis = subprocess.Popen(get_argv("vis", arg_defs, config))

    # Wait for the visualisations to be made.
    p_stuttermodelvis.wait()
    p_bgrawvis.wait()
    p_seqconvert1.wait()
    p_seqconvert2.wait()
    p_profilevis.wait()
#run_ref_database_analysis


def run_case_sample_analysis(arg_defs, config):
    # Read input/output file configuration.
    in_library = ini_require_option(config, NAME, "in-library")
    if in_library == "-":
        raise ValueError("The library file cannot be named '-'")
    in_stuttermodel = ini_try_get_option(config, NAME, "in-stuttermodel")
    in_bgprofiles = ini_try_get_option(config, NAME, "in-bgprofiles")
    in_sample_raw = ini_require_option(config, NAME, "in-sample-raw")
    if in_sample_raw == "-":
        raise ValueError("The raw sample input file cannot be named '-'")
    tag = get_tag(
        in_sample_raw,
        regex_arg(ini_try_get_option(config, NAME, "tag-expr", DEF_TAG_EXPR)),
        ini_try_get_option(config, NAME, "tag-format", DEF_TAG_FORMAT))
    store_predictions = parse_bool_arg(
        NAME, "store-predictions", ini_try_get_option(
            config, NAME, "store-predictions", "False"))
    out_bgpredict = tag + "-bgpredict.txt" if store_predictions else "-"
    out_bgmerge = tag + "-bgmerge.txt" if store_predictions else "-"
    combine_strands = ini_try_get_option(config, NAME, "combine-strands", "")


    # Overwrite any explicit I/O configuration.
    if not config.has_section("tssv"):
        config.add_section("tssv")
    config.set("tssv", "library", in_library)
    config.set("tssv", "infile", in_sample_raw)
    config.set("tssv", "outfile", "-")
    config.set("tssv", "report", tag + "-stats.txt")
    if in_stuttermodel is not None:
        if not config.has_section("bgpredict"):
            config.add_section("bgpredict")
        config.set("bgpredict", "stuttermodel", in_stuttermodel)
        config.set("bgpredict", "seqs", "-")
        config.set("bgpredict", "outfile", out_bgpredict)
        config.set("bgpredict", "library", in_library)
        config.set("bgpredict", "combine-strands", combine_strands)
        if in_bgprofiles is not None:
            if not config.has_section("bgmerge"):
                config.add_section("bgmerge")
            config.set("bgmerge", "infiles", " ".join(
                "'" + x.replace("'", "\'") + "'" for x in
                    (in_bgprofiles, out_bgpredict)))
            config.set("bgmerge", "outfile", out_bgmerge)
            config.set("bgmerge", "library", in_library)
    if in_stuttermodel is not None or in_bgprofiles is not None:
        if not config.has_section("bgcorrect"):
            config.add_section("bgcorrect")
        config.set("bgcorrect", "profiles",
            out_bgmerge if in_bgprofiles is not None and
                in_stuttermodel is not None else
            in_bgprofiles if in_bgprofiles is not None else
            out_bgpredict)
        config.set("bgcorrect", "infile", "-")
        config.remove_option("bgcorrect", "infiles")
        config.set("bgcorrect", "outfile", "-")
        config.remove_option("bgcorrect", "outfiles")
        config.set("bgcorrect", "library", in_library)
        config.set("bgcorrect", "combine-strands", combine_strands)
    if not config.has_section("seqconvert"):
        config.add_section("seqconvert")
    if not config.has_option("seqconvert", "sequence-format"):
        config.set("seqconvert", "sequence-format", "allelename")
    config.set("seqconvert", "infile", "-")
    config.remove_option("seqconvert", "infiles")
    config.set("seqconvert", "outfile", "-")
    config.remove_option("seqconvert", "outfiles")
    config.set("seqconvert", "library", in_library)
    config.set("seqconvert", "marker-column", "marker")
    config.set("seqconvert", "allele-column", "sequence")
    if not config.has_section("samplestats"):
        config.add_section("samplestats")
    config.set("samplestats", "infile", "-")
    config.remove_option("samplestats", "infiles")
    config.set("samplestats", "outfile", tag + ".csv")
    config.remove_option("samplestats", "outfiles")
    if not config.has_section("vis"):
        config.add_section("vis")
    config.set("vis", "type", "sample")
    config.set("vis", "infile", "-")
    config.set("vis", "outfile", tag + ".html")
    if not config.has_option("vis", "title"):
        config.set("vis", "title", tag)

    # Start with TSSV.
    p_tssv = subprocess.Popen(get_argv("tssv", arg_defs, config),
        stdout=subprocess.PIPE)

    tmpfile = [None, None, False]  # Dirname, filename, isFIFO.
    try:
        # If a stutter model is available, do BGPredict.
        next_pipe = p_tssv.stdout
        if in_stuttermodel is not None:
            if config.get("bgcorrect", "profiles") == "-":
                # Set up FIFO if BGCorrect gets two pipes in!
                tmpfile[0] = tempfile.mkdtemp()
                tmpfile[1] = os.path.join(tmpfile[0], "fifo.tmp")
                try:
                    os.mkfifo(tmpfile[1])
                    tmpfile[2] = True  # Is a FIFO.
                except AttributeError:
                    pass  # Tempfile instead of FIFO on Windows.
                config.set("bgcorrect", "profiles", tmpfile[1])
                if in_bgprofiles is None:
                    config.set("bgpredict", "outfile", tmpfile[1])
                    out_bgpredict = tmpfile[1]
                else:
                    config.set("bgmerge", "outfile", tmpfile[1])
                    out_bgmerge = tmpfile[1]

            # Read TSSV output into buffer.
            sample_data = p_tssv.communicate()[0]

            # Send data to BGPredict.
            p_bgpredict = subprocess.Popen(
                get_argv("bgpredict", arg_defs, config),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE if out_bgpredict == "-" else None)
            thread = threading.Thread(
                target=write_to_stream,
                args=(p_bgpredict.stdin, sample_data))
            thread.daemon = True
            thread.start()

            if out_bgpredict == "-":
                next_pipe = p_bgpredict.stdout
            else:
                next_pipe = subprocess.PIPE
                if not (tmpfile[2] and in_bgprofiles is None):
                    # Wait for BGPredict to write its output file.
                    # (Not waiting if it is writing to a FIFO.)
                    p_bgpredict.wait()

            # If BGEstimate profiles are available as well, do BGMerge.
            if in_bgprofiles is not None:
                p_bgmerge = subprocess.Popen(
                    get_argv("bgmerge", arg_defs, config),
                    stdin=p_bgpredict.stdout if out_bgpredict == "-" else None,
                    stdout=subprocess.PIPE if out_bgmerge == "-" else None)
                if out_bgpredict == "-":
                    # We won't read BGPredict's output pipe.
                    p_bgpredict.stdout.close()
                if out_bgmerge == "-":
                    next_pipe = p_bgmerge.stdout
                else:
                    next_pipe = subprocess.PIPE
                    if not tmpfile[2]:
                        # Wait for BGMerge to write its output file.
                        # (Not waiting if it is writing to a FIFO.)
                        p_bgmerge.wait()

        # Do BGCorrect with BGPredict and/or BGEstimate profiles.
        if in_stuttermodel is not None or in_bgprofiles is not None:
            p_bgcorrect = subprocess.Popen(
                get_argv("bgcorrect", arg_defs, config),
                stdin=next_pipe,
                stdout=subprocess.PIPE)

            if next_pipe is subprocess.PIPE:
                # Send data to BGCorrect.
                thread = threading.Thread(
                    target=write_to_stream,
                    args=(p_bgcorrect.stdin, sample_data))
                thread.daemon = True
                thread.start()
            else:
                # We won't read the pipe.
                next_pipe.close()
            next_pipe = p_bgcorrect.stdout

        # Seqconvert.
        p_seqconvert = subprocess.Popen(
            get_argv("seqconvert", arg_defs, config),
            stdin=next_pipe,
            stdout=subprocess.PIPE)
        next_pipe.close()  # We won't read the pipe.

        # Run Samplestats on Seqconvert output.
        p_samplestats = subprocess.Popen(
            get_argv("samplestats", arg_defs, config),
            stdin=subprocess.PIPE)

        # Visualise Seqconvert output with Vis.
        p_vis = subprocess.Popen(get_argv("vis", arg_defs, config),
            stdin=subprocess.PIPE)

        # Send Seqconvert output to Samplestats and Vis.
        tee(p_seqconvert.stdout, p_samplestats.stdin, p_vis.stdin)

        # Wait for the last tools to finish.
        p_tssv.wait()
        if in_stuttermodel is not None:
            p_bgpredict.wait()
        if in_stuttermodel is not None and in_bgprofiles is not None:
            p_bgmerge.wait()
        if in_stuttermodel is not None or in_bgprofiles is not None:
            p_bgcorrect.wait()
        p_seqconvert.wait()
        p_samplestats.wait()
        p_vis.wait()
    finally:
        if tmpfile[1] is not None and os.path.exists(tmpfile[1]):
            os.remove(tmpfile[1])
        if tmpfile[0] is not None and os.path.exists(tmpfile[0]):
            os.rmdir(tmpfile[0])
#run_case_sample_analysis


def add_ini_arg(argument, value, config):
    if value is None or value == False:
        return
    if isinstance(value, list):
        value = " ".join(
            "'" + x.replace("'", "\'") + "'" for x in value)
    else:
        value = str(value)
    v = ini_try_get_option(config, NAME, argument)
    if v is None:
        config.set(NAME, argument, value)
    elif value != v:
        raise ValueError(
            "different values for '%s' given in INI file and on the command "
            "line: '%s' and '%s'" % (argument, v, value))
#add_ini_arg


def get_config(args):
    config = read_ini(args.config)
    for argument in args.__dict__:
        add_ini_arg(argument.replace("_","-"), getattr(args, argument), config)
    return config
#get_config


def run_ini(config):
    analysis_functions = {
        "reference-sample": run_ref_sample_analysis,
        "reference-database": run_ref_database_analysis,
        "case-sample": run_case_sample_analysis
    }
    analysis = ini_require_option(config, NAME, "analysis")
    if analysis not in analysis_functions:
        raise ValueError("Unknown analysis type '%s', specify one of: %s" %
            (analysis, ", ".join(analysis_functions)))
    analysis_functions[analysis](get_arguments(get_tools(True)), config)
#run_ini


def write_ini(args):
    tools = get_tools(False)
    arglists = get_arguments(tools)
    for name in ANALYSIS_TOOLS[args.analysis]:
        if name == NAME:
            fmt = "%%-%is" % max(map(len, ANALYSIS_ARGS[args.analysis] or [""]))
        else:
            fmt = "%%-%is" % reduce(max,
                (len(x[0]) for x in arglists[name] if x[0] not in HIDDEN_ARGS), 0)
        if fmt == "%-0s":
            continue  # No arguments.
        args.config.write("[%s]\n%s\n\n" % (name,
            format_help(tools[name].__doc__.split("\n\n", 1)[0].strip())))
        for arg in arglists[name]:
            if name == NAME:
                if arg[0] not in ANALYSIS_ARGS[args.analysis]:
                    continue
                if getattr(args, arg[0].replace("-","_")) is not None:
                    arg[1]["default"] = getattr(args, arg[0].replace("-","_"))
            elif arg[0] in HIDDEN_ARGS:
                continue  # Argument is being overridden by this tool.
            if "help" in arg[1]:
                args.config.write("%s\n" % format_help(arg[1]["help"], arg[1]))
            if "default" in arg[1]:
                args.config.write(
                    "%s = %s\n\n" % (fmt % arg[0], arg[1]["default"]))
            else:
                args.config.write("%s =\n\n" % (fmt % arg[0]))
        args.config.write("\n")
#write_ini


def add_arguments(parser):
    parser.add_argument("config", metavar="INI",
        help="pipeline configuration file; if it does not exist, a new file "
             "with default settings will be created")
    parser.add_argument("-a", "--analysis", metavar="ANALYSIS",
        choices=("reference-sample", "reference-database", "case-sample"),
        help="controls which predefined analysis pipeline will be run; "
             "'reference-sample' runs a single sample's FastA/FastQ file "
             "through TSSV and Stuttermark to prepare it for the reference-"
             "database analysis; 'reference-database' runs a collection of "
             "reference samples through Allelefinder, BGEstimate, and "
             "Stuttermodel to create a reference database of systemic noise; "
             "'case-sample' runs a single sample's FastA/FastQ file through "
             "TSSV, BGPredict, BGCorrect, and Samplestats")
    group = parser.add_argument_group("sample tag parsing options",
        "these options are used to extract sample tags (names) from their "
        "file names; for details about REGEX syntax and capturing groups, "
        "check https://docs.python.org/howto/regex")
    group.add_argument("-e", "--tag-expr", metavar="REGEX",
        help="regular expression that captures (using one or more "
             "capturing groups) the sample tags from the file names; by "
             "default, the entire file name except for its extension (if "
             "any) is captured")
    group.add_argument("-f", "--tag-format", metavar="EXPR",
        help="format of the sample tags produced; a capturing group "
             "reference like '\\n' refers to the n-th capturing group in "
             "the regular expression specified with -e/--tag-expr (the "
             "default of '\\1' simply uses the first capturing group); "
             "with a single sample, you can enter the sample tag here "
             "explicitly")
    group = parser.add_argument_group("input/output file options",
        "words in [brackets] indicate applicable analysis types; "
        "all of these values can also be specified in the [%s] section of "
        "the INI file" % NAME)
    group.add_argument("-l", "--in-library", metavar="LIBRARY",
        help="library file containing marker definitions")
    group.add_argument("-s", "--in-sample-raw", metavar="FASTA",
        help="[ref-sample, case-sample] FastA or FastQ file containing raw "
             "sequence data of the sample")
    group.add_argument("-m", "--in-stuttermodel", metavar="STUT",
        help="[case-sample] file containing a trained stutter model")
    group.add_argument("-p", "--in-bgprofiles", metavar="PROFILES",
        help="[case-sample] file containing noise profiles from BGEstimate")
    group.add_argument("-r", "--store-predictions", action="store_true",
        help="[case-sample] if this option is specified, output files named "
             "'sampletag-bgpredict.txt' and 'sampletag-bgmerge.txt' will be "
             "created if applicable; these files contain predicted stutter "
             "amounts for the sequences in the sample based on the given "
             "stutter model")
    group.add_argument("-S", "--in-samples", metavar="SAMPLE", nargs="+",
        help="[ref-database] file names of reference sample data files "
             "('.csv' output files of the 'reference-sample' analysis)")
    group.add_argument("-A", "--in-allelelist", metavar="ALLELEFILE",
        help="[ref-database] file containing a list of the true alleles of "
             "each sample; if not given, Allelefinder will be run as part of "
             "the pipeline to create this file; it is ESSENTIAL that you "
             "check the correctness and completeness of the allele list")
    group.add_argument("-P", "--prefix",
        help="[ref-database] if specified, all output file names are prefixed "
             "with this value")
    group.add_argument("-C", "--combine-strands", action="store_true",
        help="[ref-database, case-sample] if specified, noise analysis will be performed on the "
             "total number of reads, instead of separately for either strand")
#add_arguments


def run(args):
    # Export debug flag to the global scope.
    global debug
    debug = args.debug

    if args.config == "-":
        raise ValueError("The pipeline configuration file cannot be named '-'")
    if os.path.exists(args.config):
        # Run analysis pipeline using configuration file.
        args.config = open(args.config, "rt", encoding="UTF-8")
        run_ini(get_config(args))
        args.config.close()
    elif args.analysis is not None:
        # Create default configuration file for pipeline.
        print("Configuration file '%s' does not exist, creating default "
              "configuration file for %s analysis..." %
              (args.config, args.analysis))
        args.config = open(args.config, "wt", encoding="UTF-8")
        write_ini(args)
        args.config.close()
    else:
        raise ValueError(
            "Configuration file '%s' does not exist; if you wish to create a "
            "default configuration file, please specify the type of analysis "
            "pipeline using the -a/--analysis option" % args.config)
#run