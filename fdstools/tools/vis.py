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
Create a data visualisation web page or Vega graph specification.

With no optional arguments specified, a self-contained web page (HTML
file) is produced.  You can open this file in a web browser to view
interactive visualisations of your data.  The web page contains a file
selection element which can be used to select the data to be visualised.

Visualisations make use of the Vega JavaScript library
(https://vega.github.io).  The required JavaScript libraries (Vega and
D3) are embedded in the generated HTML file.  With the -O/--online
option specified, the HTML file will instead link to the latest version
of these libraries on the Internet.

Vega supports generating visualisations on the command line.  By
default, FDSTools produces a full-featured HTML file.  Specify the
-V/--vega option if you wish to obtain a bare Vega graph specification
(a JSON file) instead.  You can pass this file through Vega to generate
a PNG or SVG image file.

If an input file is specified, the visualisation will be set up
specifically to visualise the contents of this file.  To this end, the
entire file contents are embedded in the generated visualisation.
"""
import argparse
import json
import os.path
import pkgutil
import re
import sys

from errno import EPIPE

from ..lib.cli import pos_int_arg

__version__ = "1.1.0"


# Default values for parameters are specified below.

# Default minimum number of reads to require.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS = 5

# Default minimum amount of reads to require, as a percentage of the
# highest allele of each marker.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT_OF_MAX = 0.5

# Default minimum amount of reads to require, as a percentage of the
# total number of reads each marker.
# This value can be overridden by the -S command line option.
_DEF_THRESHOLD_PCT_OF_SUM = 0.0

# Default minimum number of reads per orientation to require.
# This value can be overridden by the -s command line option.
_DEF_THRESHOLD_ORIENTATION = 0

# Default minimum number of reads to mark as allele in Samplevis.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS_ALLELE = 30

# Default minimum percentage of reads w.r.t. the highest allele of the
# marker to mark as allele in Samplevis.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT_OF_MAX_ALLELE = 2.

# Default minimum percentage of reads w.r.t. the marker's total number
# of reads to mark as allele in Samplevis.
# This value can be overridden by the -p command line option.
_DEF_THRESHOLD_PCT_OF_SUM_ALLELE = 1.5

# Default minimum percentage of correction to mark as allele in Samplevis.
# This value can be overridden by the -c command line option.
_DEF_THRESHOLD_CORRECTION_ALLELE = 0

# Default minimum number of recovered reads as a percentage of the
# number of reads after correction to mark as allele in Samplevis.
# This value can be overridden by the -r command line option.
_DEF_THRESHOLD_RECOVERY_ALLELE = 0

# Default minimum number of reads per strand to mark as allele in Samplevis.
# This value can be overridden by the -Z command line option.
_DEF_THRESHOLD_ORIENTATION_ALLELE = 0

# Default percentage of reads on one strand to mark as bias.
# This value can be overridden by the -B command line option.
_DEF_THRESHOLD_BIAS = 0.0

# Default width of bars in bar graphs, in pixels.
# This value can be overridden by the -b command line option.
_DEF_BAR_WIDTH = 15

# Default amount of padding between subgraphs, in pixels.
# This value can be overridden by the -p command line option.
_DEF_SUBGRAPH_PADDING = 70

# Default maximum length of sequences to display.
# This value can be overridden by the -x command line option.
_DEF_MAX_SEQ_LEN = 70

# Default graph width in pixels.
# This value can be overridden by the -w command line option.
_DEF_WIDTH = 600

# Default graph height in pixels.
# This value can be overridden by the -H command line option.
_DEF_HEIGHT = 400

# Default amount of jitter.
# This value can be overridden by the -j command line option.
_DEF_JITTER = 0.25

# Default marker name matching regular expression.
# This value can be overridden by the -M command line option.
_DEF_MARKER_NAME = ""

# Default repeat unit matching regular expression.
# This value can be overridden by the -U command line option.
_DEF_UNIT = ""

# Default data file that Vega will read when -V/--vega is specified
# without providing data to embed in the file.
# It is currently impossible to override this value.
_DEF_DATA_FILENAME = "data.csv"


_PAT_LIBRARIES = re.compile("<!--\s*BEGIN_LIBRARIES\s*-->\s*(.*?)\s*"
                            "<!--\s*END_LIBRARIES\s*-->", flags=re.DOTALL)
_PAT_LOAD_SCRIPT = re.compile("<!--\s*BEGIN_LOAD_SCRIPT\s*-->\s*(.*?)\s*"
                              "<!--\s*END_LOAD_SCRIPT\s*-->", flags=re.DOTALL)

_SCRIPT_BEGIN = '<script type="text/javascript">'
_SCRIPT_END = "</script>"

_EXTERNAL_LIBRARIES = ("vis/d3.min.js", "vis/vega.min.js")


def set_signal_value(spec, signalname, value):
    if "signals" not in spec:
        return False
    for signal in spec["signals"]:
        if signal["name"] == signalname:
            signal["init"] = value
            return True
    return False
#set_signal_value


def set_axis_scale(spec, scalename, value):
    success = False
    for marks in spec["marks"]:
        if "scales" not in marks:
            continue
        for scale in marks["scales"]:
            if scale["name"] != scalename:
                continue
            scale["type"] = value;
            success = True
    return success
#set_axis_scale


def create_visualisation(vistype, infile, infile2, outfile, vega, tidy, min_abs, min_pct_of_max,
                         min_pct_of_sum, min_per_strand, bias_threshold, bar_width, padding,
                         marker, width, height, log_scale, repeat_unit, no_alldata, no_aggregate,
                         no_ce_length_sort, max_seq_len, jitter, title, allele_min_abs,
                         allele_min_pct_of_max, allele_min_pct_of_sum, allele_min_correction,
                         allele_min_recovery, allele_min_per_strand):
    # Get graph spec.
    spec = json.loads(
        pkgutil.get_data("fdstools", "vis/%svis/%svis.json" % (vistype, vistype)).decode())
    if infile is not None:
        # Embed the given data file.
        spec["data"][0]["values"] = infile.read()
    elif vega:
        # Vega should load data from somewhere in headless mode.
        del spec["data"][0]["values"]
        spec["data"][0]["url"] = _DEF_DATA_FILENAME
    if vistype in ("profile", "stuttermodel") and infile2 is not None:
        # Embed given raw data points file.
        spec["data"][1]["values"] = infile2.read()

    # Apply width, height, and padding settings.
    spec["width"] = width
    if vistype == "stuttermodel":
        set_signal_value(spec, "graphheight", height)
        set_signal_value(spec, "jitter", jitter)
    elif vistype == "allele":
        spec["height"] = height
    else:
        set_signal_value(spec, "barwidth", bar_width)
    if vistype != "allele":
        if vistype != "bganalyse":
            set_signal_value(spec, "subgraphoffset", padding)
        set_signal_value(spec, "filter_marker", marker)

    # Apply type-specific settings.
    if vistype == "stuttermodel":
        set_signal_value(spec, "filter_unit", repeat_unit)
        set_signal_value(spec, "show_all_data", False if no_alldata else True)
    elif vistype == "sample" or vistype == "bgraw" or vistype == "profile":
        set_signal_value(spec, "amplitude_threshold", min_abs)
        set_signal_value(spec, "amplitude_pct_threshold", min_pct_of_max)
    if vistype == "profile":
        set_signal_value(spec, "low", 0.001 if log_scale else 0)
    elif vistype == "sample":
        set_signal_value(spec, "orientation_threshold", min_per_strand)
        set_signal_value(spec, "bias_threshold", bias_threshold)
        set_signal_value(spec, "amplitude_markerpct_threshold", min_pct_of_sum)
        set_signal_value(spec, "show_other", not no_aggregate)
        set_signal_value(spec, "sort_str_by_length", not no_ce_length_sort)
        set_signal_value(spec, "max_seq_len", max_seq_len)
        set_signal_value(spec, "allele_amplitude_threshold", allele_min_abs)
        set_signal_value(spec, "allele_amplitude_pct_threshold", allele_min_pct_of_max)
        set_signal_value(spec, "allele_amplitude_markerpct_threshold", allele_min_pct_of_sum)
        set_signal_value(spec, "allele_correction_threshold", allele_min_correction)
        set_signal_value(spec, "allele_recovery_threshold", allele_min_recovery)
        set_signal_value(spec, "allele_orientation_threshold", allele_min_per_strand)

    # Apply axis scale settings.
    if vistype != "stuttermodel" and vistype != "allele":
        if not log_scale:
            set_axis_scale(spec, "x", "linear")
        elif vistype == "sample" or vistype == "bganalyse":
            set_axis_scale(spec, "x", "sqrt")
        else:
            set_axis_scale(spec, "x", "log")

    # Add title if available.
    if title is None and infile is not None and infile != sys.stdin:
        try:
            title = os.path.splitext(os.path.basename(infile.name))[0]
        except AttributeError:
            pass
    if title:
        spec["data"][0]["fdstools_filename"] = title;

    # Stringify spec.
    if tidy:
        spec = json.dumps(spec, indent=2, separators=(",", ": "))
    else:
        spec = json.dumps(spec, separators=(",", ":"))

    if vega:
        # Vega graph spec is all that we were asked for.
        outfile.write(spec)
        return

    # Creating a fully self-contained HTML visualisation instead.
    html = pkgutil.get_data("fdstools", "vis/%svis/index.html" % vistype).decode()
    match = _PAT_LOAD_SCRIPT.search(html)
    if match:
        html = "".join([
            html[:match.start(1)],
            _SCRIPT_BEGIN,
            "var graph_spec=",
            spec,
            ";onLoadSpec(",
            "true" if infile is not None else "false",
            "" if vistype not in ("profile", "stuttermodel") else
                ", true" if infile2 is not None else ", false",
            ");",
            _SCRIPT_END,
            html[match.end(1):]])

    # Replace external libraries with inline libraries.
    match = _PAT_LIBRARIES.search(html)
    if match:
        parts = [html[:match.start(1)]]
        for library in _EXTERNAL_LIBRARIES:
            parts += [_SCRIPT_BEGIN, pkgutil.get_data("fdstools", library).decode(), _SCRIPT_END]
        parts.append(html[match.end(1):])
        html = "".join(parts)

    outfile.write(html)
#create_visualisation


def add_arguments(parser):
    parser.add_argument("type", metavar="TYPE",
        choices=("sample", "profile", "bgraw", "stuttermodel", "allele", "bganalyse"),
        help="the type of data to visualise; use 'sample' to visualise "
             "sample data files and bgcorrect output; use 'profile' to "
             "visualise background noise profiles obtained with bgestimate, "
             "bghomstats, and bgpredict; use 'bgraw' to visualise raw "
             "background noise data obtained with bghomraw; use "
             "'stuttermodel' to visualise models of stutter obtained from "
             "stuttermodel; 'bganalyse' to visualise data obtained from "
             "bganalyse; use 'allele' to visualise the allele list obtained "
             "from allelefinder")
    parser.add_argument("infile", metavar="IN", nargs="?", default=sys.stdin,
        help="file containing the data to embed in the visualisation file; if "
             "not specified, HTML visualisation files will contain a file "
             "selection control, and Vega visualisation files will load data "
             "from a file called '%s'" % _DEF_DATA_FILENAME)
    parser.add_argument("outfile", metavar="OUT", nargs="?",
        type=argparse.FileType("tw", encoding="UTF-8"), default=sys.stdout,
        help="file to write output to (default: write to stdout)")
    parser.add_argument("-V", "--vega", action="store_true",
        help="by default, a full-featured HTML file offering an interactive "
             "visualisation is created; if this option is specified, only a "
             "bare Vega graph specification (JSON file) is produced instead")
    parser.add_argument("-t", "--tidy", action="store_true",
        help="tidily indent the generated JSON")
    parser.add_argument("-T", "--title",
        help="prepend the given value to the title of HTML visualisations "
             "(default: prepend name of data file if given)")

    visgroup = parser.add_argument_group("visualisation options",
        description="words in [brackets] indicate applicable visualisation types")
    visgroup.add_argument("-n", "--min-abs", metavar="N", type=float, default=_DEF_THRESHOLD_ABS,
        help="[sample, profile, bgraw] only show sequences with this minimum "
             "number of reads (default: %(default)s)")
    visgroup.add_argument("-m", "--min-pct-of-max", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT_OF_MAX,
        help="[sample, profile, bgraw] for sample: only show sequences with "
             "at least this percentage of the number of reads of the highest "
             "allele of a marker; for profile and bgraw: at least this "
             "percentage of the true allele (default: %(default)s)")
    visgroup.add_argument("-S", "--min-pct-of-sum", metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT_OF_SUM,
        help="[sample] only show sequences with at least this percentage of "
             "the total number of reads of a marker (default: %(default)s)")
    visgroup.add_argument("-s", "--min-per-strand", metavar="N", type=float,
        default=_DEF_THRESHOLD_ORIENTATION,
        help="[sample] only show sequences with this minimum number of reads "
             "for both orientations (forward/reverse) (default: %(default)s)")
    visgroup.add_argument("-B", "--bias-threshold", metavar="N", type=float,
        default=_DEF_THRESHOLD_BIAS,
        help="[sample] mark sequences that have less than this percentage of "
             "reads on one strand (default: %(default)s)")
    visgroup.add_argument("-c", "--no-ce-length-sort", action="store_true",
        help="[sample] if specified, do not sort STR alleles by length")
    visgroup.add_argument("-M", "--marker", metavar="MARKER", default=_DEF_MARKER_NAME,
        help="[sample, profile, bgraw, stuttermodel, bganalyse] only show "
             "graphs for the markers that contain the given value in their "
             "name; separate multiple values with spaces; prepend any value "
             "with '=' for an exact match (default: show all markers)")
    visgroup.add_argument("-U", "--repeat-unit", metavar="UNIT", default=_DEF_UNIT,
        help="[stuttermodel] only show graphs for the repeat units that "
             "contain the given value; separate multiple values with spaces; "
             "prepend any value with '=' for an exact match (default: show "
             "all repeat units)")
    visgroup.add_argument("-A", "--no-alldata", action="store_true",
        help="[stuttermodel] if specified, show only marker-specific fits")
    visgroup.add_argument("-a", "--no-aggregate", action="store_true",
        help="[sample] if specified, do not replace filtered sequences with a "
             "per-marker aggregate 'Other sequences' entry")
    visgroup.add_argument("-I", "--input2", dest="infile2", metavar="FILE",
        type=argparse.FileType("tr", encoding="UTF-8"),
        help="[profile, stuttermodel] raw data points file to overlay on the "
             "background noise profiles or stutter model graphs (as obtained "
             "from bghomraw or the -r/--raw-outfile option of stuttermodel); "
             "if not specified, HTML visualisation files will contain a file "
             "selection control")

    dispgroup = parser.add_argument_group("display options")
    dispgroup.add_argument("-L", "--log-scale", action="store_true",
        help="[sample, profile, bgraw, bganalyse] use logarithmic scale (for "
             "sample and bganalyse: square root scale) instead of linear scale")
    dispgroup.add_argument("-b", "--bar-width", metavar="N", type=pos_int_arg,
        default=_DEF_BAR_WIDTH,
        help="[sample, profile, bgraw, bganalyse] width of the bars in pixels "
             "(default: %(default)s)")
    dispgroup.add_argument("-p", "--padding", metavar="N", type=pos_int_arg,
        default=_DEF_SUBGRAPH_PADDING,
        help="[sample, profile, bgraw, stuttermodel] amount of padding (in "
             "pixels) between graphs of different markers/alleles (default: %(default)s)")
    dispgroup.add_argument("-w", "--width", metavar="N", type=pos_int_arg, default=_DEF_WIDTH,
        help="[sample, profile, bgraw, stuttermodel, bganalyse, allele] width "
             "of the graph area in pixels (default: %(default)s)")
    dispgroup.add_argument("-H", "--height", metavar="N", type=pos_int_arg, default=_DEF_HEIGHT,
        help="[stuttermodel, allele] height of the graph area in pixels (default: %(default)s)")
    dispgroup.add_argument("-x", "--max-seq-len", metavar="N", type=pos_int_arg,
        default=_DEF_MAX_SEQ_LEN,
        help="[sample] truncate long sequences to this number of characters "
             "(default: %(default)s)")
    dispgroup.add_argument("-j", "--jitter", metavar="N", type=float, default=_DEF_JITTER,
        help="[stuttermodel] apply this amount of jitter to raw data points "
             "(between 0 and 1, default: %(default)s)")

    allelegroup = parser.add_argument_group("allele calling options",
        "for sample visualisations only; sequences that match the -C or -Y "
        "option (or both) and all of the other settings are marked as 'allele'")
    allelegroup.add_argument("-N", "--allele-min-abs", metavar="N", type=float,
        default=_DEF_THRESHOLD_ABS_ALLELE,
        help="the minimum number of reads (default: %(default)s)")
    allelegroup.add_argument("-X", "--allele-min-pct-of-max", metavar="PCT",
        type=float, default=_DEF_THRESHOLD_PCT_OF_MAX_ALLELE,
        help="the minimum percentage of reads w.r.t. the highest allele of "
             "the marker (default: %(default)s)")
    allelegroup.add_argument("-Q", "--allele-min-pct-of-sum", metavar="PCT",
        type=float, default=_DEF_THRESHOLD_PCT_OF_SUM_ALLELE,
        help="the minimum percentage of reads w.r.t. the marker's total "
             "number of reads (default: %(default)s)")
    allelegroup.add_argument("-C", "--allele-min-correction", metavar="N",
        type=float, default=_DEF_THRESHOLD_CORRECTION_ALLELE,
        help="the minimum change in read count due to correction by e.g., "
             "bgcorrect (default: %(default)s)")
    allelegroup.add_argument("-Y", "--allele-min-recovery", metavar="N",
        type=float, default=_DEF_THRESHOLD_RECOVERY_ALLELE,
        help="the minimum number of reads that was recovered thanks to "
             "noise correction (by e.g., bgcorrect), as a percentage of the "
             "total number of reads after correction (default: %(default)s)")
    allelegroup.add_argument("-Z", "--allele-min-per-strand", metavar="N",
        type=float, default=_DEF_THRESHOLD_ORIENTATION_ALLELE,
        help="the minimum number of reads in both orientations (default: %(default)s)")
#add_arguments


def run(args):
    if args.infile == "-":
        args.infile = sys.stdin
    if args.infile != sys.stdin and args.outfile == sys.stdout and not os.path.exists(args.infile):
        # One filename given, and it does not exist.  Assume outfile.
        args.outfile = open(args.infile, "wt", encoding="UTF-8")
        args.infile = sys.stdin
    elif args.outfile.isatty():
        raise ValueError("Please specify a file name to write the %s to." %
            ("Vega graph specification (JSON format)" if args.vega else "HTML document"))

    if args.infile != sys.stdin:
        # Open the specified input file.
        args.infile = open(args.infile, "rt", encoding="UTF-8")
    elif args.infile.isatty():
        # No input given.  Produce an empty visualisation.
        args.infile = None

    try:
        create_visualisation(args.type, args.infile, args.infile2, args.outfile, args.vega,
                             args.tidy, args.min_abs, args.min_pct_of_max, args.min_pct_of_sum,
                             args.min_per_strand, args.bias_threshold, args.bar_width,
                             args.padding, args.marker, args.width, args.height, args.log_scale,
                             args.repeat_unit, args.no_alldata, args.no_aggregate,
                             args.no_ce_length_sort, args.max_seq_len, args.jitter, args.title,
                             args.allele_min_abs, args.allele_min_pct_of_max,
                             args.allele_min_pct_of_sum, args.allele_min_correction,
                             args.allele_min_recovery, args.allele_min_per_strand)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#run
