#!/usr/bin/env python
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

Vega also supports generating visualisations on the command line, which
is useful if you wish to generate graphs automatically in your analysis
pipeline.  Specify the -V/--vega option to obtain a bare Vega graph
specification (a JSON file) instead of the full-featured HTML file.  You
can pass this file through Vega to generate a PNG or SVG image file.

If an input file is specified, the visualisation will be set up
specifically to visualise the contents of this file.  To this end, the
entire file contents are embedded in the generated visualisation.
"""
import argparse
import sys
import json
import re
import os
import cgi

from pkg_resources import resource_stream, resource_string

from ..lib import pos_int_arg

__version__ = "0.1dev"


# Default values for parameters are specified below.

# Default minimum number of reads to require.
# This value can be overridden by the -n command line option.
_DEF_THRESHOLD_ABS = 15

# Default minimum amount of reads to require, as a percentage of the
# highest allele of each marker.
# This value can be overridden by the -m command line option.
_DEF_THRESHOLD_PCT = 0.5

# Default width of bars in bar graphs, in pixels.
# This value can be overridden by the -b command line option.
_DEF_BAR_WIDTH = 15

# Default amount of padding between subgraphs, in pixels.
# This value can be overridden by the -p command line option.
_DEF_SUBGRAPH_PADDING = 70

# Default graph width in pixels.
# This value can be overridden by the -w command line option.
_DEF_WIDTH = 600

# Default marker name matching regular expression.
# This value can be overridden by the -M command line option.
_DEF_MARKER_REGEX = ".*"

# Default data file that Vega will read when -V/--vega is specified
# without providing data to embed in the file.
# It is currently impossible to override this value.
_DEF_DATA_FILENAME = "data.csv"



_PAT_LOAD_SCRIPT = re.compile("<!--\s*BEGIN_LOAD_SCRIPT\s*-->\s*(.*?)\s*"
                              "<!--\s*END_LOAD_SCRIPT\s*-->", flags=re.DOTALL)
_PAT_LIBRARIES = re.compile("<!--\s*BEGIN_LIBRARIES\s*-->\s*(.*?)\s*"
                            "<!--\s*END_LIBRARIES\s*-->", flags=re.DOTALL)

_SCRIPT_BEGIN = '<script type="text/javascript">'
_SCRIPT_END = '</script>'

_EXTERNAL_LIBRARIES = ("vis/d3.min.js", "vis/vega.min.js")



def set_data_formula_transform_value(spec, dataname, fieldname, value):
    for data in spec["data"]:
        if data["name"] != dataname:
            continue
        if "transform" not in data:
            return False
        for transform in data["transform"]:
            if (transform["type"] == "formula" and
                    transform["field"] == fieldname):
                transform["expr"] = str(value)
                return True
    return False
#set_data_formula_transform_value


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


def create_visualisation(vistype, infile, outfile, vega, online, tidy, min_abs,
                         min_pct, bar_width, padding, marker, width,
                         log_scale):
    # Get graph spec.
    spec = json.load(resource_stream(
        "fdstools", "vis/%svis/%svis.json" % (vistype, vistype)))
    if infile is not None:
        # Embed the given data file.
        spec["data"][0]["values"] = infile.read()
    elif vega:
        # Vega should load data from somewhere in headless mode.
        del spec["data"][0]["values"]
        spec["data"][0]["url"] = _DEF_DATA_FILENAME

    # Apply settings.
    spec["width"] = width
    set_data_formula_transform_value(spec, "yscale", "barwidth", bar_width)
    set_data_formula_transform_value(spec, "yscale", "subgraphoffset", padding)
    set_data_formula_transform_value(
        spec, "table", "filter_marker", "'" + marker + "'")
    if vistype == "sample" or vistype == "bgraw":
        set_data_formula_transform_value(
            spec, "table", "amplitude_threshold", min_abs)
        set_data_formula_transform_value(
            spec, "table", "amplitude_pct_threshold", min_pct)
    elif vistype == "profile":
        set_data_formula_transform_value(
            spec, "table", "filter_threshold", min_pct)
        set_data_formula_transform_value(
            spec, "table", "low", "0.001" if log_scale else "0")
    if not log_scale:
        set_axis_scale(spec, "x", "linear")
    elif vistype == "sample":
        set_axis_scale(spec, "x", "sqrt")
    else:
        set_axis_scale(spec, "x", "log")

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
    html = resource_string("fdstools", "vis/%svis/index.html" % vistype)
    match = _PAT_LOAD_SCRIPT.search(html)
    if match:
        html = "".join([html[:match.start(1)],
                        _SCRIPT_BEGIN,
                        "var graph_spec=",
                        spec,
                        ";onLoadSpec(",
                        "true" if infile is not None else "false",
                        ");",
                        _SCRIPT_END,
                        html[match.end(1):]])

    if not online:
        # Replace external libraries with inline libraries.
        match = _PAT_LIBRARIES.search(html)
        if match:
            parts = [html[:match.start(1)]]
            for library in _EXTERNAL_LIBRARIES:
                parts += [_SCRIPT_BEGIN,
                          resource_string("fdstools", library),
                          _SCRIPT_END]
            parts.append(html[match.end(1):])
            html = "".join(parts)

    outfile.write(html)    
#create_visualisation


def add_arguments(parser):
    parser.add_argument('type', metavar="TYPE",
        choices=("sample", "profile", "bgraw"),
        help="the type of data to visualise; use 'sample' to visualise "
             "sample data files and BGCorrect output; use 'profile' to "
             "visualise background noise profiles obtained with BGEstimate, "
             "BGHomStats, and BGPredict; use 'bgraw' to visualise raw "
             "background noise data obtained with BGHomRaw")
    parser.add_argument('infile', metavar="IN", nargs="?",
        help="file containing the data to embed in the visualisation file; if "
             "not specified, HTML visualisation files will contain a file "
             "selection control, and Vega visualisation files will load data "
             "from a file called '%s'" % _DEF_DATA_FILENAME)
    parser.add_argument('outfile', metavar="OUT", nargs="?",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="file to write output to (default: write to stdout)")
    parser.add_argument('-V', '--vega', action="store_true",
        help="an HTML file containing an interactive visualisation is created "
             "by default; if this option is specified, a Vega graph "
             "specification (JSON file) is produced instead")
    parser.add_argument('-O', '--online', action="store_true",
        help="when generating an HTML visualisation file, required JavaScript "
             "libraries (D3 and Vega) are embedded in the file; if this "
             "option is specified, the HTML file will instead link to these "
             "libraries on the Internet, thereby always using the latest "
             "versions of D3 and Vega")
    parser.add_argument('-t', '--tidy', action="store_true",
        help="tidily indent the generated JSON")

    visgroup = parser.add_argument_group("visualisation options",
        description="words in [brackets] indicate applicable visualisation "
                    "types")
    visgroup.add_argument('-n', '--min-abs', metavar="N", type=pos_int_arg,
        default=_DEF_THRESHOLD_ABS,
        help="[sample, bgraw] only show sequences with this minimum number of "
             "reads (default: %(default)s)")
    visgroup.add_argument('-m', '--min-pct', metavar="PCT", type=float,
        default=_DEF_THRESHOLD_PCT,
        help="[sample, profile, bgraw] for sample: only show sequences with "
             "at least this percentage of the number of reads of the highest "
             "allele of a marker; for profile and bgraw: at least this "
             "percentage of the true allele (default: %(default)s)")
    visgroup.add_argument('-M', '--marker', metavar="REGEX",
        default=_DEF_MARKER_REGEX,
        help="[sample, profile, bgraw] only show graphs for the markers that "
             "match the given regular expression; the default value "
             "'%(default)s' matches any marker name")
    visgroup.add_argument('-L', '--log-scale', action="store_true",
        help="[sample, profile, bgraw] use logarithmic scale (for sample: "
             "square root scale) instead of linear scale")
    visgroup.add_argument('-b', '--bar-width', metavar="N", type=pos_int_arg,
        default=_DEF_BAR_WIDTH,
        help="[sample, profile, bgraw] width of the bars in pixels (default: "
             "%(default)s)")
    visgroup.add_argument('-p', '--padding', metavar="N", type=pos_int_arg,
        default=_DEF_SUBGRAPH_PADDING,
        help="[sample, profile, bgraw] amount of padding (in pixels) between "
             "graphs of different markers/alleles (default: %(default)s)")
    visgroup.add_argument('-w', '--width', metavar="N", type=pos_int_arg,
        default=_DEF_WIDTH,
        help="[sample, profile, bgraw] width of the graph area in pixels "
             "(default: %(default)s)")
#add_arguments


def run(args):
    if args.infile == "-" and not sys.stdin.isatty():
        # User appears to want to pipe data in.
        args.infile = sys.stdin
    if (args.infile is not None and args.outfile == sys.stdout
            and not os.path.exists(args.infile)):
        # One filename given, and it does not exist.  Assume outfile.
        args.outfile = open(args.infile, 'w')
        args.infile = None

    if args.outfile.isatty():
        raise ValueError("Please specify a file name to write the %s to." %
                         ("Vega graph specification (JSON format)" if args.vega
                          else "HTML document"))

    if args.infile is not None and args.infile != sys.stdin:
        # Open the specified input file.
        args.infile = open(args.infile, 'r')

    create_visualisation(args.type, args.infile, args.outfile, args.vega,
                         args.online, args.tidy, args.min_abs, args.min_pct,
                         args.bar_width, args.padding, args.marker, args.width,
                         args.log_scale)
#run
