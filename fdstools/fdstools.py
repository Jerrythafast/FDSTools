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

import argparse
import importlib
import logging
import pkgutil
import re
import sys
import textwrap

from pathlib import Path

from . import tools, usage, version


logging.basicConfig(
    format="[%(name)s][%(levelname)s] %(message)s",
    level=logging.DEBUG if "-d" in sys.argv or "--debug" in sys.argv else logging.INFO)
logging.captureWarnings(True)
LOG = logging.getLogger(__name__)


# Map tool names to argument parsers.
TOOL_SUBPARSERS = {}


class _HelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If valid 'tool' argument is given, call that tool's help.
        if values and values[0] in TOOL_SUBPARSERS:
            TOOL_SUBPARSERS[values[0]].print_help()
        else:
            parser.print_help()
        sys.exit()
    #__call__
#_HelpAction


class _VersionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If valid 'tool' argument is given, print that tool's version.
        if values and values[0] in TOOL_SUBPARSERS:
            parser = TOOL_SUBPARSERS[values[0]]
        formatter = parser._get_formatter()
        formatter.add_text(parser.version)
        parser.exit(message=formatter.format_help())
    #__call__
#_VersionAction


class _HelpFormatter(argparse.HelpFormatter):
    _pat_paragraph_delim = re.compile("\n\n+")
    def _fill_text(self, text, width, indent):
        # Reflow (wrap) description text, but maintain paragraphs.
        return "\n\n".join(
            textwrap.fill(self._whitespace_matcher.sub(" ", p).strip(), width,
                          initial_indent=indent, subsequent_indent=indent)
            for p in self._pat_paragraph_delim.split(text))
    #_fill_text

    def _format_args(self, action, default_metavar):
        if action.dest in ("help", "version"):
            return ""  # Don't display '...' for -h and -v.
        return super()._format_args(action, default_metavar)
    #_format_args

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        # Show metavar once for each option.
        return ", ".join(action.option_strings) + " " + self._format_args(
            action, self._get_default_metavar_for_optional(action))
#_HelpFormatter


def main():
    """
    Main entry point.
    """
    prog = Path(__file__).stem
    parser = argparse.ArgumentParser(formatter_class=_HelpFormatter, prog=prog, add_help=False,
        description=usage[0], exit_on_error=LOG.getEffectiveLevel()>=logging.INFO)
    parser.version = version(prog)
    parser.add_argument("-h", "--help", action=_HelpAction,
        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
        help="show this help message, or help for the specified TOOL, and exit")
    parser.add_argument("-v", "--version", action=_VersionAction,
        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
        help="show version number and exit")
    parser.add_argument("-d", "--debug", action="store_true",
        help="if specified, additional debug output is given")
    subparsers = parser.add_subparsers(title="available tools", dest="tool", metavar="TOOL",
        required=True, help="specify which tool to run")

    # Import each tool and add its arguments to a subparser.
    prefix = tools.__name__ + "."
    for importer, name, ispkg in pkgutil.iter_modules(tools.__path__):
        try:
            module = importlib.import_module(prefix + name)
        except Exception as error:
            LOG.error("FDSTools failed to load '%s': %s", name, error)
            continue
        subparser = subparsers.add_parser(
            name,
            formatter_class=_HelpFormatter,
            help=module.__doc__.split("\n\n", 1)[0],
            description=module.__doc__)
        TOOL_SUBPARSERS[name] = subparser
        subparser.add_argument("-v", "--version", action=_VersionAction, default=argparse.SUPPRESS,
            nargs=0, help="show version number and exit")
        subparser.add_argument("-d", "--debug", action="store_true", default=argparse.SUPPRESS,
            help="if specified, additional debug output is given")
        try:
            subparser.version = version(prog, name, module.__version__)
            module.add_arguments(subparser)
            subparser.set_defaults(func=module.run)
        except Exception as error:
            LOG.error("FDSTools failed to configure '%s': %s", name, error)
            error_to_reraise = error
            def reraise_error(args):
                raise error_to_reraise
            subparser.set_defaults(func=reraise_error)

    # Assume the user wants help if they just type 'fdstools'.
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    LOG.debug("FDSTools version %s", parser.version)
    LOG.debug("sys.argv=%r", sys.argv)

    # Parse the arguments.
    args, unknowns = parser.parse_known_args()
    if unknowns:
        # Politely inform the user about unknown arguments, then exit.
        TOOL_SUBPARSERS[args.tool].error(
            "The following arguments are not known. Please check spelling "
            "and argument order: '%s'." % "', '".join(unknowns))

    if args.debug and LOG.getEffectiveLevel() > logging.DEBUG:
        # User included "d" with more options in a single argument.
        logging.getLogger().setLevel(logging.DEBUG)
        LOG.debug("FDSTools version %s", parser.version)
        LOG.debug("sys.argv=%r", sys.argv)
        LOG.debug("LOG.level was lowered to DEBUG after parsing arguments")

    if LOG.isEnabledFor(logging.DEBUG):
        import warnings
        warnings.simplefilter("default")

    # Run the tool.
    try:
        args.func(args)
    except Exception as error:
        if args.debug:
            raise
        if isinstance(error, BrokenPipeError):
            # This is usually caused by piping through 'head'.
            try:
                LOG.warning("Stopped early because the output was closed.")
            except:
                pass
        else:
            TOOL_SUBPARSERS[args.tool].error(error)
#main


if __name__ == "__main__":
    main()
