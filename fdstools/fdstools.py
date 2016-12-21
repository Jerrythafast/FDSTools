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

import argparse, pkgutil, os, re, textwrap
#import cProfile  # Imported only if the -d/--debug option is specified
import tools

from . import usage, version


# Map tool names to argument parsers.
__tools__ = {}


class _HelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If valid 'tool' argument is given, call that tool's help.
        if values and values[0] in __tools__:
            __tools__[values[0]].print_help()
            __tools__[values[0]].exit()
        else:
            parser.print_help()
            parser.exit()
    #__call__
#_HelpAction


class _VersionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If valid 'tool' argument is given, print that tool's version.
        if values and values[0] in __tools__:
            parser = __tools__[values[0]]
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
#_HelpFormatter


def main():
    """
    Main entry point.
    """
    prog = os.path.splitext(os.path.basename(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=_HelpFormatter, prog=prog,
                                     add_help=False, description=usage[0])
    parser.version = version(prog)
    parser.add_argument('-h', '--help', action=_HelpAction,
                        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
                        help="show this help message, or help for the "
                             "specified TOOL, and exit")
    parser.add_argument('-v', "--version", action=_VersionAction,
                        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
                        help="show version number and exit")
    parser.add_argument('-d', "--debug", action="store_true",
                        help="if specified, additional debug output is given")
    subparsers = parser.add_subparsers(title='available tools', dest='tool',
                                       metavar='TOOL', help="specify which "
                                       "tool to run")
    subparsers.required = True

    prefix = tools.__name__ + "."
    for importer, name, ispkg in pkgutil.iter_modules(
            path=[os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "tools")]):
        module = importer.find_module(prefix + name).load_module(prefix + name)
        subparser = subparsers.add_parser(
            name,
            formatter_class=_HelpFormatter,
            help=module.__doc__.split("\n\n", 1)[0],
            description=module.__doc__,
            version=version(prog, name, module.__version__))
        __tools__[name] = subparser
        subparser.add_argument('-d', "--debug", action="store_true",
            default=argparse.SUPPRESS,
            help="if specified, additional debug output is given")
        module.add_arguments(subparser)
        subparser.set_defaults(func=module.run)
    try:
        args, unknowns = parser.parse_known_args()
    except Exception as error:
        parser.error(error)
    try:
        if unknowns:
            # Politely inform the user about unknown arguments.
            __tools__[args.tool].error(
                "The following arguments are not known. Please check spelling "
                "and argument order: '%s'." % "', '".join(unknowns))
        if args.debug:
            import cProfile
            cProfile.runctx(
                "args.func(args)", globals(), locals(), sort="tottime")
        else:
            args.func(args)
    except Exception as error:
        if args.debug:
            raise
        __tools__[args.tool].error(error)
#main


if __name__ == "__main__":
    main()
