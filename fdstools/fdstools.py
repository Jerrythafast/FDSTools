#!/usr/bin/env python

import argparse, pkgutil, os
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


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(add_help=False, description=usage[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.version = version(parser.prog)
    parser.add_argument('-d', "--debug", action="store_true",
                        help="if specified, debug output is printed to stdout")
    parser.add_argument('-v', "--version", action=_VersionAction,
                        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
                        help="show version number and exit")
    parser.add_argument('-h', '--help', action=_HelpAction,
                        default=argparse.SUPPRESS, nargs=argparse.REMAINDER,
                        help="show this help message, or help for the "
                             "specified TOOL, and exit")
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
            name, help=module.__doc__, description=module.__doc__,
            version=version(parser.prog, name, module.__version__))
        __tools__[name] = subparser
        module.add_arguments(subparser)
        subparser.set_defaults(func=module.run)
        subparser.add_argument('-d', "--debug", action="store_true",
            help="if specified, debug output is printed to stdout")
    try:
        args = parser.parse_args()
    except Exception as error:
        parser.error(error)
    try:
        args.func(args)
    except Exception as error:
        if args.debug:
            raise
        __tools__[args.tool].error(error)
#main


if __name__ == "__main__":
    main()
