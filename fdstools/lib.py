#!/usr/bin/env python

import re

# Patterns that match entire sequences.
PAT_SEQ_RAW = re.compile("^[ACGT]*$")
PAT_SEQ_TSSV = re.compile("^(?:[ACGT]+\(\d+\))*$")
PAT_SEQ_ALLELENAME = re.compile(
    "^(?:(?:\d+(?:\.\d+)?_(?:[ACGT]+\[\d+\])+)|[XY])"  # nn_ACGT[nn]...
    "(?:_[-+]\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # _+3A>
        "(?!(?P=a))(?:[ACTG]+|-))*$")  # Portion of variants after '>'.

# Pattern that matches blocks of TSSV-style sequences.
PAT_TSSV_BLOCK = re.compile("([ACTG]+)\((\d+)\)")


def detect_sequence_format(seq):
    """Return format of seq.  One of 'raw', 'tssv', or 'allelename'."""
    if not seq:
        raise ValueError("Empty sequence.")
    if PAT_SEQ_RAW.match(seq):
        return 'raw'
    if PAT_SEQ_TSSV.match(seq):
        return 'tssv'
    if PAT_SEQ_ALLELENAME.match(seq):
        return 'allelename'
    raise ValueError("Unrecognised sequence format.")
#detect_sequence_format


def get_column_ids(column_names, *names):
    """Find all names in column_names and return their indices."""
    result = []
    for name in names:
        try:
            result.append(column_names.index(name))
        except ValueError:
            raise Exception("Column not found in input file: %s" % name)
    if len(result) == 1:
        return result[0]
    return tuple(result)
#get_column_ids


def pos_int_arg(value):
    """Convert str to int, raise ArgumentTypeError if not positive."""
    if not value.isdigit() or not int(value):
        raise argparse.ArgumentTypeError(
            "invalid positive int value: '%s'" % value)
    return int(value)
#pos_int_arg


def print_db(text, debug):
    """Print text if debug is True."""
    if debug:
        print(text)
#print_db
