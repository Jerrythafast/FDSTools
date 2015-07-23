#!/usr/bin/env python

import re, sys

from ConfigParser import RawConfigParser, MissingSectionHeaderError
from collections import defaultdict

# Patterns that match entire sequences.
PAT_SEQ_RAW = re.compile("^[ACGT]*$")
PAT_SEQ_TSSV = re.compile("^(?:[ACGT]+\(\d+\))*$")
PAT_SEQ_ALLELENAME = re.compile(
    "^(?:(?:(?:CE)?\d+(?:\.\d+)?_(?:[ACGT]+\[\d+\])+)|[XY])"  # n_ACT[m]
    "(?:_[-+]\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # _+3A>
        "(?!(?P=a))(?:[ACTG]+|-))*$")  # Portion of variants after '>'.

# Pattern that matches blocks of TSSV-style sequences and allele names.
PAT_TSSV_BLOCK = re.compile("([ACGT]+)\((\d+)\)")
PAT_ALLELENAME_BLOCK = re.compile("([ACGT]+)\[(\d+)\]")

# Pattern that matches a single prefix/suffix variant.
PAT_VARIANT = re.compile(
    "^([-+]\d+)(?:\.1)?"  # Position number
    "(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # From part
    "(?!(?P=a))([ACTG]+|-)$")  # To part

# Patterns that match (parts of) an STR definition.
PAT_STR_DEF = re.compile(
    "^(?:[ACGT]+\s+\d+\s+\d+(?:\s+[ACGT]+\s+\d+\s+\d+)*)?$")
PAT_STR_DEF_BLOCK = re.compile("([ACGT]+)\s+(\d+)\s+(\d+)")

# Pattern to split a comma-, semicolon-, or space-separated list.
PAT_SPLIT = re.compile("[,; ]+")


def call_variants(template, sequence, reverse_indices=False, cache=True,
                  debug=False):
    """
    Perform a global alignment of sequence to template and return a
    list of variants detected.  All variants are given as substitutions
    in the form posX>Y, where the first base in the template is pos=1.
    Set reverse_indices to True to count from right to left instead.
    Insertions and deletions are pos.1->Y and posX>-, respectively.

    By default, the results of this function are cached.  Set cache to
    False to suppress caching the result and reduce memory usage.

    Setting debug to True will cause the alignment matrices to be
    printed to sys.stdout.  Be aware that they can be quite large.
    """
    # Saving the results in a cache to avoid repeating alignments.
    try:
        return call_variants.cache[template, sequence, reverse_indices]
    except KeyError:
        pass

    row_offset = len(template) + 1
    matrix_match = [0] * row_offset * (len(sequence)+1)
    matrix_gap1 = [-sys.maxint-1] * row_offset * (len(sequence)+1)
    matrix_gap2 = [-sys.maxint-1] * row_offset * (len(sequence)+1)

    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    GAP_OPEN_SCORE = -10
    GAP_EXTEND_SCORE = -1

    for i in range(len(matrix_match)):
        x = i % row_offset
        y = i / row_offset

        # Initialisation of first row and column.
        if x == 0 or y == 0:
            if x != 0:
                # Top row.
                matrix_gap1[i] = GAP_OPEN_SCORE + GAP_EXTEND_SCORE * (x - 1)
                matrix_match[i] = matrix_gap1[i]
            elif y != 0:
                # Left column.
                matrix_gap2[i] = GAP_OPEN_SCORE + GAP_EXTEND_SCORE * (y - 1)
                matrix_match[i] = matrix_gap2[i]
            continue

        matrix_gap1[i] = max(
            matrix_match[i-1] + GAP_OPEN_SCORE,
            matrix_gap1[i-1] + GAP_EXTEND_SCORE)
        matrix_gap2[i] = max(
            matrix_match[i-row_offset] + GAP_OPEN_SCORE,
            matrix_gap2[i-row_offset] + GAP_EXTEND_SCORE)

        if template[x-1] == sequence[y-1]:
            match = MATCH_SCORE
        else:
            match = MISMATCH_SCORE

        matrix_match[i] = max(
            matrix_match[i-1-row_offset] + match,
            matrix_gap1[i],
            matrix_gap2[i])

    if debug:
        print("GAP1")
        for i in range(0, len(matrix_gap1), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap1[i:i+row_offset]))
        print("GAP2")
        for i in range(0, len(matrix_gap2), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_gap2[i:i+row_offset]))
        print("Match")
        for i in range(0, len(matrix_match), row_offset):
            print(("%5i" * row_offset) % tuple(matrix_match[i:i+row_offset]))


    # Backtracking.
    variants = []
    variant_template = 0
    variant_sequence = 0
    i = len(matrix_match) - 1
    while i >= 0:
        x = i % row_offset
        y = i / row_offset

        if matrix_gap1[i] == matrix_match[i]:
            # Go horizontally.  Deletion.
            variant_template += 1
            i -= 1
            continue

        if matrix_gap2[i] == matrix_match[i]:
            # Go vertically.  Insertion.
            variant_sequence += 1
            i -= row_offset
            continue

        # Only backtracking diagonally if a gap is out of the question.
        # Go diagonally.  Either match or mismatch.
        if i == 0 or template[x - 1] == sequence[y - 1]:
            # Match.  Flush variants.
            if variant_template or variant_sequence:
                if variant_template == 0:
                    # Insertions: "-131.1->C" instead of "-130->C".
                    variants.append("%+i.1->%s" % (
                        x - int(reverse_indices) * row_offset,
                        sequence[y:y+variant_sequence]))
                else:
                    variants.append("%+i%s>%s" % (
                        (x + 1) - int(reverse_indices) * row_offset,
                        template[x:x+variant_template],
                        sequence[y:y+variant_sequence] or "-"))
                variant_template = 0
                variant_sequence = 0
        else:
            # Start/extend mismatch.
            variant_template += 1
            variant_sequence += 1
        i -= 1 + row_offset

    # If reverse_indices=False, we need to reverse the output instead.
    if not reverse_indices:
        variants.reverse()

    # Store the result in the cache.
    if cache:
        call_variants.cache[template, sequence, reverse_indices] = variants
    return variants
#call_variants
call_variants.cache = {}


def mutate_sequence(sequence, variants):
    """Apply the given variants to the given sequence."""
    if not sequence and len(variants) > 1:
        raise ValueError("With an empty sequence, only a single variant is "
                         "possible: an insertion '+0.1->x' or '-1.1->x'.")
    sequence = list(sequence)
    for variant in variants:
        vm = PAT_VARIANT.match(variant)
        if vm is None:
            raise ValueError("Unrecognised variant '%s'" % variant)
        pos = int(vm.group(1))
        old = vm.group(2)
        new = vm.group(3)
        if old == "-":
            old = ""
        if new == "-":
            new = ""
        if pos < 0:
            pos += len(sequence) + 1
        if pos == 0 and len(sequence) == 0:
            # Insertion into empty sequence.
            return new
        if old != "".join(sequence[pos-1:pos+len(old)-1]):
            raise ValueError("Incorrect original sequence in variant '%s'; "
                             "should be '%s'!" % (variant,
                               "".join(sequence[pos-1:pos+len(old)-1]) or "-"))
        sequence[pos-1:pos+len(old)-1] = [""] * len(old)
        if pos:
            sequence[pos-1] += new
        else:
            # Insertion at the beginning of the sequence
            sequence[0] = new + sequence[0]
    return "".join(sequence)
#mutate_sequence


def parse_library(handle):
    try:
        return parse_library_ini(handle)
    except MissingSectionHeaderError:
        # Not an ini file.
        pass
    handle.seek(0)
    return parse_library_tsv(handle)
#parse_library


def parse_library_tsv(handle):
    """
    Parse a TSSV library file (tab-separated values format).

    The provided file should contain at least four columns: marker name,
    left flanking sequence, right flanking sequence, and STR definition.

    Return a nested dict with top-level keys "flanks" and "regex".
    """
    library = {
      "flanks": {},
      "regex": {}
    }
    for line in handle:
        line = map(lambda x: x.strip(), line.rstrip("\r\n").split("\t"))
        if line == [""]:
            continue
        if len(line) < 4:
            raise ValueError(
                "Invalid library file: encountered line with %i columns, "
                "need at least 4" % len(line))
        marker = line[0]
        if PAT_SEQ_RAW.match(line[1]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" %
                             (line[1], marker))
        if PAT_SEQ_RAW.match(line[2]) is None:
            raise ValueError("Flanking sequence '%s' of marker %s is invalid" %
                             (line[2], marker))
        if PAT_STR_DEF.match(line[3]) is None:
            raise ValueError("STR definition '%s' of marker %s is invalid" %
                             (line[3], marker))
        library["flanks"][marker] = line[1:3]
        library["regex"][marker] = re.compile("".join(["^"] + map(lambda x:
            "(%s){%s,%s}" % x, PAT_STR_DEF_BLOCK.findall(line[3])) + ["$"]))
    return library
#parse_library_tsv


def parse_library_ini(handle):
    library = {
      "flanks": {},
      "prefix": {},
      "suffix": {},
      "regex": {},
      "length_adjust": defaultdict(int),
      "block_length": defaultdict(lambda: 4),
      "aliases": {}
    }
    markers = set()

    ini = RawConfigParser()
    ini.optionxform = str
    ini.readfp(handle)
    for section in ini.sections():
        for marker in ini.options(section):
            value = ini.get(section, marker)
            if section == "flanks":
                values = PAT_SPLIT.split(value)
                if len(values) != 2:
                    raise ValueError(
                        "For marker %s, %i flanking sequences were given,"
                        "need exactly 2" % (marker, len(values)))
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Flanking sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["flanks"][marker] = values
                markers.add(marker)
            elif section == "prefix":
                values = PAT_SPLIT.split(value)
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Prefix sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["prefix"][marker] = list(set(values))
                markers.add(marker)
            elif section == "suffix":
                values = PAT_SPLIT.split(value)
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Suffix sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["suffix"][marker] = list(set(values))
                markers.add(marker)
            elif section == "length_adjust":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Length adjustment '%s' of marker %s is not a valid "
                        "integer" % (value, marker))
                library["length_adjust"][marker] = value
                markers.add(marker)
            elif section == "block_length":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Block length '%s' of marker %s is not a valid integer"
                        % (value, marker))
                library["block_length"][marker] = value
                markers.add(marker)
            elif section == "aliases":
                values = PAT_SPLIT.split(value)
                if len(values) != 3:
                    raise ValueError("Alias %s does not have 3 values, but %i"
                                     % (marker, len(values)))
                if PAT_SEQ_RAW.match(values[1]) is None:
                        raise ValueError(
                            "Alias sequence '%s' of alias %s is invalid" %
                            (values[1], marker))
                library["aliases"][marker] = {
                    "marker": values[0],
                    "sequence": values[1],
                    "name": values[2]
                }
                markers.add(marker)
            elif section == "repeat":
                if PAT_STR_DEF.match(value) is None:
                    raise ValueError(
                        "STR definition '%s' of marker %s is invalid" %
                        (value, marker))
                library["regex"][marker] = value
                markers.add(marker)

    # Compile regular expressions.
    # NOTE: The libconvert tool expects "(seq){num,num}" blocks ONLY!
    # TODO: Should a single prefix/suffix be required (i.e., seq{1,1})?
    #       Then also update libconvert when converting to TSSV format.
    for marker in markers:
        parts = []
        if marker in library["prefix"]:
            parts += map(lambda x: "(%s){0,1}" % x, library["prefix"][marker])
        if marker in library["aliases"]:
            parts.append("(%s){0,1}" % library["aliases"][marker]["sequence"])
        elif marker in library["regex"]:
            parts += map(lambda x: "(%s){%s,%s}" % x,
                         PAT_STR_DEF_BLOCK.findall(library["regex"][marker]))
        if marker in library["suffix"]:
            parts += map(lambda x: "(%s){0,1}" % x, library["suffix"][marker])
        if parts:
            library["regex"][marker] = re.compile(
                "".join(["^"] + parts + ["$"]))
    return library
#parse_library_ini


def detect_sequence_format(seq):
    """Return format of seq.  One of 'raw', 'tssv', or 'allelename'."""
    if not seq:
        raise ValueError("Empty sequence")
    if PAT_SEQ_RAW.match(seq):
        return 'raw'
    if PAT_SEQ_TSSV.match(seq):
        return 'tssv'
    if PAT_SEQ_ALLELENAME.match(seq):
        return 'allelename'
    raise ValueError("Unrecognised sequence format")
#detect_sequence_format


def convert_sequence_tssv_raw(seq):
    return "".join(map(
        lambda block: block[0] * int(block[1]),
        PAT_TSSV_BLOCK.findall(seq)))
#convert_sequence_tssv_raw


def convert_sequence_raw_tssv(seq, library, marker, return_alias=False):
    # Try to match this marker's pattern, or any of its aliases.
    match = None
    if "aliases" in library:
        for alias in library["aliases"]:
            if (library["aliases"][alias]["marker"] == marker and
                    alias in library["regex"]):
                match = library["regex"][alias].match(seq)
                if match is not None:
                    marker = alias
                    break
    if match is None and marker in library["regex"]:
        match = library["regex"][marker].match(seq)

    if match is None:
        # TODO: still try prefix(1)middle(1)suffix(1) for this marker!
        seq = "%s(1)" % seq
    else:
        seq = reduce(
            lambda a, b:
                (a[0] + "%s(%i)" % (b[0], (b[1]-a[1])/len(b[0])),
                 b[1]),
            reduce(
                lambda x, y:
                    x if y == (None, -1) else
                    x[:-1] + [y] if x[-1][0] == y[0] else
                    x + [y],
                map(
                    lambda z: (match.group(z), match.regs[z][1]),
                    range(1, len(match.regs))),
                [("", 0)]))[0]
    return (seq, marker) if return_alias else seq
#convert_sequence_raw_tssv


def convert_sequence_allelename_tssv(seq, library, marker):
    # Check whether there is an alias for this sequence.
    if "aliases" in library:
        for alias in library["aliases"]:
            if library["aliases"][alias]["marker"] == marker and (
                    seq == library["aliases"][alias]["name"] or
                    seq.startswith(library["aliases"][alias]["name"] + "_")):
                marker = alias
                seq = "".join([
                    "0_",
                    library["aliases"][alias]["sequence"] + "[1]",
                    seq[len(library["aliases"][alias]["name"]):]])
                break

    allele = seq.split("_")

    # Get and mutate prefix and suffix.
    if "prefix" in library and marker in library["prefix"]:
        prefix = library["prefix"][marker][0]
    else:
        prefix = ""
    if "suffix" in library and marker in library["suffix"]:
        suffix = library["suffix"][marker][0]
    else:
        suffix = ""
    variants = [[], []]
    for variant in allele[2:]:
        if variant[0] == "-":
            if not prefix:
                raise ValueError("Encountered prefix variant '%s', but marker "
                                 "'%s' has no prefix!" % (variant, marker))
            variants[0].append(variant)
        elif variant[0] == "+":
            if not suffix:
                raise ValueError("Encountered suffix variant '%s', but marker "
                                 "'%s' has no suffix!" % (variant, marker))
            variants[1].append(variant)
        else:
            raise ValueError("Unrecognised variant '%s'" % variant)
    if variants[0]:
        prefix = mutate_sequence(prefix, variants[0])
    if variants[1]:
        suffix = mutate_sequence(suffix, variants[1])

    blocks = []
    if prefix:
        blocks.append((prefix, 1))
    for block in PAT_ALLELENAME_BLOCK.findall(allele[1]):
        blocks.append((block[0], int(block[1])))
    if suffix:
        blocks.append((suffix, 1))
    return "".join(map(lambda block: "%s(%i)" % block, blocks))
#convert_sequence_allelename_tssv


def convert_sequence_raw_allelename(seq, library, marker):
    # We actually convert raw->allelename via TSSV format.
    seq, alias = convert_sequence_raw_tssv(seq, library, marker, True)
    blocks = PAT_TSSV_BLOCK.findall(
        convert_sequence_raw_tssv(seq, library, marker))

    # Generate prefix/suffix variants.
    # TODO: Handle missing/unrecognised prefix/suffix gracefully.
    length = 0
    variants = []
    try:
        if "prefix" in library and marker in library["prefix"]:
            prefix = library["prefix"][marker][0]
            if prefix:
                if blocks[0][1] != "1":
                    raise ValueError("Repeated prefix")
                if prefix != blocks[0][0]:
                    variants += call_variants(prefix, blocks[0][0], True)
                length += len(blocks[0][0]) - len(prefix)
                blocks = blocks[1:]
        if "suffix" in library and marker in library["suffix"]:
            suffix = library["suffix"][marker][0]
            if suffix:
                if blocks[-1][1] != "1":
                    raise ValueError("Repeated suffix")
                if suffix != blocks[-1][0]:
                    variants += call_variants(suffix, blocks[-1][0], False)
                length += len(blocks[-1][0]) - len(suffix)
                blocks = blocks[:-1]
    except IndexError:
        raise ValueError("Missing prefix/suffix in '%s' allele '%s'!" %
            (marker, seq))

    # We are ready to return the allele name of aliases.
    if alias != marker:
        return "_".join([library["aliases"][marker]["name"]] + variants)

    # Compute CE allele number for the other alleles.
    # TODO: perhaps produce a more intelligent name if there is exactly
    #       one alias with the same length
    blocknames = []
    if "block_length" in library:
        blocksize = library["block_length"][marker]
    else:
        blocksize = 4
    if "length_adjust" in library:
        length -= library["length_adjust"][marker]
    for block in blocks:
        blocknames.append("%s[%s]" % (block[0], block[1]))
        length += len(block[0]) * int(block[1])

    allelename = "CE" + str(length / blocksize)
    if length % blocksize:
        allelename += "." + str(length % blocksize)
    return "_".join([allelename, "".join(blocknames)] + variants)
#convert_sequence_raw_allelename


def ensure_sequence_format(seq, to_format, from_format=None, library=None,
                           marker=None):
    """Convert seq to 'raw', 'tssv', or 'allelename' format."""
    known_formats = ["raw", "tssv", "allelename"]
    if to_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" %
                         (to_format, known_formats))
    if from_format is None:
        from_format = detect_sequence_format(seq) if seq else "raw"
    elif from_format not in known_formats:
        raise ValueError("Unknown format '%s', choose from %s" %
                         (from_format, known_formats))

    # No conversion needed?
    if to_format == from_format:
        return seq

    # From TSSV to raw sequence is easy.
    # We'll need a library and marker name for anything else.
    if (library is None or marker is None) and (from_format != "tssv" or
            to_format != "raw"):
        raise ValueError("Sequence needs to be converted from %s to %s, this "
                          "conversion requires a library file" %
                          (from_format, to_format))

    # Perform conversions.
    if from_format == "allelename":
        seq = convert_sequence_allelename_tssv(seq, library, marker)
    if to_format == "tssv":
        if from_format == "raw":
            return convert_sequence_raw_tssv(seq, library, marker)
        return seq
    if from_format != "raw":
        seq = convert_sequence_tssv_raw(seq)
    if to_format == "raw":
        return seq
    return convert_sequence_raw_allelename(seq, library, marker)
#ensure_sequence_format


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
