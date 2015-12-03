#!/usr/bin/env python

import re, sys, argparse, random
#import numpy as np  # Imported only when calling nnls()

from ConfigParser import RawConfigParser, MissingSectionHeaderError
from StringIO import StringIO

# Patterns that match entire sequences.
PAT_SEQ_RAW = re.compile("^[ACGT]*$")
PAT_SEQ_TSSV = re.compile("^(?:[ACGT]+\(\d+\))*$")
PAT_SEQ_ALLELENAME_STR = re.compile(  # First line: n_ACT[m] or alias.
    "^(?:(?:(?:CE)?-?\d+(?:\.\d+)?_(?:[ACGT]+\[\d+\])*)|((?!_).+?))"
    "(?:_[-+]\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"  # _+3A>
        "(?!(?P=a))(?:[ACGT]+|-))*$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_SNP = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "\d+(?:\.1)?(?P<a>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
        "(?!(?P=a))(?:[ACGT]+|-))+$")  # Portion of variants after '>'.
PAT_SEQ_ALLELENAME_MT = re.compile(
    "^REF$|^(?:(?:(?<=^)|(?<!^) )"  # 'REF' or space-separated variants.
    "(?:-?\d+\.\d+[ACGT]|(?P<a>[ACGT])?\d+(?(a)(?!(?P=a)))(?:[ACGT-]|del)))+$")

# Patterns that match blocks of TSSV-style sequences and allele names.
PAT_TSSV_BLOCK = re.compile("([ACGT]+)\((\d+)\)")
PAT_ALLELENAME_BLOCK = re.compile("([ACGT]+)\[(\d+)\]")
PAT_ALIAS = re.compile("^(?!_).+$")

# Patterns that match a single variant.
PAT_VARIANT_STR = re.compile(
    "^(?P<pos>[-+]\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_SNP = re.compile(
    "^(?P<pos>\d+)(?:\.(?P<ins>1))?"
    "(?P<old>(?:(?<=\.1)-)|(?<!\.1)[ACGT]+)>"
    "(?!(?P=old))(?P<new>[ACGT]+|-)$")
PAT_VARIANT_MT = re.compile(
    "^(?P<old>(?P<a>[ACGT])|-?)"
    "(?P<pos>\d+)(?(a)|(?:\.(?P<ins>\d+))?)"
    "(?P<new>[ACGT-]|del)$")

# Patterns that match (parts of) an STR definition.
PAT_STR_DEF = re.compile("^(?:(?:(?<=^)|(?<!^)\s+)[ACGT]+\s+\d+\s+\d+)*$")
PAT_STR_DEF_BLOCK = re.compile("([ACGT]+)\s+(\d+)\s+(\d+)")

# Pattern to split a comma-, semicolon-, or space-separated list.
PAT_SPLIT = re.compile("\s*[,; \t]\s*")

# Pattern that matches a chromosome name/number.
PAT_CHROMOSOME = re.compile(
    "^(?:[Cc][Hh][Rr](?:[Oo][Mm])?)?([1-9XYM]|1\d|2[0-2])$")

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

# IUPAC Table of complementary bases.
COMPL = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G", "R": "Y", "Y": "R",
         "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D",
         "a": "t", "t": "a", "u": "a", "g": "c", "c": "g", "r": "y", "y": "r",
         "k": "m", "m": "k", "b": "v", "v": "b", "d": "h", "h": "d"}


def get_genome_pos(location, x, invert=False):
    """Get the genome position of the x-th base in a sequence."""
    if invert:
        offset = 0
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif pos <= x <= location[i]:
                # x is in the current range
                break
            else:
                offset += location[i]-pos+1
        else:
            if len(location) % 2:
                raise ValueError("Position %i is outside sequence range" % x)
        return offset + x - pos
    else:
        for i in range(1, len(location)):
            if i % 2:
                # Starting position.
                pos = location[i]
            elif location[i]-pos < x:
                # x is after this ending position
                x -= location[i]-pos+1
            else:
                # x is before this ending position
                break
        return pos + x
#get_genome_pos


def call_variants(template, sequence, location="suffix", cache=True,
                  debug=False):
    """
    Perform a global alignment of sequence to template and return a
    list of variants detected.  The format (nomenclature) of the
    returned variants depends on the location argument.

    If location is "suffix" (the default), all variants are given as
    substitutions in the form posX>Y, where the first base in the
    template is pos=1.  With location set to "prefix", bases are counted
    from right to left instead.  Insertions and deletions are written as
    pos.1->Y and posX>-, respectively.

    If location is a tuple ("M", position) with any integer for the
    position, variants are written following the mtDNA nomenclature
    guidelines.  The given position is that of the first base in the
    template.

    If location is a tuple ("chromosome name", position), a
    NotImplementedError is raised.

    By default, the results of this function are cached.  Set cache to
    False to suppress caching the result and reduce memory usage.

    Setting debug to True will cause the alignment matrices to be
    printed to sys.stdout.  Be aware that they can be quite large.
    """
    # Saving the results in a cache to avoid repeating alignments.
    try:
        return call_variants.cache[template, sequence, location]
    except KeyError:
        cache_key = location

    row_offset = len(template) + 1
    matrix_match = [0] * row_offset * (len(sequence)+1)
    matrix_gap1 = [-sys.maxint-1] * row_offset * (len(sequence)+1)
    matrix_gap2 = [-sys.maxint-1] * row_offset * (len(sequence)+1)

    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    GAP_OPEN_SCORE = -10
    GAP_EXTEND_SCORE = -1
    variant_format = "%i%s>%s"

    if location == "prefix":
        location = ("prefix", -len(template))
    elif location == "suffix":
        # Include plus signs for position numbers.
        variant_format = "%+i%s>%s"
        location = ("suffix", 1)
    elif type(location) != tuple or len(location) < 2:
        raise ValueError("Unknown location %r. It should be 'prefix', "
            "'suffix', or a tuple (chromosome, position [, endpos])" %
            location)
    elif location[0] == "M":
        # No need to avoid gaps in mtDNA notation.
        GAP_OPEN_SCORE = -1

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
                if location[0] == "M":
                    # MtDNA variants are one-base-at-a-time.
                    for j in range(
                            max(variant_template, variant_sequence)-1, -1, -1):
                        variants.append("%s%i%s%s" % (
                            template[x+j] if j < variant_template else "",#"-",
                            get_genome_pos(
                                location, x + min(j, variant_template-1)),
                            ".%i" % (j-variant_template+1)
                                if j >= variant_template else "",
                            sequence[y+j] if j < variant_sequence else "del"))
                elif variant_template == 0:
                    # Insertions: "-131.1->C" instead of "-130->C".
                    variants.append(variant_format % (
                        get_genome_pos(location, x - 1),
                        ".1-",
                        sequence[y:y+variant_sequence]))
                else:
                    variants.append(variant_format % (
                        get_genome_pos(location, x),
                        template[x:x+variant_template],
                        sequence[y:y+variant_sequence] or "-"))
                variant_template = 0
                variant_sequence = 0
        else:
            # Start/extend mismatch.
            variant_template += 1
            variant_sequence += 1
        i -= 1 + row_offset

    # Variants were called from right to left.  Reverse their order.
    if location[0] != "prefix":
        variants.reverse()

    # Store the result in the cache.
    if cache:
        call_variants.cache[template, sequence, cache_key] = variants
    return variants
#call_variants
call_variants.cache = {}


def mutate_sequence(seq, variants, location=None):
    """Apply the given variants to the given sequence."""
    if type(location) != tuple or len(location) < 2:
        pattern = PAT_VARIANT_STR
        location = (None, 0)
    elif location[0] == "M":
        pattern = PAT_VARIANT_MT
        location = (location[0], location[1]-1) + tuple(location[2:])
    else:
        pattern = PAT_VARIANT_SNP
        location = (location[0], location[1]-1) + tuple(location[2:])

    seq = [[]] + [[base] for base in seq]
    for variant in variants:
        vm = pattern.match(variant)
        if vm is None:
            raise ValueError("Unrecognised variant '%s'" % variant)
        pos = int(vm.group("pos"))
        ins = int(vm.group("ins") or 0)
        old = vm.group("old")
        new = vm.group("new")
        if old == "-":
            old = ""
        if new == "-" or new == "del":
            new = ""
        if pos < 0:
            pos += len(seq) + 1
        pos = get_genome_pos(location, pos, True)
        if pos < 0 or (pos == 0 and not ins) or pos >= len(seq):
            raise ValueError(
                "Position of variant '%s' is outside sequence range" %
                    (variant))
        if (not ins and old and old != "".join("".join(x[:1])
                for x in seq[pos:pos+len(old)])):
            raise ValueError(
                "Incorrect original sequence in variant '%s'; should be '%s'!"
                % (variant, "".join("".join(x[:1])
                    for x in seq[pos:pos+len(old)])))
        elif not ins and not old:
            # MtDNA substitution with reference base omitted.
            old = "".join("".join(x[:1]) for x in seq[pos:pos+len(new)])
        if not ins:
            # Remove old bases, retaining those inserted between/after.
            seq[pos:pos+len(old)] = [
                [""] + x[1:] for x in seq[pos:pos+len(old)]]
            # Place new entirely in the position of the first old base.
            seq[pos][0] = new
        else:
            # Insert new exactly ins positions after pos.
            while len(seq[pos]) <= ins:
                seq[pos].append("")
            seq[pos][ins] = new
    return "".join("".join(x) for x in seq)
#mutate_sequence


def parse_library(libfile, stream=False):
    try:
        if not stream:
            libfile = sys.stdin if libfile == "-" else open(libfile, "r")
        if libfile == sys.stdin:
            # Can't seek on pipes, so read it into a buffer first.
            libfile = StringIO(sys.stdin.read())
        try:
            library = parse_library_ini(libfile)
            if not stream:
                libfile.close()
            return library
        except MissingSectionHeaderError:
            # Not an ini file.
            pass
        libfile.seek(0)
        library = parse_library_tsv(libfile)
        if not stream and libfile != sys.stdin:
            libfile.close()
        return library
    except ValueError as err:
        raise argparse.ArgumentTypeError(err)
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
      "regex": {},
      "regex_middle": {}
    }
    for line in handle:
        line = [x.strip() for x in line.rstrip("\r\n").split("\t")]
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
        library["regex_middle"][marker] = re.compile("".join(
            "(%s){%s,%s}" % x for x in PAT_STR_DEF_BLOCK.findall(line[3])))
        library["regex"][marker] = re.compile(
            "".join(["^", library["regex_middle"][marker].pattern, "$"]))
    return library
#parse_library_tsv


def parse_library_ini(handle):
    library = {
      "flanks": {},
      "prefix": {},
      "suffix": {},
      "regex": {},
      "regex_middle": {},
      "nostr_reference": {},
      "genome_position": {},
      "length_adjust": {},
      "block_length": {},
      "max_expected_copies": {},
      "aliases": {}
    }
    markers = set()

    ini = RawConfigParser()
    ini.optionxform = str
    ini.readfp(handle)
    for section in ini.sections():
        for marker in ini.options(section):
            value = ini.get(section, marker)
            section_low = section.lower()
            if section_low == "flanks":
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
            elif section_low == "prefix":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "A prefix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Prefix sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["prefix"][marker] = values
                markers.add(marker)
            elif section_low == "suffix":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "A suffix was defined for non-STR marker %s" % marker)
                values = PAT_SPLIT.split(value)
                for value in values:
                    if PAT_SEQ_RAW.match(value) is None:
                        raise ValueError(
                            "Suffix sequence '%s' of marker %s is invalid" %
                            (value, marker))
                library["suffix"][marker] = values
                markers.add(marker)
            elif section_low == "genome_position":
                values = PAT_SPLIT.split(value)
                chromosome = PAT_CHROMOSOME.match(values[0])
                if chromosome is None:
                    raise ValueError(
                        "Invalid chromosome '%s' for marker %s." %
                        (values[0], marker))
                pos = [chromosome.group(1)]
                for i in range(1, len(values)):
                    try:
                        pos.append(int(values[i]))
                    except:
                        raise ValueError(
                            "Position '%s' of marker %s is not a valid integer"
                            % (values[i], marker))
                    if not i % 2 and pos[-2] >= pos[-1]:
                        raise ValueError(
                            "End position %i of marker %s must be higher than "
                            "corresponding start position %i" %
                            (pos[-1], marker, pos[-2]))
                if len(values) == 1:
                    pos.append(1)
                library["genome_position"][marker] = tuple(pos)
                markers.add(marker)
            elif section_low == "length_adjust":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Length adjustment '%s' of marker %s is not a valid "
                        "integer" % (value, marker))
                library["length_adjust"][marker] = value
                markers.add(marker)
            elif section_low == "block_length":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Block length '%s' of marker %s is not a valid integer"
                        % (value, marker))
                library["block_length"][marker] = value
                markers.add(marker)
            elif section_low == "max_expected_copies":
                try:
                    value = int(value)
                except:
                    raise ValueError(
                        "Maximum number of expected copies '%s' of marker %s "
                        "is not a valid integer" % (value, marker))
                library["max_expected_copies"][marker] = value
                markers.add(marker)
            elif section_low == "aliases":
                values = PAT_SPLIT.split(value)
                if len(values) != 3:
                    raise ValueError("Alias %s does not have 3 values, but %i"
                                     % (marker, len(values)))
                if PAT_SEQ_RAW.match(values[1]) is None:
                    raise ValueError(
                        "Alias sequence '%s' of alias %s is invalid" %
                        (values[1], marker))
                if PAT_ALIAS.match(values[2]) is None:
                    raise ValueError(
                        "Allele name '%s' of alias %s is invalid" %
                        (values[2], marker))
                library["aliases"][marker] = {
                    "marker": values[0],
                    "sequence": values[1],
                    "name": values[2]
                }
                markers.add(marker)
            elif section_low == "repeat":
                if marker in library["nostr_reference"]:
                    raise ValueError(
                        "Marker %s was encountered in both [repeat] and "
                        "[no_repeat] sections" % marker)
                if PAT_STR_DEF.match(value) is None:
                    raise ValueError(
                        "STR definition '%s' of marker %s is invalid" %
                        (value, marker))
                library["regex"][marker] = value
                markers.add(marker)
            elif section_low == "no_repeat":
                if marker in library["regex"]:
                    raise ValueError(
                        "Marker %s was encountered in both [repeat] and "
                        "[no_repeat] sections" % marker)
                if marker in library["prefix"] or marker in library["suffix"]:
                    raise ValueError(
                        "A prefix or suffix was defined for non-STR marker %s"
                        % marker)
                if PAT_SEQ_RAW.match(value) is None:
                    raise ValueError(
                        "Reference sequence '%s' of marker %s is invalid" %
                        (value, marker))
                library["nostr_reference"][marker] = value
                markers.add(marker)

    # Sanity check: prohibit prefix/suffix for aliases of non-STRs.
    for alias in library["aliases"]:
        if library["aliases"][alias]["marker"] in library["nostr_reference"] \
                and (alias in library["prefix"] or alias in library["suffix"]):
            raise ValueError(
                "A prefix or suffix was defined for alias %s of non-STR "
                "marker %s" % (alias, library["aliases"][alias]["marker"]))

    # Sanity check: end position of marker should reflect ref length.
    for marker in library["genome_position"]:
        if marker not in library["nostr_reference"]:
            continue
        pos = library["genome_position"][marker]
        reflength = len(library["nostr_reference"][marker])
        length = 0
        for i in range(2, len(pos), 2):
            length += pos[i] - pos[i-1] + 1
        if reflength < length or (len(pos) % 2 and reflength != length):
            raise ValueError(
                "Length of reference sequence of marker %s is %i bases, but "
                "genome positions add up to %i bases" %
                (marker, reflength, length))

    # Compile regular expressions.
    # NOTE: The libconvert tool expects "(seq){num,num}" blocks ONLY!
    # TODO: Should a single prefix/suffix be required (i.e., seq{1,1})?
    #       Then also update libconvert when converting to TSSV format.
    for marker in markers:
        parts = []
        partsm = []
        if marker in library["prefix"]:
            parts += ("(%s){0,1}" % x for x in library["prefix"][marker])
        if marker in library["aliases"]:
            parts.append("(%s){0,1}" % library["aliases"][marker]["sequence"])
            partsm.append("(%s){0,1}" % library["aliases"][marker]["sequence"])
        elif marker in library["regex"]:
            partsm = ["(%s){%s,%s}" % x for x in
                      PAT_STR_DEF_BLOCK.findall(library["regex"][marker])]
            parts += partsm
        if marker in library["suffix"]:
            parts += ("(%s){0,1}" % x for x in library["suffix"][marker])
        if parts:
            library["regex"][marker] = re.compile(
                "".join(["^"] + parts + ["$"]))
        if partsm:
            library["regex_middle"][marker] = re.compile("".join(partsm))
    return library
#parse_library_ini


def load_profiles(profilefile, library=None):
    column_names = profilefile.readline().rstrip("\r\n").split("\t")
    (colid_marker, colid_allele, colid_sequence, colid_fmean, colid_rmean,
     colid_tool) = get_column_ids(column_names, "marker", "allele", "sequence",
        "fmean", "rmean", "tool")

    profiles = {}
    for line in profilefile:
        line = line.rstrip("\r\n").split("\t")
        if line == [""]:
            continue
        marker = line[colid_marker]
        if marker not in profiles:
            profiles[marker] = {
                "m": set(),  # To be replaced by its length below.
                "n": set(),  # To be replaced by its length below.
                "seqs": [],
                "forward": {},  # To be replaced by a list below.
                "reverse": {},  # To be replaced by a list below.
                "tool": {}  # To be replaced by a list below.
                }
        allele = ensure_sequence_format(line[colid_allele], "raw",
            library=library, marker=marker)
        sequence = ensure_sequence_format(line[colid_sequence], "raw",
            library=library, marker=marker)
        if (allele, sequence) in profiles[marker]["forward"]:
            raise ValueError(
                "Invalid background noise profiles file: encountered "
                "multiple values for marker '%s' allele '%s' sequence '%s'" %
                (marker, allele, sequence))
        profiles[marker]["forward"][allele,sequence] = float(line[colid_fmean])
        profiles[marker]["reverse"][allele,sequence] = float(line[colid_rmean])
        profiles[marker]["tool"][allele, sequence] = line[colid_tool]
        profiles[marker]["m"].update((allele, sequence))
        profiles[marker]["n"].add(allele)

    # Check completeness and reorder true alleles.
    for marker in profiles:
        profiles[marker]["seqs"] = list(profiles[marker]["n"]) + \
            list(profiles[marker]["m"]-profiles[marker]["n"])
        profiles[marker]["n"] = len(profiles[marker]["n"])
        profiles[marker]["m"] = len(profiles[marker]["m"])
        newprofiles = {"forward": [], "reverse": []}
        tools = []
        for i in range(profiles[marker]["n"]):
            allele = profiles[marker]["seqs"][i]
            for direction in newprofiles:
                newprofiles[direction].append([0] * profiles[marker]["m"])
            tools.append([""] * profiles[marker]["m"])
            for j in range(profiles[marker]["m"]):
                sequence = profiles[marker]["seqs"][j]
                if (allele, sequence) in profiles[marker]["forward"]:
                    for direction in newprofiles:
                        newprofiles[direction][i][j] = \
                            profiles[marker][direction][allele, sequence]
                    tools[i][j] = profiles[marker]["tool"][allele, sequence]
        profiles[marker]["forward"] = newprofiles["forward"]
        profiles[marker]["reverse"] = newprofiles["reverse"]
        profiles[marker]["tool"] = tools

    return profiles
#load_profiles


def regex_longest_match(pattern, subject):
    """Return the longest match of the pattern in the subject string."""
    match = None
    pos = 0
    while pos < len(subject):
        m = pattern.search(subject, pos)
        if m is None:
            break
        if match is None or m.end()-m.start() > match.end()-match.start():
            match = m
        pos = m.start() + 1
    return match
#regex_longest_match


def detect_sequence_format(seq):
    """Return format of seq.  One of 'raw', 'tssv', or 'allelename'."""
    if not seq:
        raise ValueError("Empty sequence")
    if PAT_SEQ_RAW.match(seq):
        return 'raw'
    if PAT_SEQ_TSSV.match(seq):
        return 'tssv'
    if PAT_SEQ_ALLELENAME_STR.match(seq) or PAT_SEQ_ALLELENAME_MT.match(seq) \
            or PAT_SEQ_ALLELENAME_SNP.match(seq):
        return 'allelename'
    raise ValueError("Unrecognised sequence format")
#detect_sequence_format


def convert_sequence_tssv_raw(seq):
    return "".join(block[0] * int(block[1])
                   for block in PAT_TSSV_BLOCK.findall(seq))
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
    if match is not None:
        parts = ((match.group(i), match.end(i)) for i in range(1, 1 if
            match.lastindex is None else match.lastindex+1) if match.group(i))

    # Use heuristics if the sequence does not match the pattern.
    else:

        # Find explictily provided prefix and/or suffix if present.
        pre_suf = ["", ""]
        if "prefix" in library and marker in library["prefix"]:
            for prefix in library["prefix"][marker]:
                if seq.startswith(prefix):
                    pre_suf[0] = prefix
                    seq = seq[len(prefix):]
                    break
        if "suffix" in library and marker in library["suffix"]:
            for suffix in library["suffix"][marker]:
                if seq.endswith(suffix):
                    pre_suf[1] = suffix
                    seq = seq[:-len(suffix)]
                    break

        # Find longest match of middle pattern.
        middle = [(seq, len(pre_suf[0])+len(seq))] if seq else []
        if middle and marker in library["regex_middle"]:
            match = regex_longest_match(library["regex_middle"][marker], seq)
            if match is not None and match.end()-match.start():

                # If this allele does not match the prefix of this
                # marker, but the canonical prefix of the marker ends
                # with the same sequence as the start of our match, we
                # move that portion of the match into the prefix.
                # Then, we do the same thing with the suffix.
                matched = match.group()
                start = match.start()
                end = match.end()
                modified = False
                if (not pre_suf[0] and "prefix" in library
                        and marker in library["prefix"]):
                    ref = library["prefix"][marker][0]
                    i = min(len(ref), len(matched))
                    while i > 0:
                        if ref.endswith(matched[:i]):
                            start += i
                            matched = matched[i:]
                            modified = True
                            break
                        i -= 1
                if (not pre_suf[1] and "suffix" in library
                        and marker in library["suffix"]):
                    ref = library["suffix"][marker][0]
                    i = min(len(ref), len(matched))
                    while i > 0:
                        if ref.startswith(matched[-i:]):
                            end -= i
                            matched = matched[:-i]
                            modified = True
                            break
                        i -= 1
                if modified:
                    from_start = start-match.start()
                    from_end = match.end()-end
                    middle = reduce(
                        lambda x, i: (
                            x[0] + [match.group(i)] *
                                ((match.end(i)-x[1])/len(match.group(i))),
                            match.end(i)) if match.group(i) else x,
                        range(1, match.lastindex+1), ([], match.start()))[0]
                    while from_start:
                        if from_start < len(middle[0]):
                            middle[0] = middle[0][from_start:]
                            break
                        else:
                            from_start -= len(middle[0])
                            middle = middle[1:]
                    while from_end:
                        if from_end < len(middle[-1]):
                            middle[-1] = middle[-1][:-from_end]
                            break
                        else:
                            from_end -= len(middle[-1])
                            middle = middle[:-1]
                    if middle:
                        middle = reduce(
                            lambda x, y: (x[:-1] if x[-1][0] == y else x) +
                                [(y, x[-1][1]+len(y))], middle[1:],
                                [(middle[0],
                                  start+len(middle[0])+len(pre_suf[0]))])

                else:
                    # No trickery with prefix or suffix was done.
                    middle = [(match.group(i), match.end(i)+len(pre_suf[0]))
                        for i in range(1, match.lastindex+1) if match.group(i)]

                pre_suf[0] += seq[:start]
                pre_suf[1] = seq[end:] + pre_suf[1]
                seq = matched

        # Now construct parts.
        parts = []
        if pre_suf[0]:
            parts.append((pre_suf[0], len(pre_suf[0])))
        parts += middle
        if pre_suf[1]:
            parts.append((pre_suf[1], sum(map(len,pre_suf))+len(seq)))

    seq = reduce(
        lambda a, b: (a[0] + "%s(%i)" % (b[0], (b[1]-a[1])/len(b[0])), b[1]),
        reduce(
            lambda x, y: x[:-1] + [y] if x[-1][0] == y[0] else x + [y],
            parts,
            [("", 0)]))[0]
    return (seq, marker) if return_alias else seq
#convert_sequence_raw_tssv


def convert_sequence_allelename_tssv(seq, library, marker):
    # Check whether there is an alias for this sequence.
    alias_of = None
    if "aliases" in library:
        for alias in library["aliases"]:
            if library["aliases"][alias]["marker"] == marker and (
                    seq == library["aliases"][alias]["name"] or
                    seq.startswith(library["aliases"][alias]["name"] + "_")):
                alias_of = marker
                marker = alias
                seq = "".join([
                    "0_",
                    library["aliases"][alias]["sequence"] + "[1]",
                    seq[len(library["aliases"][alias]["name"]):]])
                break

    nameformat = None
    if PAT_SEQ_ALLELENAME_MT.match(seq) is not None:
        nameformat = "MtDNA"
    elif PAT_SEQ_ALLELENAME_SNP.match(seq) is not None:
        nameformat = "SNP"
    if nameformat is not None:
        # MtDNA and SNP markers.
        try:
            reference = library["nostr_reference"][marker]
        except KeyError:
            raise ValueError(
                "%s allele '%s' found for marker %s, but "
                "no reference sequence was found in the library" %
                (nameformat, seq, marker))
        if seq == "REF":
            return reference + "(1)"
        return mutate_sequence(reference, seq.split(),
            library["genome_position"].get(marker,
                ("M" if nameformat == "MtDNA" else "", 1))) + "(1)"

    # Note: aliases of mtDNA and SNP markers end up here as well.
    # It should NOT look like an alias now, however.
    match = PAT_SEQ_ALLELENAME_STR.match(seq)
    if match is None or match.group(1) is not None:
        raise ValueError("Invalid allele name '%s' encountered!" % seq)

    allele = seq.split("_")

    # Get and mutate prefix and suffix.
    prefix = ""
    suffix = ""
    if "prefix" in library:
        if marker in library["prefix"]:
            prefix = library["prefix"][marker][0]
        elif alias_of is not None and alias_of in library["prefix"]:
            prefix = library["prefix"][alias_of][0]
    if "suffix" in library:
        if marker in library["suffix"]:
            suffix = library["suffix"][marker][0]
        elif alias_of is not None and alias_of in library["suffix"]:
            suffix = library["suffix"][alias_of][0]
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
    return "".join("%s(%i)" % block for block in blocks)
#convert_sequence_allelename_tssv


def convert_sequence_raw_allelename(seq, library, marker):
    # We actually convert raw->allelename via TSSV format.
    seq, alias = convert_sequence_raw_tssv(seq, library, marker, True)
    blocks = PAT_TSSV_BLOCK.findall(seq)

    if marker in library["nostr_reference"]:
        # Handle non-STR markers here.
        if alias != marker:
            return library["aliases"][alias]["name"]
        if not blocks:
            # Oh dear, empty sequence... Primer dimer?
            blocks = (("",),)
        if library["nostr_reference"][marker] == blocks[0][0]:
            return "REF"
        return " ".join(
            call_variants(library["nostr_reference"][marker], blocks[0][0],
                library["genome_position"].get(marker, "suffix")))

    # Find prefix and suffix.
    prefix = suffix = this_prefix = this_suffix = ""
    remaining_blocks = len(blocks)
    if "prefix" in library:
        if alias in library["prefix"]:
            prefix = library["prefix"][alias][0]
        elif marker in library["prefix"]:
            prefix = library["prefix"][marker][0]
        if prefix and remaining_blocks > 0 and blocks[0][1] == "1":
            remaining_blocks -= 1
    if "suffix" in library:
        if alias in library["suffix"]:
            suffix = library["suffix"][alias][0]
        elif marker in library["suffix"]:
            suffix = library["suffix"][marker][0]
        if suffix and remaining_blocks > 0 and blocks[-1][1] == "1":
            remaining_blocks -= 1
    if remaining_blocks > 0 and prefix and blocks[0][1] == "1":
        this_prefix = blocks[0][0]
        blocks = blocks[1:]
    if remaining_blocks > 0 and suffix and blocks[-1][1] == "1":
        this_suffix = blocks[-1][0]
        blocks = blocks[:-1]

    # Generate prefix/suffix variants.
    length = 0
    variants = []
    if prefix != this_prefix:
        variants += call_variants(prefix, this_prefix, "prefix")
        length += len(this_prefix) - len(prefix)
    if suffix != this_suffix:
        variants += call_variants(suffix, this_suffix, "suffix")
        length += len(this_suffix) - len(suffix)

    # We are ready to return the allele name of aliases.
    if alias != marker:
        return "_".join([library["aliases"][alias]["name"]] + variants)

    # Compute CE allele number for the other alleles.
    # TODO: perhaps produce a more intelligent name if there is exactly
    #       1 alias with the same length, or only 1 alias sequence is
    #       contained somewhere within the allele.
    blocknames = []
    blocksize = library.get("block_length", {}).get(marker, 4)
    length -= library.get("length_adjust", {}).get(marker, 0)
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
    known_formats = ("raw", "tssv", "allelename")
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


def reverse_complement(sequence):
    """Return the reverse complement of the given DNA sequence."""
    return "".join(COMPL[x] if x in COMPL else x for x in reversed(sequence))
#reverse_complement


def nnls(A, C, B=None, max_iter=200, min_change=0.0001, debug=False):
    """
    Solve for B in A * B = C in the least squares sense, s.t. B >= 0.

    Hint: call nnls(B.T, C.T).T to solve for A.

    Algorithm has converged if the sum of squared error has decreased
    by less than a factor of min_change in one iteration.  If debug is
    True, print the sum of squared error to sys.stdout after each
    iteration.

    This code was partially adopted from nimfa.methods.factorization.bd.
    """
    import numpy as np
    if B is None:
        B = np.matrix(np.zeros([A.shape[1], C.shape[1]]))
    E = A.T * A
    F = A.T * C
    prev_score = cur_score = sys.float_info.max
    for i in range(max_iter):
        for n in range(B.shape[0]):
            nn = list(range(n)) + list(range(n + 1, B.shape[0]))
            tmp = (F[n, :] - E[n, nn] * B[nn, :]) / E[n, n]
            tmp[np.isnan(tmp)] = 0
            tmp[tmp < 0] = 0
            B[n, :] = tmp
        prev_score = cur_score
        cur_score = np.square(C - A * B).sum()
        score_change = (prev_score-cur_score)/prev_score

        if debug:
            if i:
                print("%4i %15.6f %15.6f %6.2f" % (i, cur_score,
                    prev_score-cur_score, 100*score_change))
            else:
                print("%4i %15.6f" % (i, cur_score))

        if not cur_score or score_change < min_change:
            # We have converged.
            break

    return B
#nnls


def adjust_stats(value, stats=None):
    """
    Adjust the given stats in place with the given observed value and
    return the adjusted stats as well.  If no stats dict is given,
    create a new stats dict with the following initial values:
    {"n": 1, "min": value, "max": value, "mean": value, "m2": 0.0,
     "variance": 0.0}
    """
    value += 0.0
    if not stats:
        return {"n": 1, "min": value, "max": value, "mean": value, "m2": 0.0,
                "variance": 0.0}
    stats["n"] += 1
    delta = value - stats["mean"]
    stats["mean"] += delta / stats["n"]
    stats["m2"] += delta * (value - stats["mean"])
    try:
        stats["variance"] = stats["m2"] / (stats["n"] - 1)
        stats["min"] = min(stats["min"], value)
        stats["max"] = max(stats["max"], value)
    except ZeroDivisionError:
        stats["variance"] = 0
        stats["min"] = value
        stats["max"] = value
    return stats
#adjust_stats


def get_repeat_pattern(seq):
    """Return compiled regular expression that finds repeats of seq."""
    return re.compile("".join(             # For AGAT, one obtains:
        ["(?:" * (len(seq)-1)] +           # (?:(?:(?:
        ["%s)?" % x for x in seq[1:]] +    # G)?A)?T)?
        ["(?:", seq, ")+"] +               # (?AGAT)+
        ["(?:%s" % x for x in seq[:-1]] +  # (?:A(?:G(?:A
        [")?" * (len(seq)-1)]))            # )?)?)?
#get_repeat_pattern


def read_sample_data_file(infile, data, annotation_column=None, seqformat=None,
                          library=None, default_marker=None):
    """Add data from infile to data dict as [marker, allele]=reads."""
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_allele, colid_forward, colid_reverse = \
        get_column_ids(column_names, "allele", "forward", "reverse")

    # Get marker name column if it exists.
    colid_name = get_column_ids(column_names, "name", optional=True)

    # Also try to get annotation column if we have one.
    if annotation_column is not None:
        try:
            colid_annotation = get_column_ids(column_names, annotation_column)
        except:
            annotation_column = None

    found_alleles = []
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        marker = line[colid_name] if colid_name is not None else default_marker
        allele = line[colid_allele] if seqformat is None \
            else ensure_sequence_format(line[colid_allele], seqformat,
                                        library=library, marker=marker)
        if (annotation_column is not None and
                line[colid_annotation].startswith("ALLELE")):
            found_alleles.append((marker, allele))
        data[marker, allele] = map(int,
            (line[colid_forward], line[colid_reverse]))

    return found_alleles
#read_sample_data_file


def reduce_read_counts(data, limit_reads):
    sum_reads = 0
    for markerallele in data:
        sum_reads += sum(data[markerallele])
    if sum_reads <= limit_reads:
        return

    remove = sorted(random.sample(xrange(sum_reads), sum_reads - limit_reads))
    i = 0
    seen = 0
    while i < len(remove) and seen > remove[i]:
        # Skip the reads filtered out above.
        i += 1
    for markerallele in data:
        for direction in (0, 1):
            seen += data[markerallele][direction]
            while i < len(remove) and seen > remove[i]:
                data[markerallele][direction] -= 1
                i += 1
#reduce_read_counts


def get_sample_data(tags_to_files, callback, allelelist=None,
                    annotation_column=None, seqformat=None, library=None,
                    marker=None, homozygotes=False, limit_reads=sys.maxint,
                    drop_samples=0):
    if drop_samples:
        sample_tags = tags_to_files.keys()
        for tag in random.sample(xrange(len(sample_tags)),
                                 int(len(sample_tags) * drop_samples)):
            del tags_to_files[sample_tags[tag]]

    for tag in tags_to_files:
        data = {}
        alleles = set()
        for infile in tags_to_files[tag]:
            infile = sys.stdin if infile == "-" else open(infile, "r")
            alleles.update(read_sample_data_file(infile, data,
                annotation_column, seqformat, library, marker))
            if infile != sys.stdin:
                infile.close()
        if limit_reads < sys.maxint:
            reduce_read_counts(data, limit_reads)
        if allelelist is not None:
            if tag not in allelelist:
                allelelist[tag] = {}
            for markerx, allele in alleles:
                if markerx not in allelelist[tag]:
                    allelelist[tag][markerx] = set()
                allelelist[tag][markerx].add(allele)
            if marker:
                if marker in allelelist[tag]:
                    allelelist[tag] = {marker: allelelist[tag][marker]}
                else:
                    allelelist[tag] = {}
            if homozygotes:
                for markerx in allelelist[tag].keys():
                    if len(allelelist[tag][markerx]) > 1:
                        del allelelist[tag][markerx]
        callback(tag, data)
#get_sample_data


def get_column_ids(column_names, *names, **optional):
    """Find all names in column_names and return their indices."""
    result = []
    for name in names:
        try:
            result.append(column_names.index(name))
        except ValueError:
            if "optional" in optional and optional["optional"]:
                result.append(None)
            else:
                raise ValueError("Column not found in input file: %s" % name)
    if len(result) == 1:
        return result[0]
    return tuple(result)
#get_column_ids


def parse_allelelist(allelelist, convert=None, library=None):
    """Read allele list from open file handle."""
    column_names = allelelist.readline().rstrip("\r\n").split("\t")
    colid_sample, colid_marker, colid_allele = get_column_ids(column_names,
        "sample", "marker", "allele")
    alleles = {}
    for line in allelelist:
        line = line.rstrip("\r\n").split("\t")
        sample = line[colid_sample]
        marker = line[colid_marker]
        allele = line[colid_allele]
        if convert is not None:
            allele = ensure_sequence_format(allele, convert, library=library,
                                            marker=marker)
        if sample not in alleles:
            alleles[sample] = {}
        if marker not in alleles[sample]:
            alleles[sample][marker] = set()
        alleles[sample][marker].add(allele)
    return alleles
#parse_allelelist


def pos_int_arg(value):
    """Convert str to int, raise ArgumentTypeError if not positive."""
    if not value.isdigit() or not int(value):
        raise argparse.ArgumentTypeError(
            "invalid positive int value: '%s'" % value)
    return int(value)
#pos_int_arg


def regex_arg(value):
    """Compile value into a regular expression."""
    try:
        return re.compile(value)
    except re.error as err:
        raise argparse.ArgumentTypeError(err)
#regex_arg


def add_allele_detection_args(parser):
    group = parser.add_argument_group("allele detection options")
    group.add_argument('-a', '--allelelist', metavar="ALLELEFILE",
        type=argparse.FileType('r'),
        help="file containing a list of the true alleles of each sample "
             "(e.g., obtained from allelefinder)")
    group.add_argument('-c', '--annotation-column', metavar="COLNAME",
        help="name of a column in the sample files, which contains a value "
             "beginning with 'ALLELE' for the true alleles of the sample")
#add_allele_detection_args


def add_random_subsampling_args(parser):
    group = parser.add_argument_group("random subsampling options (advanced)")
    group.add_argument('-Q', '--limit-reads', metavar="N", type=pos_int_arg,
        default=sys.maxint,
        help="simulate lower sequencing depth by randomly dropping reads down "
             "to this maximum total number of reads for each sample")
    group.add_argument('-x', '--drop-samples', metavar="N", type=float,
        default=0, help="randomly drop this fraction of input samples")
#add_random_subsampling_args


def add_sequence_format_args(parser, default_format=None, force=False):
    group = parser.add_argument_group("sequence format options")
    if force:
        group.set_defaults(sequence_format=default_format)
    else:
        group.add_argument('-F', '--sequence-format', metavar="FORMAT",
            choices=("raw", "tssv", "allelename"),
            default=default_format,
            help="convert sequences to the specified format: one of "
                 "%(choices)s (default: " + (
                 "no conversion" if default_format is None else default_format)
                 + ")")
    group.add_argument('-l', '--library', metavar="LIBRARY",
        type=parse_library,
        help="library file for sequence format conversion")
#add_sequence_format_args


def add_input_output_args(parser, single_in=False, batch_support=False,
                          report_out=False):
    """Add arguments for opening sample files to the given parser."""
    # Input file options group.
    if not single_in:
        parser.add_argument('infiles', nargs='*', metavar="FILE",
            default=["-"],
            help="the sample data file(s) to process (default: read from "
                 "stdin)")
    elif not batch_support:
        parser.add_argument('infile', nargs='?', metavar="IN",
            default="-",
            help="the sample data file to process (default: read from stdin)")
    else:
        mutex = parser.add_argument_group(
                    "input file options").add_mutually_exclusive_group()
        mutex.add_argument('infile', nargs='?', metavar="IN",
            default="-",
            help="single sample data file to process (default: read from "
                 "stdin)")
        mutex.add_argument("-i", "--input", dest="infiles", nargs="+",
            metavar="IN",
            help="multiple sample data files to process (use with "
                 "-o/--output)")

    # Output file options group.
    group = parser.add_argument_group("output file options")
    if batch_support and single_in:
        mutex = group.add_mutually_exclusive_group()
        mutex.add_argument('outfile', nargs='?', metavar="OUT",
            default=sys.stdout,
            help="the file to write the output to (default: write to stdout)")
        mutex.add_argument('-o', '--output', dest="outfiles", nargs="+",
            metavar="OUT",
            help="list of names of output files to match with input files "
                 "specified with -i/--input, or a format string to construct "
                 "file names from sample tags; e.g., the default value is "
                 "'\\1-%s.out', which expands to 'sampletag-%s.out'" %
                    ((parser.prog.rsplit(" ", 1)[-1],)*2))
    elif batch_support:
        group.add_argument('-o', '--output', dest="outfiles", nargs="+",
            metavar="OUT",
            default=[sys.stdout],
            help="a single file name to write all output to (default: write "
                 "to stdout) OR a list of names of output files to match with "
                 "input files OR a format string to construct file names from "
                 "sample tags; e.g., the value '\\1-%s.out' expands to "
                 "'sampletag-%s.out'" % ((parser.prog.rsplit(" ", 1)[-1],)*2))
    else:
        group.add_argument('-o', '--output', dest="outfile", metavar="FILE",
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="file to write output to (default: write to stdout)")
    if report_out:
        group.add_argument('-R', '--report', metavar="FILE",
            type=argparse.FileType('w'),
            default=sys.stderr,
            help="file to write a report to (default: write to stderr)")

    # Sample tag parsing options group.
    group = parser.add_argument_group("sample tag parsing options",
        "for details about REGEX syntax and capturing groups, check "
        "https://docs.python.org/howto/regex")
    group.add_argument('-e', '--tag-expr', metavar="REGEX", type=regex_arg,
        default=DEF_TAG_EXPR,
        help="regular expression that captures (using one or more capturing "
             "groups) the sample tags from the file names; by default, the "
             "entire file name except for its extension (if any) is captured")
    group.add_argument('-f', '--tag-format', metavar="EXPR",
        default=DEF_TAG_FORMAT,
        help="format of the sample tags produced; a capturing group reference "
             "like '\\n' refers to the n-th capturing group in the regular "
             "expression specified with -e/--tag-expr (the default of '\\1' "
             "simply uses the first capturing group); with a single sample, "
             "you can enter the sample tag here explicitly")
#add_input_output_args


def get_tag(filename, tag_expr, tag_format):
    """Return formatted sample tag from filename using regex."""
    try:
        return tag_expr.search(filename).expand(tag_format)
    except:
        return filename
#get_tag


def get_input_output_files(args, single=False, batch_support=False):
    if single and not batch_support:
        # One infile, one outfile.  Return 2-tuple (infile, outfile).
        if args.infile == "-" and sys.stdin.isatty():
            return False  # No input specified.
        return args.infile, args.outfile


    if not single and not batch_support:
        # N infiles, one outfile.  Return 2-tuple ({tag: infiles}, out).
        if args.infiles == ["-"] and sys.stdin.isatty():
            return False  # No input specified.

        tags_to_files = {}
        for infile in args.infiles:
            tag = get_tag(infile, args.tag_expr, args.tag_format)
            try:
                tags_to_files[tag].append(infile)
            except KeyError:
                tags_to_files[tag] = [infile]
        return tags_to_files, args.outfile


    if single and batch_support:
        # N infiles, N outfiles.  Return generator of (tag, [ins], out).
        # Each yielded tuple should cause a separate run of the tool.
        infiles = args.infiles if "infiles" in args \
                  and args.infiles is not None else [args.infile]
        if infiles == ["-"] and sys.stdin.isatty():
            return False  # No input specified.

        outfiles = args.outfiles if "outfiles" in args \
                   and args.outfiles is not None else [args.outfile]
        if len(outfiles) > 1 and len(outfiles) != len(infiles):
            raise ValueError(
                "Number of input files (%i) is not equal to number of output "
                "files (%i)." % (len(infiles), len(outfiles)))

        tags = [get_tag(infile, args.tag_expr, args.tag_format)
                for infile in infiles]

        if len(outfiles) == 1:
            outfile = outfiles[0]
            if outfile == sys.stdout and len(set(tags)) == 1:
                # Write output of single sample to stdout.
                return ((tag, infiles, outfile) for tag in set(tags))

            # Write output of each sample to its own outfile.
            if outfile == sys.stdout:
                outfile = DEF_OUTFILE_FORMAT
            return ((tag,
                    [infiles[i] for i in range(len(tags)) if tags[i]==tag],
                    open(outfile.replace("\\1", tag).replace("\\2",
                         args.tool), "w")) for tag in set(tags))

        # Link each output file to each input file.
        # Treating files with the same sample tag as separate samples.
        return ((tags[i], [infiles[i]], open(outfiles[i], 'w'))
                for i in range(len(tags)))

    if not single and batch_support:
        # N infiles, one or N outfiles.
        # If one outfile, return ({tag: [infiles]}, outfile).
        # If N outfiles, return generator of (tag, [infiles], outfile).
        raise NotImplementedError(
            "Multi-input with optional multi-output not supported yet.")
#get_input_output_files


def print_db(text, debug):
    """Print text if debug is True."""
    if debug:
        print(text)
#print_db
