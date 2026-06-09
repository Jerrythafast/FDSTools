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

"""
Demultiplex dual-indexed paired-end reads from different samples.

Writes a separate pair of GZipped FastQ files for each sample, based
on a table of expected barcodes. The read barcodes are extracted from
the input read headers. Reads are assigned to the sample whose true
barcode most likely resulted in the read barcode, using a log-scaled
table of base miscall probabilities. If the read barcode is most likely
the result of one of the given artefacts, it is ignored, but the other
barcode from the same pair can still be used if it is unique.

Reads will be classified as 'Undetermined' if the score difference
between the two most likely barcodes is less than the value specified
with the -Y/--filter-diff option. In addition, a maximum score can be
specified with the -Z/--filter-score option.

Scores are calculated on absolute log10 scale, which means that 0 is a
perfect match, and are averaged over the number of called bases in the
barcode read (i.e., ignoring 'N' nucleotides). Therefore, the filter
values supplied on the command line are essentially per-base values.

The built-in 'xleap' scoring matrix contains absolute log10 basecall
probabilities for the XLEAP SBS chemistry, empirically determined from a
series of NextSeq 1000 runs. Using this matrix, the output scores can be
interpreted as absolute log10 likelihoods, and the score difference is a
log10 likelihood ratio.
"""
import gzip
import re

from collections import defaultdict
from contextlib import ExitStack, contextmanager
from functools import lru_cache
from itertools import chain
from pathlib import Path
from queue import Queue
from threading import Thread

from ..lib.cli import pos_int_arg, regex_arg
from ..lib.io import get_column_ids, parse_read_pairs, std_or_open, write_sequence_record
from ..lib.seq import PAT_SEQ_RAW, reverse_complement

__version__ = "1.0.0"


# Default barcode extraction pattern.
# This value can be overridden by the -p command line option.
_DEF_BARCODE_PATTERN = re.compile(r":N:0:([ACGTN]{6,})\+([ACGTN]{6,})")

# Default per-base log10 LR filter value.
# This value can be overridden by the -Y command line option.
_DEF_FILTER_DIFF = 0.5

# Default per-base absolute log10 likelihood filter value.
# This value can be overridden by the -Z command line option.
_DEF_FILTER_SCORE = 100.


# Number of unique read barcode results to keep in cache.
CACHE_SIZE = 1_000_000

# Barcode sequences of Illumina's Indexed PhiX product.
# https://knowledge.illumina.com/instrumentation/general/instrumentation-general-reference_material-list/000009867
PHIX_BARCODE_LENGTH = 10
PHIX_BARCODES = {
    "_PhiX-Index1": ("CGACCTAACG", "GCCGGCGTAT"),
    "_PhiX-Index2": ("TAGTTCGGTA", "TAAGTGACAG"),
    "_PhiX-Index3": ("ACCGGCCGTA", "CAATGGCGCC"),
    "_PhiX-Index4": ("AATCACCAGC", "GTAATTACTG"),
    "_PhiX-Index5": ("AGCGAATTAG", "TCTGATCGCC")}

SCORE_MATRICES = {
    # Traditional scoring matrix that simply counts mismatches.
    # We count mismatches as 5, such that with the default
    # -Y/--filter-diff value, samples can be demultiplexed as
    # long as the second-best match has at least two more
    # mismatches than the best one.
    # When interpreted as absolute log10 probabilities,
    # this corresponds to a 0.03% error rate.
    "mismatch": {
        ("A", "A"): 0.,
        ("C", "A"): 5.,
        ("G", "A"): 5.,
        ("T", "A"): 5.,
        ("A", "C"): 5.,
        ("C", "C"): 0.,
        ("G", "C"): 5.,
        ("T", "C"): 5.,
        ("A", "G"): 5.,
        ("C", "G"): 5.,
        ("G", "G"): 0.,
        ("T", "G"): 5.,
        ("A", "T"): 5.,
        ("C", "T"): 5.,
        ("G", "T"): 5.,
        ("T", "T"): 0.},

    # Absolute log10 basecall probabilities observed in NextSeq data.
    "xleap": {
        ("A", "A"): 0.0009,
        ("A", "C"): 3.2587,
        ("A", "G"): 2.8521,
        ("A", "T"): 3.7211,
        ("C", "A"): 2.0568,
        ("C", "C"): 0.0054,
        ("C", "G"): 2.8872,
        ("C", "T"): 2.6333,
        ("G", "A"): 3.3076,
        ("G", "C"): 3.9644,
        ("G", "G"): 0.0003,
        ("G", "T"): 3.7711,
        ("T", "A"): 3.5858,
        ("T", "C"): 3.1835,
        ("T", "G"): 2.7266,
        ("T", "T"): 0.0012}}



def calculate_score(matrix, barcode, read):
    if len(barcode) != len(read):
        raise ValueError(f"Input barcode {barcode} is different length than read barcode {read}")
    return sum(matrix.get(pair, 0.) for pair in zip(barcode, read))
#calculate_score


def score_all(matrix, barcodes, read):
    return sorted((calculate_score(matrix, barcode, read), barcode) for barcode in barcodes)
#score_all


class Demultiplexer:
    def __init__(self, tags_barcodes, left_artefacts, right_artefacts,
                 matrix, filter_diff, filter_score):
        left_barcodes_to_tags = defaultdict(list)
        right_barcodes_to_tags = defaultdict(list)
        for tag, (left_barcode, right_barcode) in tags_barcodes.items():
            left_barcodes_to_tags[left_barcode].append(tag)
            right_barcodes_to_tags[right_barcode].append(tag)
        self.left_barcodes = tuple(left_barcodes_to_tags.keys() | left_artefacts)
        self.right_barcodes = tuple(right_barcodes_to_tags.keys() | right_artefacts)
        self.tags_barcodes = tags_barcodes
        self.left_artefacts = left_artefacts
        self.right_artefacts = right_artefacts
        self.left_barcodes_to_tags = left_barcodes_to_tags
        self.right_barcodes_to_tags = right_barcodes_to_tags
        self.matrix = matrix
        self.filter_diff = filter_diff
        self.filter_score = filter_score
    #__init__

    @lru_cache(maxsize=CACHE_SIZE)
    def find_likely_tag(self, left_barcode_read, right_barcode_read):
        # Calculate scores for both sides separately.
        left_scores = score_all(self.matrix, self.left_barcodes, left_barcode_read)
        right_scores = score_all(self.matrix, self.right_barcodes, right_barcode_read)
        left_failed = left_scores[0][1] in self.left_artefacts
        right_failed = right_scores[0][1] in self.right_artefacts

        if left_failed and right_failed:
            # Failed to read barcode at both sides.
            score = next_score = scored_bases = 0.0
            matched_tag = None
        elif left_failed:
            # Left side failed, see if right side barcode is unique.
            scored_bases = sum(1 for base in right_barcode_read if base in "ACGT")
            score, right_barcode = right_scores[0]
            matched_tags = self.right_barcodes_to_tags[right_barcode]
            if len(matched_tags) > 1:
                next_score = score
                matched_tag = None
            else:
                next_score = right_scores[1][0]
                matched_tag = matched_tags[0]
        elif right_failed:
            # Right side failed, see if left side barcode is unique.
            scored_bases = sum(1 for base in left_barcode_read if base in "ACGT")
            score, left_barcode = left_scores[0]
            matched_tags = self.left_barcodes_to_tags[left_barcode]
            if len(matched_tags) > 1:
                next_score = score
                matched_tag = None
            else:
                next_score = left_scores[1][0]
                matched_tag = matched_tags[0]
        else:
            # Determine the highest two paired scores.
            scored_bases = sum(1 for b in left_barcode_read + right_barcode_read if b in "ACGT")
            candidate_tags = self.left_barcodes_to_tags[left_scores[0][1]] + \
                self.right_barcodes_to_tags[right_scores[0][1]]
            left_scores = {tag: score
                for score, barcode in left_scores
                for tag in self.left_barcodes_to_tags[barcode]}
            right_scores = {tag: score
                for score, barcode in right_scores
                for tag in self.right_barcodes_to_tags[barcode]}
            (score, matched_tag), (next_score, next_tag) = sorted(
                (left_scores[tag] + right_scores[tag], tag) for tag in self.tags_barcodes)[:2]
            if matched_tag not in candidate_tags:
                # Index hopping has likely occurred.
                matched_tag = None

        # Apply score filter.
        if scored_bases:
            score /= scored_bases
            next_score /= scored_bases
            if score > self.filter_score or next_score - score < self.filter_diff:
                matched_tag = None
        else:
            matched_tag = None

        return matched_tag, score, next_score
    #find_likely_tag
#Demultiplexer


def add_indexed_phix(tags_barcodes, outfiles):
    """Add indexed PhiX to tags_barcodes, send them to Undetermined."""
    left_length, right_length = map(len, next(iter(tags_barcodes.values())))
    if max(left_length, right_length) > PHIX_BARCODE_LENGTH:
        raise ValueError(
            f"The -P/--filter-idx-phix option only supports barcodes "
            f"of up to {PHIX_BARCODE_LENGTH} nucleotides")
    for tag, (left_barcode, right_barcode) in PHIX_BARCODES.items():
        tags_barcodes[tag] = [left_barcode[:left_length], right_barcode[:right_length]]
        outfiles[tag] = outfiles[None]
#add_indexed_phix


def writer_thread(queue):
    while task := queue.get():
        if task is None:
            return
        write_sequence_record(*task)
#writer_thread


def demultiplex(read_generator, outfiles, tags_barcodes, barcode_pattern, left_artefacts,
                right_artefacts, matrix, filter_phix, filter_diff, filter_score, output_scores):
    if filter_phix:
        add_indexed_phix(tags_barcodes, outfiles)

    demultiplexer = Demultiplexer(
        tags_barcodes, left_artefacts, right_artefacts, matrix, filter_diff, filter_score)

    # The gzip module apparently does not support writes from different threads,
    # even if appropriate locking is in place, so we start threads manually
    # to make sure writes to the same file always come from the same thread.
    queue_r1 = Queue(1000)
    queue_r2 = Queue(1000)
    thread_r1 = Thread(target=writer_thread, args=(queue_r1,))
    thread_r2 = Thread(target=writer_thread, args=(queue_r2,))
    thread_r1.start()
    thread_r2.start()

    try:
        for read1, read2 in read_generator:
            barcode_header = barcode_pattern.search(read1[0])
            if not barcode_header:
                raise ValueError(f"Could not extract barcode sequence from header '{read1[0]}'")

            tag, score, next_score = demultiplexer.find_likely_tag(*barcode_header.groups())

            if output_scores:
                read1 = (f"{read1[0]} {score=} {next_score=}",) + read1[1:]
                read2 = (f"{read2[0]} {score=} {next_score=}",) + read2[1:]

            # Write reads to this sample's output files.
            outfile_r1, outfile_r2 = outfiles[tag]
            queue_r1.put((outfile_r1, read1))
            queue_r2.put((outfile_r2, read2))
    finally:
        queue_r1.put(None)
        queue_r2.put(None)
        thread_r1.join()
        thread_r2.join()
#demultiplex


def read_barcodes(instream):
    """Read barcodes from Illumina SampleSheet.csv or plain tab-separated file."""
    # Skip to first non-empty line to determine format.
    for line in instream:
        if line := line.rstrip("\r\n\t ,"): break
    else:
        raise ValueError("Barcodes file is empty")

    if line.startswith("[Header]"):
        # Assume it is a Illumina SampleSheet.csv file.
        # Skip to the [BCLConvert_Data] section.
        for line in instream:
            if line.startswith("[BCLConvert_Data]"): break
        else:
            raise ValueError("Illumina SampleSheet should contain [BCLConvert_Data] section")

    if line.startswith("[BCLConvert_Data]"):
        # [BCLConvert_Data] section in Illumina SampleSheet.csv file.
        tags_barcodes = {}
        headers = instream.readline().rstrip("\r\n\t ,").split(",")
        for required in ("Sample_ID", "Index", "Index2"):
            if required not in headers:
                raise ValueError(f"[BCLConvert_Data] table should contain the {required} column")
        for line in instream:
            values = line.rstrip("\r\n\t ,").split(",")
            if not values[0] or values[0].startswith("["):
                break
            if len(values) != len(headers):
                raise ValueError(f"[BCLConvert_Data] table is malformed ({headers=} {values=})")
            values = dict(zip(headers, values))
            tags_barcodes[values["Sample_ID"]] = [values["Index"], values["Index2"]]

    else:
        # Plain tab-separated file.
        tags_barcodes = {tag: barcodes for tag, *barcodes in (
            line.rstrip("\r\n").split("\t") for line in chain([line], instream))}

    # Validate input before returning.
    if not tags_barcodes:
        raise ValueError("No barcodes specified")
    for tag, barcodes in tags_barcodes.items():
        for barcode in barcodes:
            if not barcode or not PAT_SEQ_RAW.match(barcode):
                raise ValueError(f"Barcode '{barcode}' of sample {tag} is invalid")
        if len(barcodes) != 2:
            raise ValueError(f"Got {len(barcodes)} instead of 2 barcodes for sample {tag}")
    return tags_barcodes
#read_barcodes


def read_score_matrix(instream):
    """Read tab-separated score matrix file."""
    matrix = {}
    column_names = instream.readline().rstrip("\r\n").split("\t")
    first_column = column_names[0]
    if first_column not in ("from", "to"):
        raise ValueError(
            f"Column name of first column in score matrix must be "
            f"'from' or 'to', not '{first_column}'")
    columns = list(zip("ACGT", get_column_ids(column_names[1:], "A", "C", "G", "T")))
    for line in instream:
        base_line, *values = line.rstrip("\r\n").split("\t")
        if base_line not in "ACGT":
            raise ValueError(f"Invalid nucleotide in score matrix: '{base_line}'")
        for base_column, index in columns:
            key = (base_column, base_line) if first_column == "to" else (base_line, base_column)
            matrix[key] = float(values[index])
            if matrix[key] < 0:
                raise ValueError(f"Score matrix cannot include negative value {values[index]}")
    if len(matrix) != 16:
        raise ValueError(f"Each base (A/C/G/T) should appear once in the '{first_column}' column")
    return matrix
#read_score_matrix


@contextmanager
def open_outfiles(tags, path, s0_name, truncate, zip_level):
    path.mkdir(parents=True, exist_ok=True)
    open_kwargs = {
        "mode": "wt" if truncate else "at",
        "compresslevel": min(max(1, zip_level), 9),
        "newline": "\n"}
    with ExitStack() as stack:
        result = {
            tag: (
                stack.enter_context(gzip.open(path / f"{tag}_S{i}_R1_001.fastq.gz", **open_kwargs)),
                stack.enter_context(gzip.open(path / f"{tag}_S{i}_R2_001.fastq.gz", **open_kwargs)))
            for i, tag in enumerate(tags, 1)}
        result[None] = (
            stack.enter_context(gzip.open(path / f"{s0_name}_R1_001.fastq.gz", **open_kwargs)),
            stack.enter_context(gzip.open(path / f"{s0_name}_R2_001.fastq.gz", **open_kwargs)))
        yield result
#open_outfiles


def add_arguments(parser):
    parser.add_argument("barcodefile", metavar="TAGS",
        help="input file with sample tags and barcodes, either tab-separated (without headers) "
             "or Illumina SampleSheet.csv file with [BCLConvert_Data] section (no other features "
             "are supported); both barcodes must be specified in the same orientation as expected "
             "in the headers of FastQ files unless -c/--reverse-right or -X/--xleap is specified")
    parser.add_argument("infile1", metavar="R1", nargs="?", default="-",
        help="input FastQ or FastA file, possibly GZipped, either an interleaved file or R1 "
             "(default: read from stdin)")
    parser.add_argument("infile2", metavar="R2", nargs="?",
        help="input FastQ or FastA file, possibly GZipped, R2")
    parser.add_argument("-p", "--barcode-pattern", metavar="PATTERN",
        type=regex_arg, default=_DEF_BARCODE_PATTERN,
        help=f"regular expression that extracts, from read headers, the two read barcodes in two "
             f"capturing groups (default: '{_DEF_BARCODE_PATTERN.pattern}')")
    parser.add_argument("-c", "--reverse-right", action="store_true",
        help="if specified, the right barcodes from the TAGS file are reverse-complemented "
             "(does not apply to sequences specified through the *-artefacts options below)")

    outputgroup = parser.add_argument_group("output file options")
    outputgroup.add_argument("-D", "--dir", default=Path(), type=Path,
        help="output directory (default: current working directory)")
    outputgroup.add_argument("-t", "--truncate", action="store_true",
        help="if specified, preexisting output files are truncated (default: append output)")
    outputgroup.add_argument("-z", "--zip-level", metavar="LEVEL", type=pos_int_arg, default=6,
        help="GZip compression level (1 is fastest, 9 is best, default: 6)")
    outputgroup.add_argument("-s", "--output-scores", action="store_true",
        help="if specified, append scores to output read headers")

    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--matrix", metavar="MATRIX",
        help="scoring matrix to use: a 5x5 (including headers) tab-separated file with first "
             "column named 'from' or 'to', the 16 scores are nonnegative real values; built-in "
             "files: 'mismatch' (the default) scores matches as 0 and mismatches as 5, 'xleap' "
             "uses empirical absolute log10 basecall probabilities for XLEAP SBS chemistry")
    filtergroup.add_argument("-P", "--filter-idx-phix", action="store_true",
        help="if specified, recognize indexed PhiX barcodes and sort to Undetermined")
    filtergroup.add_argument("-X", "--xleap", action="store_true",
        help="convenience option, equivalent to: -c -m xleap -a GGGGGG (with correct length); "
             "do not use for MiSeq i100 index-first data")
    filtergroup.add_argument("-Y", "--filter-diff", metavar="DIFF",
        type=float, default=_DEF_FILTER_DIFF,
        help="drop reads with score difference below this value (default: %(default)s)")
    filtergroup.add_argument("-Z", "--filter-score", metavar="SCORE",
        type=float, default=_DEF_FILTER_SCORE,
        help="drop reads with score above this value (default: %(default)s)")
    filtergroup.add_argument("-a", "--artefacts", metavar="SEQ", nargs="+", default=[],
        help="sequences to treat as known artefacts for both barcodes")
    filtergroup.add_argument("-l", "--left-artefacts", metavar="SEQ", nargs="+", default=[],
        help="sequences to treat as known artefacts for left barcode")
    filtergroup.add_argument("-r", "--right-artefacts", metavar="SEQ", nargs="+", default=[],
        help="sequences to treat as known artefacts for right barcode")
#add_arguments


def run(args):
    # Read score matrix.
    if args.matrix is None:
        matrix = SCORE_MATRICES["xleap" if args.xleap else "mismatch"]
    elif args.matrix in SCORE_MATRICES:
        matrix = SCORE_MATRICES[args.matrix]
    else:
        with std_or_open(args.matrix) as matrixfile:
            matrix = read_score_matrix(matrixfile)

    # Read barcodes.
    with std_or_open(args.barcodefile) as barcodefile:
        tags_barcodes = read_barcodes(barcodefile)

    # Make it possible to use preexisting Undetermined_S0_*.fastq.gz in
    # the output dir as input, writing our undetermined to S00 instead.
    s0_out = "Undetermined_S0"
    try:
        if (args.dir / "Undetermined_S0_R1_001.fastq.gz").samefile(args.infile1):
            s0_out = "Undetermined_S00"
        elif (args.dir / "Undetermined_S0_R2_001.fastq.gz").samefile(args.infile2):
            s0_out = "Undetermined_S00"
    except:
        pass  # No output dir yet, or reading stdin or interleaved.

    # Deduplicate and validate artefacts.
    artefacts = set(args.artefacts)
    left_artefacts = set(args.left_artefacts) | artefacts
    right_artefacts = set(args.right_artefacts) | artefacts
    for artefact in left_artefacts | right_artefacts:
        if not PAT_SEQ_RAW.match(artefact):
            raise ValueError(f"Invalid artefact sequence '{artefact}'")

    # Recognise "dark read" artefact in XLEAP SBS chemistry
    # and reverse-complement right barcodes.
    if args.xleap:
        left_length, right_length = map(len, next(iter(tags_barcodes.values())))
        left_artefacts.add("G" * left_length)
        right_artefacts.add("G" * right_length)
        args.reverse_right = True

    if args.reverse_right:
        tags_barcodes = {tag: (left, reverse_complement(right))
            for tag, (left, right) in tags_barcodes.items()}
    elif args.filter_idx_phix:
        # TODO: the hardcoded PhiX indexes contain reverse-complemented i5 sequences.
        # If reverse_right is not wanted, we will be ending up adding the wrong sequences.
        raise ValueError("-P/--filter-idx-phix requires -c/--reverse-right")

    with parse_read_pairs(args.infile1, args.infile2, as_fastq=True) as read_generator, \
         open_outfiles(tags_barcodes, args.dir, s0_out, args.truncate, args.zip_level) as outfiles:
        demultiplex(read_generator, outfiles, tags_barcodes, args.barcode_pattern,
                    left_artefacts, right_artefacts, matrix, args.filter_idx_phix,
                    args.filter_diff, args.filter_score, args.output_scores)
#run
