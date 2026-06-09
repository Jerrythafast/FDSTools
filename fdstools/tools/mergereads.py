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
Merge paired-end reads.

When the sequenced DNA fragments are shorter than the combined read
length, which is the case for most forensic STR markers sequenced with
2x300nt or 350+50nt (ForenSeq) reads on an Illumina instrument, the
read pairs will overlap. This tool will extend fragments that don't fit
entirely in one read using its mate. Base call errors in the overlap
region will be corrected when the correct base exists in either read
with a higher quality score. Therefore, running TSSV on the output of
this tool will result in cleaner DNA profiles than using only R1 or R2.

Merging ForenSeq 350+50nt reads is beneficial to fix incorrect base
calls at the end of R1, but the -1/--uncombined-r1 option (without its
optional FILE argument) is required to avoid loss of the longest alleles
in the sample.

Do not specify both the -1/--uncombined-r1 and -2/--uncombined-r2
options without their optional FILE arguments at the same time. This
would cause both reads of an uncombined read pair to be included in the
output, potentially leading to double counting in TSSV.

This tool is an adaptation of the original FLASH program [1].

References:
[1] https://doi.org/10.1093/bioinformatics/btr507
"""
import logging

from contextlib import contextmanager

from ..lib.cli import pos_int_arg
from ..lib.io import parse_read_pairs, std_or_open, write_sequence_record
from ..lib.seq import reverse_complement
from ..lib.pair_align import align_pair

__version__ = "1.0.0"

LOG = logging.getLogger(__name__)


# Default values for parameters are specified below.

# Default minimum required overlap length.
# This value can be overridden by the -m command line option.
_DEF_MIN_OVERLAP = 10

# Default maximum mismatch density.
# This value can be overridden by the -x command line option.
_DEF_MAX_MISMATCH_DENSITY = 0.25

# Sentinel value to indicate uncombined reads go to combined outfile.
COMBINED_OUTFILE = object()

# Lookup table to find the IUPAC code for a pair of nucleotides.
IUPAC_CODE = {
    "AC": "M",
    "AG": "R",
    "AT": "W",
    "AN": "A",
    "CG": "S",
    "CN": "C",
    "CT": "Y",
    "GN": "G",
    "GT": "K",
    "NN": "N",
    "NT": "T"}


def iupac_code(a, b):
    return IUPAC_CODE["".join(sorted((a, b)))]
#iupac_code


def generate_combined_read(read1, read2, overlap_pos, was_outie,
                           lowercase_overhang, ambiguous, qual_eq_threshold):
    read1_seq = read1[1]
    read2_seq = read2[1]
    read1_qual = read1[2]
    read2_qual = read2[2]

    overlap_len = len(read1_seq) - overlap_pos

    combined_seq = [""] * (len(read2_seq) + overlap_pos)
    combined_qual = [0] * (len(read2_seq) + overlap_pos)

    read1_seq_non_overlapping = read1_seq[:overlap_pos]
    if lowercase_overhang:
        combined_seq[:overlap_pos] = read1_seq_non_overlapping.lower()
    else:
        combined_seq[:overlap_pos] = read1_seq_non_overlapping
    combined_qual[:overlap_pos] = read1_qual[:overlap_pos]

    for i in range(overlap_pos, overlap_pos + overlap_len):
        j = i - overlap_pos
        if read1_seq[i] == read2_seq[j]:
            # Same base, output the higher quality value.
            combined_seq[i] = read1_seq[i]
            combined_qual[i] = read1_qual[i] if read1_qual[i] > read2_qual[j] else read2_qual[j]
        else:
            # Different base, calculate new quality value and decide which base to output.
            q1 = ord(read1_qual[i])
            q2 = ord(read2_qual[j])
            combined_qual[i] = chr(max(abs(q1 - q2), 2) + 33)
            if q1 > q2 + qual_eq_threshold:
                # Read1 had higher quality value, output that bese.
                combined_seq[i] = read1_seq[i]
            elif q1 < q2 - qual_eq_threshold:
                # Read2 had higher quality value, output that base.
                combined_seq[i] = read2_seq[j]
            elif ambiguous:
                # Equal quality scores, output ambiguous base.
                combined_seq[i] = iupac_code(read1_seq[i], read2_seq[j])
            elif was_outie:
                # Equal quality scores, output earliest base in respective read (outie pair).
                if (i <= len(read2_seq) - j and read1_seq[i] != "N") or read2_seq[j] == "N":
                    combined_seq[i] = read1_seq[i]
                else:
                    combined_seq[i] = read2_seq[j]
            else:
                # Equal quality scores, output earliest base in respective read (innie pair).
                if (i >= len(read2_seq) - j and read2_seq[j] != "N") or read1_seq[i] == "N":
                    combined_seq[i] = read2_seq[j]
                else:
                    combined_seq[i] = read1_seq[i]

    read2_seq_non_overlap = read2_seq[overlap_len:]
    if lowercase_overhang:
        combined_seq[len(read1_seq):] += read2_seq_non_overlap.lower()
    else:
        combined_seq[len(read1_seq):] += read2_seq_non_overlap
    combined_qual[len(read1_seq):] = read2_qual[overlap_len:]

    header = read2[0] if was_outie else read1[0]
    combined_seq = "".join(combined_seq)
    combined_qual = "".join(combined_qual)

    return header, combined_seq, combined_qual
#generate_combined_read


def merge_reads(read_generator, combined_out, uncombined_r1, uncombined_r2,
                min_overlap, lowercase_overhang, ambiguous, qual_eq_threshold,
                max_mismatch_density, allow_outies):
    outies = 0
    combined_pairs = 0
    uncombined_pairs = 0
    for read1, read2 in read_generator:
        read2_rc = (read2[0], reverse_complement(read2[1]), read2[2][::-1])
        overlap_pos, overlap_len, was_outie = align_pair(
            read1, read2_rc, min_overlap, max_mismatch_density, allow_outies)

        if was_outie:
            read1, read2_rc = read2_rc, read1
            outies += 1

        if overlap_pos is not None:
            combined_pairs += 1
            combined_read = generate_combined_read(
                read1, read2_rc, overlap_pos, was_outie,
                lowercase_overhang, ambiguous, qual_eq_threshold)
            write_sequence_record(combined_out, combined_read)
        else:
            uncombined_pairs += 1
            if uncombined_r1:
                write_sequence_record(uncombined_r1, read1)
            if uncombined_r2:
                write_sequence_record(uncombined_r2, read2)

    total_pairs = combined_pairs + uncombined_pairs
    LOG.info("Read combination statistics:")
    LOG.info(f"Total pairs:      {total_pairs}")
    LOG.info(f"Combined pairs:   {combined_pairs}")
    if allow_outies:
        innies = combined_pairs - outies
        LOG.info(f"    Innie pairs:  {innies} ({(innies / combined_pairs) * 100:.2f}% of combined)")
        LOG.info(f"    Outie pairs:  {outies} ({(outies / combined_pairs) * 100:.2f}% of combined)")
    LOG.info(f"Uncombined pairs: {uncombined_pairs}")
    LOG.info(f"Percent combined: {(combined_pairs / total_pairs) * 100:.2f}%")
#merge_reads


@contextmanager
def open_uncombined(name, combined_outfile, zip_level):
    if name is None:
        yield None
    elif name == COMBINED_OUTFILE:
        yield combined_outfile
    else:
        with std_or_open(name, "wt", zip_level=zip_level) as f:
            yield f
#open_uncombined


def add_arguments(parser):
    parser.add_argument("infile1", metavar="IN", nargs="?", default="-",
        help="input FastQ or FastA file, possibly GZipped, either an interleaved file or R1 "
             "(default: read from stdin)")
    parser.add_argument("infile2", metavar="IN", nargs="?",
        help="input FastQ or FastA file, possibly GZipped, R2")
    parser.add_argument("-m", "--min-overlap", metavar="N", type=pos_int_arg,
        default=_DEF_MIN_OVERLAP,
        help="minimum required overlap length (default: %(default)s)")
    parser.add_argument("-l", "--lowercase-overhang", action="store_true",
        help="output the non-overlapped portion of the merged reads in lowercase, indicating to "
             "TSSV whether the whole target (between flanks) was sequenced on both strands or not")
    parser.add_argument("-a", "--ambiguous", action="store_true",
        help="when mismatching bases have equal quality scores, output the corresponding IUPAC "
             "ambiguous base code (default: output the base that came earliest in its read)")
    parser.add_argument("-E", "--qual-eq-threshold", metavar="N", type=int, default=0,
        help="treat quality scores as equal if their absolute difference is less than or equal to "
             "N (default: %(default)s)")
    outputgroup = parser.add_argument_group("output file options")
    outputgroup.add_argument("-o", "--output", dest="outfile", metavar="FILE", default="-",
        help="file to write output to (default: write to stdout)")
    outputgroup.add_argument("-1", "--uncombined-r1", metavar="FILE", nargs="?",
        const=COMBINED_OUTFILE, default=None,
        help="if specified, output R1 from uncombined pairs (to the given FILE if given, "
             "otherwise together with the combined reads)")
    outputgroup.add_argument("-2", "--uncombined-r2", metavar="FILE", nargs="?",
        const=COMBINED_OUTFILE, default=None,
        help="if specified, output R2 from uncombined pairs (to the given FILE if given, "
             "otherwise together with the combined reads)")
    outputgroup.add_argument("-z", "--zip-level", metavar="LEVEL", type=pos_int_arg, default=0,
        help="GZip compression level to use for FastQ output files "
             "(1 is fastest, 9 is best, default: no GZip compression)")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-x", "--max-mismatch-density", metavar="X", type=float,
        default=_DEF_MAX_MISMATCH_DENSITY,
        help="maximum allowed ratio between the number of mismatched base pairs and the overlap "
             "length, ignoring 'N' bases (default: %(default)s)")
    filtergroup.add_argument("-O", "--allow-outies", action="store_true",
        help="also try combining read pairs in the 'outie' orientation, where the beginning of "
             "both reads overlap instead of their ends, possibly combining reads that still "
             "contain adapter sequence at the end")
#add_arguments


def run(args):
    with parse_read_pairs(args.infile1, args.infile2, as_fastq=True) as read_generator, \
         std_or_open(args.outfile, "wt", zip_level=args.zip_level) as outfile, \
         open_uncombined(args.uncombined_r1, outfile, args.zip_level) as uncombined_r1, \
         open_uncombined(args.uncombined_r2, outfile, args.zip_level) as uncombined_r2:
        merge_reads(read_generator, outfile, uncombined_r1, uncombined_r2,
                    args.min_overlap, args.lowercase_overhang, args.ambiguous,
                    args.qual_eq_threshold, args.max_mismatch_density, args.allow_outies)
#run
