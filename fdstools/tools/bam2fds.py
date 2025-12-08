#!/usr/bin/env python3

#
# Copyright (C) 2025 Netherlands Forensic Institute
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
Convert a BAM file with aligned reads to an FDSTools table.

For single-end data, this tool extracts primary and (optionally)
supplementary alignments.  For paired-end data, only properly paired
(flag 2 set) primary alignments are used.  The output is a list of
unique sequences found for each marker including the corresponding
counts.

The sequences thus obtained are subsequently filtered in three ways.
First, the 'expected_allele_length' section in the library file may be
used to specify hard limits on the acceptable sequence length for each
marker.  Any unexpectedly short or long sequence is removed.  Second,
any sequence with an ambiguous base (i.e., not A, C, G, or T) is
removed.  Finally, the -a/--minimum option can be used to filter out
sequences that have been seen only rarely.  When the
-A/--aggregate-filtered option is given, all filtered sequences of each
marker are aggregated and reported as 'Other sequences'.
"""
import argparse
import sys

from concurrent.futures import ProcessPoolExecutor

from ..lib.cli import add_sequence_format_args, pos_int_arg
from ..lib.seq import ensure_sequence_format
from ..lib.tssv import filter_sequences

try:
    from itertools import batched  # Python 3.12+
except ImportError:
    from itertools import islice
    def batched(iterable, n):
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch


__version__ = "1.0.0"


# Default minimum number of reads to consider.
# This value can be overridden by the -a command line option.
_DEF_MINIMUM = 2



def clip_segment_to_range(segment, start, end, with_qual=False, with_cigar=False):
    """
    Clip an aligned segment's sequence to the given coordinates and
    return it.  If with_qual and/or with_cigar is True, return a
    two-tuple or three-tuple of the requested values in the order
    (clipped_sequence, clipped_qualities, hardclipped_cigartuples).
    """
    return clip_to_range(start, end,
        segment.query_alignment_sequence,
        segment.reference_start,
        segment.reference_end,
        segment.cigartuples,
        segment.query_alignment_qualities if with_qual else None,
        with_cigar)
#clip_segment_to_range


def clip_to_range(start, end, seq, orig_start, orig_end, cigar, qual=None, output_cigar=False):
    """
    Clip a sequence that aligned to [orig_start:orig_end] to the given
    coordinates and return it.  If qual is given and/or output_cigar
    is True, return a tuple of the requested values in the order
    (clipped_sequence, clipped_qualities, hardclipped_cigartuples).
    """
    # There is a table with cigar opcodes and their meaning on page 8 of https://samtools.github.io/hts-specs/SAMv1.pdf.
    cigar = [(5 if op == 4 else op, length) for op, length in cigar]
    ref_start = orig_start
    ref_end = orig_end
    if ref_start < start:
        for i, (op, length) in enumerate(cigar):
            if op in {0, 1, 2, 3, 7, 8}:
                if op != 1:
                    # Aligned, deleted or skipped (e.g. intron in mRNA)
                    # bases: advance ref_start.
                    if length >= start - ref_start:
                        cigar[i] = (op, length - (start - ref_start))
                        length = start - ref_start
                    ref_start += length
                if op in {0, 1, 7, 8}:
                    # Aligned or inserted bases: clip them off.
                    seq = seq[length:]
                    if qual:
                        qual = qual[length:]
            if ref_start == start:
                if cigar[i][1] == 0:
                    i += 1
                break
        preexisting_clip = cigar[0][1] if cigar[0][0] == 5 else 0
        resulting_clip = start - orig_start + preexisting_clip
        cigar = [(5, resulting_clip)] + cigar[i:] if resulting_clip else cigar[i:]
    if ref_end > end:
        for i, (op, length) in enumerate(reversed(cigar), 1):
            if op in {0, 1, 2, 3, 7, 8}:
                if op != 1:
                    # Aligned, deleted or skipped (e.g. intron in mRNA)
                    # bases: advance ref_end.
                    if length >= ref_end - end:
                        cigar[-i] = (op, length - (ref_end - end))
                        length = ref_end - end
                    ref_end -= length
                if op in {0, 1, 7, 8}:
                    # Aligned or inserted bases: clip them off.
                    seq = seq[:-length]
                    if qual:
                        qual = qual[:-length]
            if ref_end == end:
                if cigar[-i][1] == 0:
                    i += 1
                break
        preexisting_clip = cigar[-1][1] if cigar[-1][0] == 5 else 0
        resulting_clip = orig_end - end + preexisting_clip
        cigar = cigar[:len(cigar) - i + 1]
        if resulting_clip:
            cigar += [(5, resulting_clip)]
    if qual is None:
        if output_cigar:
            return seq, cigar
        return seq
    if output_cigar:
        return seq, qual, cigar
    return seq, qual
#clip_to_range


def overlap_consensus(seq1, seq2, qual1, qual2, cigar1, cigar2):
    """
    Calculate a consensus from two sequences that aligned to the same
    reference location, return the two-tuple (sequence, cigar).
    """
    cigar1 = [(op, length) for op, length in cigar1 if op not in {4, 5}]
    cigar2 = [(op, length) for op, length in cigar2 if op not in {4, 5}]
    if cigar1 != cigar2:
        # TODO: Handle this better; perhaps also use MAPQ?
        return seq1, cigar1
        raise ValueError(f"Different cigars for overlap: {cigar1} and {cigar2}")
    if len(set(map(len, (seq1, seq2, qual1, qual2)))) != 1:
        raise ValueError(
            f"Inconsistent lengths for overlap: seq1={len(seq1)}, seq2={len(seq2)}, "
            f"qual1={len(qual1)}, qual2={len(qual2)}")
    if seq1 == seq2:
        return seq1, cigar1
    seq = "".join(b1 if q1 >= q2 else b2 for b1, b2, q1, q2 in zip(seq1, seq2, qual1, qual2))
    return seq, cigar1
#overlap_consensus


def process_marker_reads(bam, reported_range, include_supplementary):
    """
    Yield tuples of (sequence, is_reverse) for reads that aligned to
    the given reported range.
    """
    chr, start, *_, end = reported_range.location
    start -= 1
    pair_cache = {}
    for segment in bam.fetch("chr" + chr, start, end):
        if segment.is_paired:
            yield from process_paired_read(segment, start, end, pair_cache)
        else:
            yield from process_unpaired_read(segment, start, end, include_supplementary)
    # TODO: If we want to do anything with incomplete pairs,
    #       look into pair_cache here.
#process_marker_reads


def process_unpaired_read(segment, start, end, include_supplementary):
    """
    Yield a tuple of (sequence, is_reverse) if the aligned segment is
    primary (or supplementary if allowed) and spans the given range.
    """
    if segment.is_unmapped:
        return  # TODO: could send to TSSV!
    if segment.is_supplementary and not include_supplementary:
        return
    if segment.is_secondary:
        return  # Cannot use this as it does not contain the sequence.
    if segment.reference_start <= start and segment.reference_end >= end:
        yield clip_segment_to_range(segment, start, end), segment.is_reverse
#process_unpaired_read


def process_paired_read(segment, start, end, pair_cache):
    """
    If segment is not a primary alignment, return immediately.
    If segment.query_name is not in pair_cache, add it and return.
    Else, pop segment.query_name from the pair_cache and yield a tuple
    of (sequence, is_reverse) if the aligned segment forms a proper
    pair with its mate, and overlaps the given range completely.
    """
    if segment.is_unmapped:
        return  # TODO: could send to TSSV!
    if segment.is_supplementary:
        return  # TODO: This requires some further studying.
    if segment.is_secondary:
        return  # Cannot use this as it does not contain the sequence.
    if not segment.is_proper_pair:
        return  # TODO: Mate mapped far away, as 'outie', or not at all.
    if segment.query_name not in pair_cache:
        pair_cache[segment.query_name] = segment
        return  # Wait for its mate.

    mate = pair_cache.pop(segment.query_name)
    if segment.is_read2:
        segment, mate = mate, segment
    ref_start1 = segment.reference_start
    ref_end1 = segment.reference_end
    ref_start2 = mate.reference_start
    ref_end2 = mate.reference_end
    ref_start_min = min(ref_start1, ref_start2)
    ref_start_max = max(ref_start1, ref_start2)
    ref_end_min = min(ref_end1, ref_end2)
    ref_end_max = max(ref_end1, ref_end2)
    if ref_start_max > ref_end_min:
        # The read pair doesn't overlap (and doesn't touch).
        return  # TODO: could treat them as separate reads.
    if ref_start_min > start or ref_end_max < end:
        return  # Does not fully span marker.

    # Make a consensus of the read pair.
    seq, cigar = clip_segment_to_range(
        segment if ref_start1 < ref_start2 else mate, ref_start_min, ref_start_max,
        with_cigar=True)
    overlap1, o1q, o1c = clip_segment_to_range(
        segment, ref_start_max, ref_end_min, with_qual=True, with_cigar=True)
    overlap2, o2q, o2c = clip_segment_to_range(
        mate, ref_start_max, ref_end_min, with_qual=True, with_cigar=True)
    overlap_seq, overlap_cigar = overlap_consensus(overlap1, overlap2, o1q, o2q, o1c, o2c)
    seq += overlap_seq
    cigar += overlap_cigar
    seq2, cigar2 = clip_segment_to_range(
        segment if ref_end1 > ref_end2 else mate, ref_end_min, ref_end_max, with_cigar=True)
    seq += seq2
    cigar += cigar2
    cigar = [c for i, c in enumerate(cigar) if c[0] != 5 or i in {0, len(cigar) - 1}]
    yield clip_to_range(start, end, seq, ref_start_min, ref_end_max, cigar), segment.is_reverse
#process_paired_read


def process_batch(infile, batch, include_supplementary):
    with AlignmentFile(infile) as alignment_file:
        results = [None] * len(batch)
        for i, (marker, reported_range) in enumerate(batch):
            sequences = {}
            for seq, is_reverse in process_marker_reads(
                    alignment_file, reported_range, include_supplementary):
                counts = sequences.get(seq, [0, 0])
                counts[is_reverse] += 1
                sequences[seq] = counts
            results[i] = (marker, sequences)
    return results
#process_batch


def print_marker_results(outfile, library, seqformat, minimum, aggregate_filtered, missing_marker_action, result):
    for marker, sequences in result:
        is_present = filter_sequences(sequences, library.get_range(marker), minimum,
                                      aggregate_filtered, missing_marker_action == "include")
        if not is_present and missing_marker_action == "halt":
            raise ValueError("Marker %s was not detected!" % marker)
        for total, seq in sorted(((sum(sequences[s]), s) for s in sequences), reverse=True):
            print("\t".join(map(str, [
                marker,
                ensure_sequence_format(seq, seqformat, from_format="raw", library=library,
                    marker=marker),
                total] + sequences[seq])), file=outfile)
#print_marker_results


def read_bam(infile, outfile, library, seqformat, include_supplementary, minimum,
             aggregate_filtered, missing_marker_action, num_threads):
    print("\t".join(("marker", "sequence", "total", "forward", "reverse")), file=outfile)
    with ProcessPoolExecutor(num_threads) as pool:
        batch_size = max(1, len(library.get_ranges()) // num_threads)
        for batch in batched(library.get_ranges().items(), batch_size):
            pool.submit(process_batch, infile, batch, include_supplementary).add_done_callback(
                lambda future: print_marker_results(outfile, library, seqformat, minimum,
                aggregate_filtered, missing_marker_action, future.result()))
#read_bam


def add_arguments(parser):
    add_sequence_format_args(parser, default_format="raw", force=False, require_library=True)
    parser.add_argument("-T", "--num-threads", metavar="THREADS", type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    parser.add_argument("infile", metavar="IN",
        help="the sample data file to process (must be a .bai indexed bam file)")
    parser.add_argument("outfile", nargs="?", metavar="OUT",
        type=argparse.FileType("tw", encoding="UTF-8"),
        default=sys.stdout, help="the file to write the output to (default: write to stdout)")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-S", "--no-supplementary", action="store_true",
        help="ignore supplementary alignments")
    filtergroup.add_argument("-a", "--minimum", metavar="N", type=pos_int_arg,
        default=_DEF_MINIMUM,
        help="report only sequences with this minimum number of reads (default: %(default)s)")
    filtergroup.add_argument("-B", "--no-aggregate-filtered", dest="aggregate_filtered",
        action="store_false",
        help="by default, sequences that have been filtered (as per the -a/--minimum option, the "
             "expected_allele_length section in the library file, as well as all sequences with "
             "ambiguous bases) will be aggregated per marker and reported as 'Other sequences'; "
             "specify this option to remove such sequences entirely")
    filtergroup.add_argument("-M", "--missing-marker-action", metavar="ACTION",
        choices=("include", "exclude", "halt"), default="include",
        help="action to take when no sequences are linked to a marker: one of "
             "%(choices)s (default: %(default)s)")
#add_arguments


def run(args):
    global AlignmentFile
    try:
        from pysam import AlignmentFile
    except ImportError:
        try:
            from bamnostic import AlignmentFile

            # Patch Bamnostic to support the query_alignment_qualities property.
            import array
            from bamnostic.core import AlignedSegment, cigar_alignment
            AlignedSegment.query_alignment_qualities = property(lambda self: array.array('B',
                (x[0] for x in cigar_alignment(self.query_qualities, self.cigar, query=True))))
        except ImportError:
            raise ValueError(
                "please install pysam (unix systems) or bamnostic to use this tool, "
                "by running 'pip install pysam' or 'pip install bamnostic'")
    read_bam(args.infile, args.outfile, args.library, args.sequence_format,
             not args.no_supplementary, args.minimum, args.aggregate_filtered,
             args.missing_marker_action, args.num_threads)
#run
