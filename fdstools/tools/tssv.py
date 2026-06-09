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
Link raw reads in a FastA or FastQ file to markers and count the number
of reads for each unique sequence.

Scans a FastA or FastQ file, finding the sequences from the 'flanks'
section in the provided library file (or, if none given, the GRCh38
reference sequences immediately flanking either end of the provided
'genome_position').  Each time a pair of flanks is found, that sequence
read is linked to the corresponding marker and the portion of the
sequence between the flanks is extracted.  The number of times each
such extracted sequence was encountered is counted, along with the
orientation (strand) in which it was found in the input file.
The output is a list of unique sequences found for each marker including
the corresponding counts.

By default, a small number of mismatches is allowed when aligning the
flanks to the reads.  This can be controlled with the -m/--mismatches
option.  Furthermore, when the portion of the sequence to which a flank
aligns is completely written in lowercase letters in the input file,
that match is discarded.  This way, FDSTools works well together with
the paired-end read merging tool FLASH, version 1.2.11/lo, which
(optionally) writes the non-overlapping portion of the reads in
lowercase [1].  Together, this ensures repetitive sequences (such as
STRs) are not truncated when the paired-end reads are merged.

The sequences thus obtained are subsequently filtered in three ways.
First, the 'expected_allele_length' section in the library file may be
used to specify hard limits on the acceptable sequence length for each
marker.  Any unexpectedly short or long sequence is removed.  Second,
any sequence with an ambiguous base (i.e., not A, C, G, or T) that does
not match with any of the fully unambiguous sequences is removed.
Finally, the -a/--minimum option can be used to filter out
sequences that have been seen only rarely.  When the
-A/--aggregate-filtered option is given, all filtered sequences of each
marker are aggregated and reported as 'Other sequences'.

This tool is an evolution of the original TSSV program [2].

References:
[1] https://github.com/Jerrythafast/FLASH-lowercase-overhang
[2] https://github.com/jfjlaros/tssv
"""
import gzip
import re

from bisect import insort
from collections import defaultdict
from contextlib import ExitStack, contextmanager
from functools import partial
from itertools import islice
from math import ceil
from pathlib import Path
from threading import Lock, Semaphore

from ..lib.cli import add_sequence_format_args, pos_int_arg
from ..lib.io import parse_reads, std_or_open, write_sequence_record
from ..lib.seq import PAT_SEQ_RAW, reverse_complement, ensure_sequence_format
from ..lib.sg_align import align
from ..lib.tssv import filter_sequences
from ..lib.util import batched, thread_safe_process_pool

__version__ = "2.3.0"


# Default values for parameters are specified below.

# Default flanking sequence length.
# This value can be overridden by the -L command line option.
_DEF_FLANK_LENGTH = 16

# Default maximum number of mismatches per nucleotide in the flanking
# sequences to allow.
# This value can be overridden by the -m command line option.
_DEF_MISMATCHES = 0.1

# Default penalty multiplier for insertions and deletions in the
# flanking sequences.
# This value can be overridden by the -n command line option.
_DEF_INDEL_SCORE = 2

# Default minimum number of reads to consider.
# This value can be overridden by the -a command line option.
_DEF_MINIMUM = 2

# N bases in the read are treated as mismatches during alignments.
IUPAC_BITS_READ = str.maketrans("ACGTUMRSWYKVHDBN", "ABDHHCEFIJLGKMN@")
IUPAC_BITS_FLANK = str.maketrans("ACGTUMRSWYKVHDBN", "124883569:<7;=>?")
IUPAC_REGEX = {
    "M": "[AC]",
    "R": "[AG]",
    "W": "[ATU]",
    "S": "[CG]",
    "Y": "[CTU]",
    "K": "[GTU]",
    "V": "[ACG]",
    "H": "[ACTU]",
    "D": "[AGTU]",
    "B": "[CGTU]",
    "N": "[ACGTU]"}



class TSSV:
    def __init__(self, reads_generator, library, flank_length, threshold, resolve_ambiguous,
                 output_qualities, distrust_base_quality, indel_score, outdir, workers,
                 deduplicate):
        # User inputs.
        self.library = library
        self.indel_score = indel_score
        self.outfiles = outdir
        self.workers = workers
        self.deduplicate = deduplicate
        self.output_qualities = output_qualities or distrust_base_quality
        self.distrust_base_quality = distrust_base_quality
        self.input = reads_generator

        self.lock = Lock()
        self.semaphore = Semaphore(2 * workers)

        # Convert library.
        refseq_store = library.get_structure_store().get_refseq_store()
        self.tssv_library = {}
        for marker, reported_range in library.get_ranges().items():
            flanks = list(reported_range.get_option("flanks", ["", ""]))
            left_min, left_max, right_min, right_max = reported_range.get_option(
                "flank_length", [flank_length] * 4)
            left_length = sorted((left_min, flank_length, left_max))[1]
            right_length = sorted((right_min, flank_length, right_max))[1]
            autoload_reference = reported_range.get_option("autoload_reference")
            chromosome, start, *_, end = reported_range.location
            if left_length > len(flanks[0]):
                # Extend left flank (to the left) with reference bases.
                flanks[0] = refseq_store.get_refseq(
                    chromosome, start - left_length, start - len(flanks[0]),
                    autoload=autoload_reference) + flanks[0]
            elif left_length < len(flanks[0]):
                # Cut bases from the left end of the left flank.
                flanks[0] = flanks[0][-left_length:]
            if right_length > len(flanks[1]):
                # Extend right flank (to the right) with reference bases.
                flanks[1] += refseq_store.get_refseq(
                    chromosome, end + len(flanks[1]) + 1, end + right_length + 1,
                    autoload=autoload_reference)
            elif right_length < len(flanks[1]):
                # Cut bases from the right end of the right flank.
                flanks[1] = flanks[1][:right_length]
            flanks[1] = reverse_complement(flanks[1])
            self.tssv_library[marker] = (
                [flank.translate(IUPAC_BITS_FLANK) for flank in flanks], (
                    ceil(len(flanks[0]) * threshold if threshold < 1 else threshold),
                    ceil(len(flanks[1]) * threshold if threshold < 1 else threshold)),
                "".join(refseq_store.get_refseq(chromosome, start, end + 1) for start, end in
                    zip(reported_range.location[1::2], reported_range.location[2::2]))
                    .translate(IUPAC_BITS_FLANK))

        if not self.tssv_library:
            raise ValueError("No markers were defined in the given library file")

        # Internal state.
        self.sequences = {marker: defaultdict(lambda: [0, 0, []]) for marker in self.tssv_library}
        if resolve_ambiguous:
            self.ambiguous_sequences = {
                marker: defaultdict(lambda: [0, 0, []]) for marker in self.tssv_library}
        else:
            self.ambiguous_sequences = self.sequences
        self.counters = {marker: {key: 0 for key in
                ("fPaired", "rPaired", "fLeft", "rLeft", "fRight", "rRight")}
            for marker in self.tssv_library}
        self.total_reads = 0
        self.unrecognised = 0
        self.cache = {}
    #__init__


    def dedup_reads(self):
        for record in self.input:
            self.lock.acquire()
            if record[1] in self.cache:
                completed, data = self.cache[record[1]]
                if completed:
                    self.process_results(record, data)
                else:
                    # Store header and (if applicable) quality scores.
                    data.append(record[0::2])
                self.lock.release()
            else:
                self.cache[record[1]] = (False, [record[0::2]])
                self.lock.release()
                yield record[1]
    #dedup_reads


    def cache_results(self, seq, results):
        """NOTE: In a concurrent context, caller must hold self.lock!"""
        for record in self.cache[seq][1]:
            self.process_results((record[0], seq) + record[1:], results)
        if self.deduplicate:
            self.cache[seq] = (True, results)
        else:
            del self.cache[seq]
    #cache_results


    def process_results(self, record, results):
        self.total_reads += 1
        recognised = 0
        for marker, matches, start_end_seq1, start_end_seq2, ambiguous in results:
            recognised |= matches
            sequences = self.ambiguous_sequences[marker] if ambiguous else self.sequences[marker]
            counters = self.counters[marker]
            if matches & 0b0001:
                counters["fLeft"] += 1
            if matches & 0b0010:
                counters["fRight"] += 1
            if matches & 0b0100:
                counters["rLeft"] += 1
            if matches & 0b1000:
                counters["rRight"] += 1

            # Search in the forward strand.
            if start_end_seq1 is not None:
                start, end, seq = start_end_seq1
                counters["fPaired"] += 1
                sequences[seq][0] += 1
                if self.output_qualities:
                    qual = record[2][start : end]
                    update_qualities(sequences[seq][2], qual)
                if self.outfiles:
                    write_sequence_record(self.outfiles["markers"][marker]["paired"], record)
            elif self.outfiles:
                if matches & 0b0001:
                    write_sequence_record(self.outfiles["markers"][marker]["noend"], record)
                if matches & 0b0010:
                    write_sequence_record(self.outfiles["markers"][marker]["nostart"], record)

            # Search in the reverse strand.
            if start_end_seq2 is not None:
                start, end, seq = start_end_seq2
                counters["rPaired"] += 1
                sequences[seq][1] += 1
                if self.output_qualities:
                    qual = record[2][start : end][::-1]
                    update_qualities(sequences[seq][2], qual)
                if self.outfiles:
                    write_sequence_record(self.outfiles["markers"][marker]["paired"], record)
            elif self.outfiles:
                if matches & 0b0100:
                    write_sequence_record(self.outfiles["markers"][marker]["noend"], record)
                if matches & 0b1000:
                    write_sequence_record(self.outfiles["markers"][marker]["nostart"], record)

        if not recognised:
            self.unrecognised += 1
            if self.outfiles:
                write_sequence_record(self.outfiles["unknown"], record)
    #process_results


    def process_ambiguous_sequences(self):
        for marker in self.tssv_library:
            if not self.ambiguous_sequences[marker]:
                continue
            proper_seqs = defaultdict(list)
            for seq, counts in self.sequences[marker].items():
                insort(proper_seqs[len(seq)], (-sum(counts[:2]), seq))
            for ambiguous_seq, ambiguous_counts in self.ambiguous_sequences[marker].items():
                pat_match = re.compile(
                    "".join(IUPAC_REGEX.get(base, base) for base in ambiguous_seq))
                for _, proper_seq in proper_seqs.get(len(ambiguous_seq), []):
                    if pat_match.fullmatch(proper_seq):
                        proper_counts = self.sequences[marker][proper_seq]
                        proper_counts[0] += ambiguous_counts[0]
                        proper_counts[1] += ambiguous_counts[1]
                        element_wise_add(proper_counts[2], ambiguous_counts[2])
                        break
                else:
                    # Failed to resolve ambiguous bases.
                    self.sequences[marker][ambiguous_seq] = ambiguous_counts
            self.ambiguous_sequences[marker].clear()  # Might free some memory.
    #process_ambiguous_sequences


    def flatten_low_quality_bases(self):
        """
        Find bases with cubic mean quality score <= self.distrust_base_quality
        in self.sequences, and try to find another sequence with higher
        read count that matches all higher-quality bases. Add the reads
        of the lower-quality sequence to it if one was found.
        """
        threshold = chr(self.distrust_base_quality + 33)
        for marker, marker_seqs in self.sequences.items():
            grouped_seqs = defaultdict(list)
            for seq, counts in marker_seqs.items():
                insort(grouped_seqs[len(seq)], (-sum(counts[:2]), seq, counts[2]))
            for length, seqs in grouped_seqs.items():
                for i in range(len(seqs) - 1, 0, -1):
                    negative_total, seq, qual = seqs[i]
                    qual = flatten_qualities(qual, -negative_total)
                    if any(q <= threshold for q in qual):
                        pat_match = re.compile(
                            "".join("[ACGTU]" if q <= threshold else b for b, q in zip(seq, qual)))
                        for other in islice(seqs, i):
                            if pat_match.fullmatch(other[1]):
                                ambiguous_counts = marker_seqs[seq]
                                proper_counts = marker_seqs[other[1]]
                                proper_counts[0] += ambiguous_counts[0]
                                proper_counts[1] += ambiguous_counts[1]
                                element_wise_add(proper_counts[2], ambiguous_counts[2])
                                del marker_seqs[seq]
                                break
    #flatten_low_quality_bases


    def process_file(self, resolve_ambiguous):
        if self.workers == 1:
            for seq in self.dedup_reads():
                self.cache_results(seq, process_sequence(
                    self.tssv_library, self.indel_score, seq))
        else:
            with thread_safe_process_pool(self.workers) as pool:
                # The work is divided into batches that
                # require about 1 million alignments each.
                # The semaphore limits the number of batches in memory.
                batch_size = 1000000 // (4 * len(self.tssv_library)) or 1
                for batch in batched(self.dedup_reads(), batch_size):
                    pool.submit(process_batch, self.tssv_library, self.indel_score, batch) \
                        .add_done_callback(self.handle_results)
                    self.semaphore.acquire()

        if resolve_ambiguous:
            self.process_ambiguous_sequences()

        if self.distrust_base_quality:
            self.flatten_low_quality_bases()

        if self.output_qualities:
            # Flatten qualities for quality column.
            for marker, marker_seqs in self.sequences.items():
                for seq_values in marker_seqs.values():
                    # The seq_values are: [forward, reverse, qualities].
                    seq_values[2] = flatten_qualities(seq_values[2], seq_values[0] + seq_values[1])

        # Count number of unique sequences per marker.
        for marker in self.tssv_library:
            self.counters[marker]["unique_seqs"] = len(self.sequences[marker])
    #process_file


    def handle_results(self, future):
        with self.lock:
            self.semaphore.release()  # Allow new work to be submitted.
            cache_results = self.cache_results
            for seq_results in future.result():
                cache_results(*seq_results)
    #handle_results


    def filter_sequences(self, aggregate_filtered, minimum, missing_marker_action):
        for marker, marker_sequences in self.sequences.items():
            is_present = filter_sequences(marker_sequences, self.library.get_range(marker),
                minimum, aggregate_filtered, missing_marker_action == "include")
            if not is_present and missing_marker_action == "halt":
                raise ValueError("Marker %s was not detected!" % marker)
    #filter_sequences


    def write_sequence_tables(self, outfile, seqformat="raw"):
        if self.output_qualities:
            header = "\t".join(
                ("marker", "sequence", "total", "forward", "reverse", "quality")) + "\n"
            seq_columns = 3
        else:
            header = "\t".join(("marker", "sequence", "total", "forward", "reverse")) + "\n"
            seq_columns = 2
        outfile.write(header)
        if self.outfiles:
            self.outfiles["sequences"].write(header)
        for marker in sorted(self.sequences):
            if self.outfiles:
                self.outfiles["markers"][marker]["sequences"].write(header)
            sequences = self.sequences[marker]
            for total, seq in sorted(((sum(sequences[s][:2]), s) for s in sequences), reverse=True):
                if self.output_qualities and len(sequences[seq]) == 2:
                    # Add empty quality column to Other sequences and No data.
                    sequences[seq].append("")
                line = "\t".join(map(str, [
                    marker,
                    ensure_sequence_format(seq, seqformat, from_format="raw", library=self.library,
                        marker=marker),
                    total] + sequences[seq][:seq_columns])) + "\n"
                outfile.write(line)
                if self.outfiles:
                    self.outfiles["sequences"].write(line)
                    self.outfiles["markers"][marker]["sequences"].write(line)
    #write_sequence_tables


    def write_statistics_table(self, outfile):
        if not self.outfiles and not outfile:
            return
        header = "\t".join(("marker", "unique_seqs", "tPaired", "fPaired", "rPaired",
                            "tLeft", "fLeft", "rLeft", "tRight", "fRight", "rRight")) + "\n"
        if self.outfiles:
            self.outfiles["statistics"].write(header)
        if outfile:
            outfile.write(header)
        for marker in sorted(self.counters):
            line = "\t".join(map(str, (
                marker,
                self.counters[marker]["unique_seqs"],
                self.counters[marker]["fPaired"] + self.counters[marker]["rPaired"],
                self.counters[marker]["fPaired"],
                self.counters[marker]["rPaired"],
                self.counters[marker]["fLeft"] + self.counters[marker]["rLeft"],
                self.counters[marker]["fLeft"],
                self.counters[marker]["rLeft"],
                self.counters[marker]["fRight"] + self.counters[marker]["rRight"],
                self.counters[marker]["fRight"],
                self.counters[marker]["rRight"]))) + "\n"
            if self.outfiles:
                self.outfiles["statistics"].write(line)
            if outfile:
                outfile.write(line)

        line = "\ntotal reads\t%i\nunrecognised reads\t%i\n" % (self.total_reads,self.unrecognised)
        if self.outfiles:
            self.outfiles["statistics"].write(line)
        if outfile:
            outfile.write(line)
    #write_statistics_table
#TSSV


def update_qualities(cubed_qualities, quality_string):
    """
    Update the cubed_qualities list by adding cubed values extracted
    from quality_string. The len(cubed_qualities) should be equal to
    len(quality_string), or zero, in which case it gets resized to match.
    """
    if not cubed_qualities:
        # First quality_string, resize cubed_qualities and store values.
        cubed_qualities[:] = ((ord(q) - 33) ** 3 for q in quality_string)
    else:
        element_wise_add(cubed_qualities, ((ord(q) - 33) ** 3 for q in quality_string))
#update_qualities


def element_wise_add(dst_list, src_iterable):
    """Perform dst_list[i] += src_iterable[i] on equally-sized inputs."""
    for i, value in enumerate(src_iterable):
        dst_list[i] += value
#element_wise_add


def flatten_qualities(summed_cubed_qualities, num_reads):
    """Convert the summed_cubed_qualities to an ASCII quality string."""
    return "".join(chr(round((q / num_reads) ** (1/3)) + 33) for q in summed_cubed_qualities)
#flatten_qualities


@contextmanager
def open_outdir(path, markers, file_format, zip_level):
    if path is None:
        yield None
        return

    # Create output directories.
    for marker in markers:
        (path / marker).mkdir(parents=True, exist_ok=True)

    if zip_level:
        # Write .fa.gz or .fq.gz files with the given compression level.
        output_format = file_format.removesuffix(".gz") + ".gz"
        open_seq_out = partial(gzip.open, compresslevel=min(max(1, zip_level), 9))
    else:
        # Write .fa or .fq files without compression.
        output_format = file_format.removesuffix(".gz")
        open_seq_out = open

    with ExitStack() as stack:
        yield {
            "sequences": stack.enter_context(open(path / "sequences.csv", "wt")),
            "statistics": stack.enter_context(open(path / "statistics.csv", "wt")),
            "unknown": stack.enter_context(open_seq_out(path / f"unknown{output_format}", "wt")),
            "markers": {
                marker: {
                    "sequences": stack.enter_context(open(path / marker / "sequences.csv", "wt")),
                    "paired": stack.enter_context(
                        open_seq_out(path / marker / f"paired{output_format}", "wt")),
                    "noend": stack.enter_context(
                        open_seq_out(path / marker / f"noend{output_format}", "wt")),
                    "nostart": stack.enter_context(
                        open_seq_out(path / marker / f"nostart{output_format}", "wt"))
                } for marker in markers
            }
        }
#open_outdir


def align_flanks(reference, reference_rc, flanks, indel_score):
    left_dist, left_pos = align(
        reference, flanks[0], indel_score, bitwise=1)
    right_dist, right_pos = align(
        reference_rc, flanks[1], indel_score, bitwise=1)
    return (left_dist, left_pos), (right_dist, len(reference) - right_pos)
#align_flanks


def relative_distance(reference, sequence):
    """Return 2 * hamming_distance(ref, seq) - len(reference)."""
    len_reference = len(reference)
    if len_reference > len(sequence):
        return 2 * align(reference, sequence, 1, global_align=1, bitwise=1)[0] - len_reference
    return 2 * align(sequence, reference, 1, global_align=1, bitwise=1)[0] - len_reference
#relative_distance


def prune_matched_ranges(tssv_library, matched_ranges):
    """Removes overlapping ranges from the provided list."""
    matched_ranges.sort()
    group_scores = []
    group_start = 0
    i = 1
    while i <= len(matched_ranges):
        if i == len(matched_ranges) or all(
                (matched_ranges[i][1]-matched_ranges[i][0] < y[1]-y[0] and
                    matched_ranges[i][0] + (matched_ranges[i][1]-matched_ranges[i][0])//2 >= y[1])
                or (matched_ranges[i][1]-matched_ranges[i][0] >= y[1]-y[0] and
                    y[1] - (y[1]-y[0])//2 <= matched_ranges[i][0])
                    for y in matched_ranges[group_start:i]):
            # End of group with at least 50% pairwise overlap.
            if i - group_start > 1:
                # Overlapping group; eliminate the worst fit iteratively.
                if not group_scores:
                    group_scores = [
                        relative_distance(tssv_library[marker][2], seq.translate(IUPAC_BITS_READ))
                        for start, end, seq, strand, marker in matched_ranges[group_start : i]]
                worst_range_index = group_scores.index(max(group_scores))
                del matched_ranges[group_start + worst_range_index]
                del group_scores[worst_range_index]
                i = group_start + 1
                continue
            # Group consisted of a single range; new group starts here.
            group_scores = group_scores[1:]
            group_start = i
        i += 1
#prune_matched_ranges


def process_sequence(tssv_library, indel_score, seq):
    """Find markers in sequence."""
    seqs = (seq, reverse_complement(seq))
    seqs_upper = [seq.upper() for seq in seqs]
    seqs_bitwise =  [seq.translate(IUPAC_BITS_READ) for seq in seqs_upper]
    matched_ranges = []
    marker_matches = {}
    for marker, (flanks, thresholds, refseq) in tssv_library.items():
        algn = (
            align_flanks(seqs_bitwise[0], seqs_bitwise[1], flanks, indel_score),
            align_flanks(seqs_bitwise[1], seqs_bitwise[0], flanks, indel_score))
        matches = 0
        if algn[0][0][0] <= thresholds[0]:
            # Left marker was found in forward sequence
            cutout = seqs[0][max(0, algn[0][0][1] - len(flanks[0])) : algn[0][0][1]]
            if cutout.lower() != cutout:
                matches += 0b0001
        if algn[0][1][0] <= thresholds[1]:
            # Right marker was found in forward sequence.
            cutout = seqs[0][algn[0][1][1] : algn[0][1][1] + len(flanks[1])]
            if cutout.lower() != cutout:
                matches += 0b0010
        if algn[1][0][0] <= thresholds[0]:
            # Left marker was found in reverse sequence
            cutout = seqs[1][max(0, algn[1][0][1] - len(flanks[0])) : algn[1][0][1]]
            if cutout.lower() != cutout:
                matches += 0b0100
        if algn[1][1][0] <= thresholds[1]:
            # Right marker was found in reverse sequence.
            cutout = seqs[1][algn[1][1][1] : algn[1][1][1] + len(flanks[1])]
            if cutout.lower() != cutout:
                matches += 0b1000
        if (matches & 0b0011) == 0b0011 and algn[0][0][1] < algn[0][1][1]:
            # Matched pair in forward sequence.
            matched_ranges.append((algn[0][0][1], algn[0][1][1],
                seqs_upper[0][algn[0][0][1] : algn[0][1][1]], 0, marker))
        if (matches & 0b1100) == 0b1100 and algn[1][0][1] < algn[1][1][1]:
            # Matched pair in reverse sequence.
            matched_ranges.append((len(seq) - algn[1][1][1], len(seq) - algn[1][0][1],
                seqs_upper[1][algn[1][0][1] : algn[1][1][1]], 1, marker))
        marker_matches[marker] = matches

    # Eliminate overlapping matched pairs.
    if len(matched_ranges) > 1:
        prune_matched_ranges(tssv_library, matched_ranges)

    # Compile the results.
    results = []
    for marker, matches in marker_matches.items():
        forward = next(((a, b, seq) for a, b, seq, strand, m in matched_ranges
            if strand == 0 and m == marker), None)
        reverse = next(((a, b, seq) for a, b, seq, strand, m in matched_ranges
            if strand == 1 and m == marker), None)
        results.append((marker, matches, forward, reverse, any(
            PAT_SEQ_RAW.match(seq[2]) is None for seq in (forward, reverse) if seq and seq[2])))
    return results
#process_sequence


def process_batch(tssv_library, indel_score, batch):
    return tuple((seq, process_sequence(tssv_library, indel_score, seq)) for seq in batch)
#process_batch


def run_tssv_lite(reads_generator, outfile, reportfile, library, flank_length, seqformat, threshold,
                  minimum, resolve_ambiguous, output_qualities, distrust_base_quality,
                  aggregate_filtered, missing_marker_action, outdir, indel_score, workers,
                  deduplicate):
    tssv = TSSV(reads_generator, library, flank_length, threshold, resolve_ambiguous,
                output_qualities, distrust_base_quality, indel_score, outdir, workers, deduplicate)
    tssv.process_file(resolve_ambiguous)
    tssv.filter_sequences(aggregate_filtered, minimum, missing_marker_action)
    tssv.write_sequence_tables(outfile, seqformat)
    tssv.write_statistics_table(reportfile)
#run_tssv_lite


def add_arguments(parser):
    add_sequence_format_args(parser, default_format="raw", force=False, require_library=True)

    parser.add_argument("infile", nargs="?", metavar="IN", default="-",
        help="the sample data file to process (default: read from stdin)")
    parser.add_argument("outfile", nargs="?", metavar="OUT", default="-",
        help="the file to write the output to (default: write to stdout)")

    outputgroup = parser.add_argument_group("output file options")
    outputgroup.add_argument("-R", "--report", metavar="FILE",
        help="file to write a statistical report to")
    outputgroup.add_argument("-D", "--dir", type=Path,
        help="output directory for verbose output; when given, a subdirectory "
             "will be created for each marker, each with a separate "
             "sequences.csv file and a number of FastA/FastQ files containing "
             "unrecognised reads (unknown.fa), recognised reads "
             "(Marker/paired.fa), and reads that lack one of the flanks of a "
             "marker (Marker/noend.fa and Marker/nostart.fa)")
    outputgroup.add_argument("-z", "--zip-level", metavar="LEVEL", type=pos_int_arg, default=0,
        help="GZip compression level to use for FastA/FastQ output when -D/--dir is specified "
             "(1 is fastest, 9 is best, default: no GZip compression)")

    parser.add_argument("-L", "--flank-length", metavar="N", type=pos_int_arg,
        default=_DEF_FLANK_LENGTH,
        help="length of anchor (flanking) sequences to use (default: %(default)s)")
    parser.add_argument("-T", "--num-threads", metavar="THREADS", type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    parser.add_argument("-X", "--no-deduplicate", dest="deduplicate", action="store_false",
        help="disable deduplication of reads; by setting this option, memory usage will be "
             "reduced in expense of longer running time")
    parser.add_argument("-q", "--output-qualities", action="store_true",
        help="output cubic mean basecall quality scores")
    parser.add_argument("-Q", "--distrust-base-quality", metavar="N", type=pos_int_arg, default=0,
        help="attempt to correct basecalls with cubic mean quality "
             "below or equal to the given value; implies -q")

    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--mismatches", type=float, default=_DEF_MISMATCHES,
        help="number of mismatches (per nucleotide of flanking sequence if less than 1, else "
             "absolute) to allow in flanking sequences, rounded upward (default: %(default)s)")
    filtergroup.add_argument("-n", "--indel-score", metavar="N",
        type=pos_int_arg, default=_DEF_INDEL_SCORE,
        help="insertions and deletions in the flanking sequences are penalised this number of "
             "times more heavily than mismatches (default: %(default)s)")
    filtergroup.add_argument("-a", "--minimum", metavar="N", type=pos_int_arg,
        default=_DEF_MINIMUM,
        help="report only sequences with this minimum number of reads (default: %(default)s)")
    filtergroup.add_argument("-A", "--filter-ambiguous", dest="resolve_ambiguous",
        action="store_false",
        help="don't add sequences containing ambiguous bases to the most common compatible "
             "unambiguous sequence")
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
    with parse_reads(args.infile) as (file_format, reads_generator), \
         std_or_open(args.outfile, "wt") as outfile, \
         std_or_open(args.report, "wt") as reportfile, \
         open_outdir(args.dir, args.library.get_ranges(), file_format, args.zip_level) as outdir:

        if args.output_qualities or args.distrust_base_quality and ".fq" not in file_format:
            raise ValueError(
                f"Basecall quality values are not available in {file_format} input file")

        run_tssv_lite(reads_generator, outfile, reportfile, args.library, args.flank_length,
                      args.sequence_format, args.mismatches, args.minimum, args.resolve_ambiguous,
                      args.output_qualities, args.distrust_base_quality, args.aggregate_filtered,
                      args.missing_marker_action, outdir, args.indel_score, args.num_threads,
                      args.deduplicate)
#run
