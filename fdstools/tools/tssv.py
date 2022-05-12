#!/usr/bin/env python3

#
# Copyright (C) 2022 Jerry Hoogenboom
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
section in the provided library file.  Each time a pair of flanks is
found, that sequence read is linked to the corresponding marker and the
portion of the sequence between the flanks is extracted.  The number of
times each such extracted sequence was encountered is counted, along
with the orientation (strand) in which it was found in the input file.
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
any sequence with an ambiguous base (i.e., not A, C, G, or T) is
removed.  Finally, the -a/--minimum option can be used to filter out
sequences that have been seen only rarely.  When the
-A/--aggregate-filtered option is given, all filtered sequences of each
marker are aggregated and reported as 'Other sequences'.

This tool is an evolution of the original TSSV program [2].

References:
[1] https://github.com/Jerrythafast/FLASH-lowercase-overhang
[2] https://github.com/jfjlaros/tssv
"""
import gzip
import io
import itertools
import math
import os
import sys

from errno import EPIPE
from multiprocessing import Process, SimpleQueue
from threading import Thread, Lock

from ..lib.cli import add_sequence_format_args, add_input_output_args, pos_int_arg,\
                      get_input_output_files
from ..lib.io import try_write_pipe
from ..lib.seq import PAT_SEQ_RAW, reverse_complement, ensure_sequence_format
from ..lib.sg_align import align

__version__ = "2.1.1"


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

IUPAC_BITS = str.maketrans("ACGTUMRSWYKVHDBN", "ABDHHCEFIJLGKMNO")



class TSSV:
    def __init__(self, library, flank_length, threshold, indel_score, dirname, workers,
                 deduplicate, infile):
        # User inputs.
        self.library = library
        self.indel_score = indel_score
        self.workers = workers
        self.deduplicate = deduplicate
        self.lock = Lock()

        # Convert library.
        refseq_store = library.get_structure_store().get_refseq_store()
        self.has_iupac = False
        self.tssv_library = {}
        for marker, reported_range in library.get_ranges().items():
            flanks = list(reported_range.get_option("flanks", (flank_length, flank_length)))
            if not all(flanks):
                raise ValueError("Missing flanking sequence for marker %s" % marker)
            chromosome, start, *_, end = reported_range.location
            if isinstance(flanks[0], int):
                flanks[0] = refseq_store.get_refseq(chromosome, start - flanks[0], start)
            if isinstance(flanks[1], int):
                flanks[1] = refseq_store.get_refseq(chromosome, end + 1, end + 1 + flanks[1])
            flanks[1] = reverse_complement(flanks[1])
            translate_flank = [0 if PAT_SEQ_RAW.match(flank) else 1 for flank in flanks]
            if any(translate_flank):
                self.has_iupac = True
            pair = tuple((flank.translate(IUPAC_BITS) if trans else flank, trans)
                for flank, trans in zip(flanks, translate_flank))
            self.tssv_library[marker] = (
                pair, (
                    math.ceil(len(flanks[0]) * threshold if threshold < 1 else threshold),
                    math.ceil(len(flanks[1]) * threshold if threshold < 1 else threshold)),
                "".join(refseq_store.get_refseq(chromosome, start, end + 1) for start, end in
                    zip(reported_range.location[1::2], reported_range.location[2::2])))

        if not self.tssv_library:
            raise ValueError("No markers were defined in the given library file")

        # Open input file.
        file_format, self.input = init_sequence_file_read(infile)

        # Open output directory if we have one.
        if dirname:
            self.outfiles = prepare_output_dir(dirname, self.tssv_library, file_format)
        else:
            self.outfiles = None

        # Internal state.
        self.sequences = {marker: {} for marker in self.tssv_library}
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
            self.process_results(tuple([record[0], seq]) + record[1:], results)
        if self.deduplicate:
            self.cache[seq] = (True, results)
        else:
            del self.cache[seq]
    #cache_results


    def process_results(self, record, results):
        self.total_reads += 1
        recognised = 0
        for marker, matches, seq1, seq2 in results:
            recognised |= matches
            sequences = self.sequences[marker]
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
            if seq1 is not None:
                counters["fPaired"] += 1
                try:
                    sequences[seq1][0] += 1
                except KeyError:
                    sequences[seq1] = [1, 0]
                if self.outfiles:
                    write_sequence_record(self.outfiles["markers"][marker]["paired"], record)
            elif self.outfiles:
                if matches & 0b0001:
                    write_sequence_record(self.outfiles["markers"][marker]["noend"], record)
                if matches & 0b0010:
                    write_sequence_record(self.outfiles["markers"][marker]["nostart"], record)

            # Search in the reverse strand.
            if seq2 is not None:
                counters["rPaired"] += 1
                try:
                    sequences[seq2][1] += 1
                except KeyError:
                    sequences[seq2] = [0, 1]
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


    def process_file(self):
        if self.workers == 1:
            for seq in self.dedup_reads():
                self.cache_results(seq, process_sequence(
                    self.tssv_library, self.indel_score, self.has_iupac, seq))
        else:
            # Start worker processes.  The work is divided into tasks that
            # require about 1 million alignments each.
            done_queue = SimpleQueue()
            chunksize = 1000000 // (4 * len(self.tssv_library)) or 1
            thread = Thread(target=feeder, args=(self.dedup_reads(), self.tssv_library,
                self.indel_score, self.has_iupac, self.workers, chunksize, done_queue))
            thread.daemon = True
            thread.start()

            # Process the results as they come in.
            # Below is speed-optimized to manage as many workers as possible.
            acquire_lock = self.lock.acquire
            release_lock = self.lock.release
            cache_results = self.cache_results
            for seq, results in itertools.chain.from_iterable(iter(done_queue.get, None)):
                acquire_lock()
                cache_results(seq, results)
                release_lock()
            thread.join()

        # Count number of unique sequences per marker.
        for marker in self.tssv_library:
            self.counters[marker]["unique_seqs"] = len(self.sequences[marker])
    #process_file


    def filter_sequences(self, aggregate_filtered, minimum, missing_marker_action):
        # Aggregate the sequences that we are about to filter out.
        if aggregate_filtered:
            aggregates = {}
            for marker, sequences in self.sequences.items():
                expected_length = self.library.get_range(marker).get_option("expected_allele_length")
                for sequence, (forward, reverse) in sequences.items():
                    if not seq_pass_filt(sequence, forward + reverse, minimum, expected_length):
                        if marker not in aggregates:
                            aggregates[marker] = [0, 0]
                        aggregates[marker][0] += forward
                        aggregates[marker][1] += reverse

        # Filter out sequences with low read counts and invalid bases.
        for marker in tuple(self.sequences):
            expected_length = self.library.get_range(marker).get_option("expected_allele_length")
            self.sequences[marker] = {sequence: counts
                for sequence, counts in self.sequences[marker].items()
                if seq_pass_filt(sequence, sum(counts), minimum, expected_length)}

        # Add aggregate rows if the user requested so.
        if aggregate_filtered:
            for marker in aggregates:
                self.sequences[marker]["Other sequences"] = aggregates[marker]

        # Check presence of all markers.
        if missing_marker_action != "exclude":
            for marker in self.tssv_library:
                if not self.sequences[marker]:
                    if missing_marker_action == "include":
                        self.sequences[marker]["No data"] = [0, 0]
                    else:
                        raise ValueError("Marker %s was not detected!" % marker)
    #filter_sequences


    def write_sequence_tables(self, outfile, seqformat="raw"):
        header = "\t".join(("marker", "sequence", "total", "forward", "reverse")) + "\n"
        outfile.write(header)
        if self.outfiles:
            self.outfiles["sequences"].write(header)
        for marker in sorted(self.sequences):
            if self.outfiles:
                self.outfiles["markers"][marker]["sequences"].write(header)
            sequences = self.sequences[marker]
            for total, seq in sorted(((sum(sequences[s]), s) for s in sequences), reverse=True):
                line = "\t".join(map(str, [
                    marker,
                    ensure_sequence_format(seq, seqformat, from_format="raw", library=self.library,
                        marker=marker),
                    total] + sequences[seq])) + "\n"
                outfile.write(line)
                if self.outfiles:
                    self.outfiles["sequences"].write(line)
                    self.outfiles["markers"][marker]["sequences"].write(line)
    #write_sequence_tables


    def write_statistics_table(self, outfile):
        header = "\t".join(("marker", "unique_seqs", "tPaired", "fPaired", "rPaired",
                            "tLeft", "fLeft", "rLeft", "tRight", "fRight", "rRight")) + "\n"
        if self.outfiles:
            self.outfiles["statistics"].write(header)
        try_write_pipe(outfile, header)
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
            try_write_pipe(outfile, line)

        line = "\ntotal reads\t%i\nunrecognised reads\t%i\n" % (self.total_reads,self.unrecognised)
        if self.outfiles:
            self.outfiles["statistics"].write(line)
        try_write_pipe(outfile, line)
    #write_statistics_table


    def close_output_files(self):
        if self.outfiles:
            for key, value in self.outfiles.items():
                if key == "markers":
                    for markerfiles in value.values():
                        for markerfile in markerfiles.values():
                            markerfile.close()
                else:
                    value.close()
    #close_output_files
#TSSV


def prepare_output_dir(dir, markers, file_format):
    # Create output directories.
    for marker in markers:
        os.makedirs(os.path.join(dir, marker), exist_ok=True)

    # Open output files.
    open_seq_out = gzip.open if file_format.endswith(".gz") else open
    return {
        "sequences": open(os.path.join(dir, "sequences.csv"), "wt", encoding="UTF-8"),
        "statistics": open(os.path.join(dir, "statistics.csv"), "wt", encoding="UTF-8"),
        "unknown": open_seq_out(os.path.join(dir, "unknown" + file_format), "wt", encoding="UTF-8"),
        "markers": {
            marker: {
                "sequences": open(os.path.join(dir, marker, "sequences.csv"), "wt",
                    encoding="UTF-8"),
                "paired": open_seq_out(os.path.join(dir, marker, "paired" + file_format), "wt",
                    encoding="UTF-8"),
                "noend": open_seq_out(os.path.join(dir, marker, "noend" + file_format), "wt",
                    encoding="UTF-8"),
                "nostart": open_seq_out(os.path.join(dir, marker, "nostart" + file_format), "wt",
                    encoding="UTF-8"),
            } for marker in markers
        }
    }
#prepare_output_dir


def align_pair(reference, reference_rc, pair, indel_score):
    left_dist, left_pos = align(
        reference[pair[0][1]], pair[0][0], indel_score, bitwise=pair[0][1])
    right_dist, right_pos = align(
        reference_rc[pair[1][1]], pair[1][0], indel_score, bitwise=pair[1][1])
    return (left_dist, left_pos), (right_dist, len(reference[0]) - right_pos)
#align_pair


def relative_distance(reference, sequence):
    """Return (pct_changed, -len(reference))."""
    len_reference = len(reference)
    if len_reference > len(sequence):
        return (align(reference, sequence, 1, global_align=1)[0] / len_reference, -len_reference)
    return (align(sequence, reference, 1, global_align=1)[0] / len_reference, -len_reference)
#relative_distance


def prune_matched_ranges(tssv_library, matched_ranges):
    """Removes overlapping ranges from the provided list."""
    if len(matched_ranges) < 2:
        return

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
                    group_scores = [relative_distance(tssv_library[marker][2], seq)
                        for start, end, seq, strand, marker in matched_ranges]
                worst_range_index = group_scores.index(max(group_scores[:i - group_start]))
                del matched_ranges[group_start + worst_range_index]
                del group_scores[worst_range_index]
                i = group_start + 1
                continue
            # Group consisted of a single range; new group starts here.
            group_scores = group_scores[1:]
            group_start = i
        i += 1
#prune_matched_ranges


def process_sequence(tssv_library, indel_score, has_iupac, seq):
    """Find markers in sequence."""
    seqs = (seq, reverse_complement(seq))
    seqs_up = tuple((x, x.translate(IUPAC_BITS) if has_iupac else None)
        for x in (seqs[0].upper(), seqs[1].upper()))
    matched_ranges = []
    marker_matches = {}
    for marker, (pair, thresholds, refseq) in tssv_library.items():
        algn = (
            align_pair(seqs_up[0], seqs_up[1], pair, indel_score),
            align_pair(seqs_up[1], seqs_up[0], pair, indel_score))
        matches = 0
        if algn[0][0][0] <= thresholds[0]:
            # Left marker was found in forward sequence
            cutout = seqs[0][max(0, algn[0][0][1] - len(pair[0][0])) : algn[0][0][1]]
            if cutout.lower() != cutout:
                matches += 0b0001
        if algn[0][1][0] <= thresholds[1]:
            # Right marker was found in forward sequence.
            cutout = seqs[0][algn[0][1][1] : algn[0][1][1] + len(pair[1][0])]
            if cutout.lower() != cutout:
                matches += 0b0010
        if algn[1][0][0] <= thresholds[0]:
            # Left marker was found in reverse sequence
            cutout = seqs[1][max(0, algn[1][0][1] - len(pair[0][0])) : algn[1][0][1]]
            if cutout.lower() != cutout:
                matches += 0b0100
        if algn[1][1][0] <= thresholds[1]:
            # Right marker was found in reverse sequence.
            cutout = seqs[1][algn[1][1][1] : algn[1][1][1] + len(pair[1][0])]
            if cutout.lower() != cutout:
                matches += 0b1000
        if (matches & 0b0011) == 0b0011 and algn[0][0][1] < algn[0][1][1]:
            # Matched pair in forward sequence.
            matched_ranges.append((algn[0][0][1], algn[0][1][1],
                seqs[0][algn[0][0][1] : algn[0][1][1]], 0, marker))
        if (matches & 0b1100) == 0b1100 and algn[1][0][1] < algn[1][1][1]:
            # Matched pair in reverse sequence.
            matched_ranges.append((len(seq) - algn[1][1][1], len(seq) - algn[1][0][1],
                seqs[1][algn[1][0][1] : algn[1][1][1]], 1, marker))
        marker_matches[marker] = matches

    # Eliminate overlapping matched pairs.
    prune_matched_ranges(tssv_library, matched_ranges)

    # Compile the results.
    results = []
    for marker, matches in marker_matches.items():
        forward = [match[2] for match in matched_ranges if match[3] == 0 and match[4] == marker]
        reverse = [match[2] for match in matched_ranges if match[3] == 1 and match[4] == marker]
        results.append((marker, matches,
            forward[0] if forward else None,
            reverse[0] if reverse else None))
    return results
#process_sequence


def seq_pass_filt(sequence, reads, threshold, explen=None):
    """Return False if the sequence does not meet the criteria."""
    return (reads >= threshold and PAT_SEQ_RAW.match(sequence) is not None and
        (explen is None or explen[0] <= len(sequence) <= explen[1]))
#seq_pass_filt


def worker(tssv_library, indel_score, has_iupac, task_queue, done_queue):
    """
    Read sequences from task_queue, write findings to done_queue.
    """
    for task in iter(task_queue.get, None):
        done_queue.put(tuple(
            (seq, process_sequence(tssv_library, indel_score, has_iupac, seq)) for seq in task))
#worker


def feeder(input, tssv_library, indel_score, has_iupac, workers, chunksize, done_queue):
    """
    Start worker processes, feed them sequences from input and have them
    write their results to done_queue.
    """
    task_queue = SimpleQueue()
    processes = []
    for i in range(workers):
        process = Process(target=worker,
            args=(tssv_library, indel_score, has_iupac, task_queue, done_queue))
        process.daemon = True
        process.start()
        processes.append(process)
    while True:
        # Sending chunks of reads to the workers.
        task = tuple(itertools.islice(input, chunksize))
        if not task:
            break
        task_queue.put(task)
    for i in range(workers):
        task_queue.put(None)
    for process in processes:
        process.join()
    done_queue.put(None)
#feeder


def write_sequence_record(outfile, record):
    """
    Write a tuple of (header, sequence) to a FastA stream or
    a tuple of (header, sequence, quality) to a FastQ stream.
    """
    outfile.write((">%s\n%s\n" if len(record) == 2 else "@%s\n%s\n+\n%s\n") % record)
#write_sequence_record


def init_sequence_file_read(rawfile):
    """
    Return a 2-tuple with "fasta" or "fastq" and a generator that
    generates tuples of (header, sequence) from a FastA stream or
    tuples of (header, sequence, quality) from a FastQ stream.
    """
    magic = rawfile.peek(2)[:2]
    if magic == b"\x1f\x8b":
        # GZipped input file.
        infile = gzip.open(rawfile, "rt", encoding="UTF-8")
        filetype = ".gz"
    elif magic and magic[0] in b">@":
        # FastA/FastQ input file.
        infile = io.TextIOWrapper(rawfile, encoding="UTF-8")
        filetype = ".fa" if magic[0] == b">" else ".fq"
    else:
        raise ValueError("Input file is not a GZip, FastQ or FastA file")

    firstchar = infile.read(1)
    if not firstchar or firstchar not in ">@":
        raise ValueError("Input GZip file does not contain a FastQ or FastA file")
    if filetype == ".gz":
        filetype = ".fa.gz" if firstchar == ">" else ".fq.gz"

    def genreads():
        state = 0
        for line in infile:
            if state == 1:  # Put most common state on top.
                if firstchar == ">" and line.startswith(">"):
                    yield (header, seq)
                    header = line[1:].strip()
                    seq = ""
                elif firstchar == "@" and line.startswith("+"):
                    qual = ""
                    state = 2
                else:
                    seq += line.strip()
            elif state == 2:
                if line.startswith("@") and len(qual) >= len(seq):
                    yield (header, seq, qual)
                    header = line[1:].strip()
                    seq = ""
                    state = 1
                else:
                    qual += line.strip()
            elif state == 0:
                header = line.strip()
                seq = ""
                state = 1
        yield (header, seq) if state == 1 else (header, seq, qual)

    return filetype, genreads()
#init_sequence_file_read


def run_tssv_lite(infile, outfile, reportfile, library, flank_length, seqformat, threshold,
                  minimum, aggregate_filtered, missing_marker_action, dirname, indel_score,
                  workers, no_deduplicate):
    tssv = TSSV(library, flank_length, threshold, indel_score, dirname, workers,
        not no_deduplicate, infile)
    tssv.process_file()
    tssv.filter_sequences(aggregate_filtered, minimum, missing_marker_action)

    try:
        tssv.write_sequence_tables(outfile, seqformat)
        tssv.write_statistics_table(reportfile)
        tssv.close_output_files()
    except IOError as e:
        if e.errno == EPIPE:
            try:
                try_write_pipe(reportfile, "Stopped early because the output was closed.\n")
            except:
                pass
            return
        raise
#run_tssv_lite


def add_arguments(parser):
    add_sequence_format_args(parser, default_format="raw", force=False, require_library=True)
    add_input_output_args(parser, single_in=True, batch_support=False, report_out=True)
    parser.add_argument("-L", "--flank-length", metavar="N", type=pos_int_arg,
        default=_DEF_FLANK_LENGTH,
        help="length of anchor (flanking) sequences to use, if not specified in the library file "
             "(default: %(default)s)")
    parser.add_argument("-D", "--dir",
        help="output directory for verbose output; when given, a subdirectory "
             "will be created for each marker, each with a separate "
             "sequences.csv file and a number of FASTA/FASTQ files containing "
             "unrecognised reads (unknown.fa), recognised reads "
             "(Marker/paired.fa), and reads that lack one of the flanks of a "
             "marker (Marker/noend.fa and Marker/nostart.fa)")
    parser.add_argument("-T", "--num-threads", metavar="THREADS", type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    parser.add_argument("-X", "--no-deduplicate", action="store_true",
        help="disable deduplication of reads; by setting this option, memory usage will be "
             "reduced in expense of longer running time")
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
    files = get_input_output_files(args, single_in=True, batch_support=False)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output of another program")
    infile = sys.stdin.buffer if files[0] == "-" else open(files[0], "rb")
    run_tssv_lite(infile, files[1], args.report, args.library, args.flank_length,
                  args.sequence_format, args.mismatches, args.minimum,
                  args.aggregate_filtered, args.missing_marker_action,
                  args.dir, args.indel_score, args.num_threads, args.no_deduplicate)
    if infile != sys.stdin.buffer:
        infile.close()
#run
