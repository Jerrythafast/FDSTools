#!/usr/bin/env python

#
# Copyright (C) 2019 Jerry Hoogenboom
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
import sys
import math
import itertools as it
import os

from multiprocessing.queues import SimpleQueue
from multiprocessing import Process
from threading import Thread, Lock
from errno import EPIPE


from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_sequence_format_args, reverse_complement, PAT_SEQ_RAW,\
                  ensure_sequence_format
from ..sg_align import align

__version__ = "2.0.0"


# Default values for parameters are specified below.

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
_DEF_MINIMUM = 1



class TSSV:
    def __init__(self, library, threshold, indel_score, dirname, workers, deduplicate, infile):
        # User inputs.
        self.library = library
        self.indel_score = indel_score
        self.workers = workers
        self.deduplicate = deduplicate
        self.lock = Lock()

        # Convert library.
        self.tssv_library = {marker: (
            (flanks[0], reverse_complement(flanks[1])), (
                int(math.ceil(len(flanks[0]) * threshold)),
                int(math.ceil(len(flanks[1]) * threshold))))
            for marker, flanks in library["flanks"].items()}

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
        self.lock.acquire()
        for record in self.cache[seq][1]:
            self.process_results(tuple([record[0], seq]) + record[1:], results)
        if self.deduplicate:
            self.cache[seq] = (True, results)
        else:
            del self.cache[seq]
        self.lock.release()
    #cache_results


    def process_results(self, record, results):
        self.total_reads += 1
        recognised = 0
        for marker, matches, seq1, seq2 in results:
            recognised |= matches
            self.counters[marker]["fLeft"] += matches & 1
            self.counters[marker]["fRight"] += matches >> 1 & 1
            self.counters[marker]["rLeft"] += matches >> 2 & 1
            self.counters[marker]["rRight"] += matches >> 3 & 1

            # Search in the forward strand.
            if seq1 is not None:
                self.counters[marker]["fPaired"] += 1
                if seq1 not in self.sequences[marker]:
                    self.sequences[marker][seq1] = [1, 0]
                else:
                    self.sequences[marker][seq1][0] += 1
                if self.outfiles:
                    write_sequence_record(self.outfiles["markers"][marker]["paired"], record)
            elif self.outfiles:
                if matches & 1:
                    write_sequence_record(self.outfiles["markers"][marker]["noend"], record)
                if matches & 2:
                    write_sequence_record(self.outfiles["markers"][marker]["nostart"], record)

            # Search in the reverse strand.
            if seq2 is not None:
                self.counters[marker]["rPaired"] += 1
                if seq2 not in self.sequences[marker]:
                    self.sequences[marker][seq2] = [0, 1]
                else:
                    self.sequences[marker][seq2][1] += 1
                if self.outfiles:
                    write_sequence_record(self.outfiles["markers"][marker]["paired"], record)
            elif self.outfiles:
                if matches & 4:
                    write_sequence_record(self.outfiles["markers"][marker]["noend"], record)
                if matches & 8:
                    write_sequence_record(self.outfiles["markers"][marker]["nostart"], record)

        if not recognised:
            self.unrecognised += 1
            if self.outfiles:
                write_sequence_record(self.outfiles["unknown"], record)
    #process_results


    def process_file(self):
        if self.workers == 1:
            for seq in self.dedup_reads():
                self.cache_results(seq, process_sequence(self.tssv_library, self.indel_score, seq))
        else:
            # Start worker processes.  The work is divided into tasks that
            # require about 1 million alignments each.
            done_queue = SimpleQueue()
            chunksize = int(1000000/(4*len(self.tssv_library))) or 1
            thread = Thread(target=feeder, args=(self.dedup_reads(), self.tssv_library,
                self.indel_score, self.workers, chunksize, done_queue))
            thread.daemon = True
            thread.start()
            for seq, results in it.chain.from_iterable(iter(done_queue.get, None)):
                self.cache_results(seq, results)
            thread.join()

        # Count number of unique sequences per marker.
        for marker in self.tssv_library:
            self.counters[marker]["unique_seqs"] = len(self.sequences[marker])
    #process_file


    def filter_sequences(self, aggregate_filtered, minimum, missing_marker_action):
        # Aggregate the sequences that we are about to filter out.
        if aggregate_filtered:
            aggregates = {}
            for marker in self.sequences:
                expected_length = self.library.get("expected_length", {}).get(marker)
                for sequence, (forward, reverse) in self.sequences[marker].items():
                    if not seq_pass_filt(sequence, forward + reverse, minimum, expected_length):
                        if marker not in aggregates:
                            aggregates[marker] = [0, 0]
                        aggregates[marker][0] += forward
                        aggregates[marker][1] += reverse

        # Filter out sequences with low read counts and invalid bases.
        for marker, sequences in self.sequences.items():
            expected_length = self.library.get("expected_length", {}).get(marker)
            self.sequences[marker] = {sequence: counts
                for sequence, counts in sequences.items()
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
                    ensure_sequence_format(seq, seqformat, "raw", self.library, marker),
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
            outfile.write(line)

        line = "\ntotal reads\t%i\nunrecognised reads\t%i\n" % (self.total_reads,self.unrecognised)
        if self.outfiles:
            self.outfiles["statistics"].write(line)
        outfile.write(line)
    #write_statistics_table
#TSSV


def prepare_output_dir(dir, markers, file_format):
    # Create output directories.
    os.makedirs(dir)
    for marker in markers:
        os.mkdir(os.path.join(dir, marker))

    # Open output files.
    return {
        "sequences": open(os.path.join(dir, "sequences.csv"), "w"),
        "statistics": open(os.path.join(dir, "statistics.csv"), "w"),
        "unknown": open(os.path.join(dir, "unknown.f" + file_format[-1]), "w"),
        "markers": {
            marker: {
                "sequences": open(os.path.join(dir, marker, "sequences.csv"), "w"),
                "paired": open(os.path.join(dir, marker, "paired.f" + file_format[-1]), "w"),
                "noend": open(os.path.join(dir, marker, "noend.f" + file_format[-1]), "w"),
                "nostart": open(os.path.join(dir, marker, "nostart.f" + file_format[-1]), "w"),
            } for marker in markers
        }
    }
#prepare_output_dir


def align_pair(reference, reference_rc, pair, indel_score=1):
    left_dist, left_pos = align(reference, pair[0], indel_score)
    right_dist, right_pos = align(reference_rc, pair[1], indel_score)
    return (left_dist, left_pos), (right_dist, len(reference) - right_pos)
#align_pair


def process_sequence(tssv_library, indel_score, seq):
    """Find markers in sequence."""
    seqs = (seq, reverse_complement(seq))
    seqs_up = map(str.upper, seqs)
    results = []
    for marker, (pair, thresholds) in tssv_library.items():
        algn = (
            align_pair(seqs_up[0], seqs_up[1], pair, indel_score),
            align_pair(seqs_up[1], seqs_up[0], pair, indel_score))
        matches = 0
        if algn[0][0][0] <= thresholds[0]:
            # Left marker was found in forward sequence
            cutout = seqs[0][max(0, algn[0][0][1]-len(pair[0])):algn[0][0][1]]
            if cutout.lower() != cutout:
                matches += 1
        if algn[0][1][0] <= thresholds[1]:
            # Right marker was found in forward sequence.
            cutout = seqs[0][algn[0][1][1] : algn[0][1][1]+len(pair[1])]
            if cutout.lower() != cutout:
                matches += 2
        if algn[1][0][0] <= thresholds[0]:
            # Left marker was found in reverse sequence
            cutout = seqs[1][max(0, algn[1][0][1]-len(pair[0])):algn[1][0][1]]
            if cutout.lower() != cutout:
                matches += 4
        if algn[1][1][0] <= thresholds[1]:
            # Right marker was found in reverse sequence.
            cutout = seqs[1][algn[1][1][1] : algn[1][1][1]+len(pair[1])]
            if cutout.lower() != cutout:
                matches += 8
        results.append((marker, matches,
            # Matched pair in forward sequence.
            seqs[0][algn[0][0][1] : algn[0][1][1]] if
                (matches & 3) == 3 and algn[0][0][1] < algn[0][1][1]
                else None,
            # Matched pair in reverse sequence.
            seqs[1][algn[1][0][1] : algn[1][1][1]] if
                (matches & 12) == 12 and algn[1][0][1] < algn[1][1][1]
                else None))
    return results
#process_sequence


def seq_pass_filt(sequence, reads, threshold, explen=None):
    """Return False if the sequence does not meet the criteria."""
    return (reads >= threshold and PAT_SEQ_RAW.match(sequence) is not None and
        (explen is None or explen[0] <= len(sequence) <= explen[1]))
#seq_pass_filt


def worker(tssv_library, indel_score, task_queue, done_queue):
    """
    Read sequences from task_queue, write findings to done_queue.
    """
    for task in iter(task_queue.get, None):
        done_queue.put(tuple((seq, process_sequence(tssv_library, indel_score, seq)) for seq in task))
#worker


def feeder(input, tssv_library, indel_score, workers, chunksize, done_queue):
    """
    Start worker processes, feed them sequences from input and have them
    write their results to done_queue.
    """
    task_queue = SimpleQueue()
    processes = []
    for i in range(workers):
        process = Process(target=worker,
            args=(tssv_library, indel_score, task_queue, done_queue))
        process.daemon = True
        process.start()
        processes.append(process)
    while 1:
        # Sending chunks of reads to the workers.
        task = tuple(it.islice(input, chunksize))
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


def init_sequence_file_read(infile):
    """
    Return a 2-tuple with "fasta" or "fastq" and a generator that
    generates tuples of (header, sequence) from a FastA stream or
    tuples of (header, sequence, quality) from a FastQ stream.
    """
    firstchar = infile.read(1)
    if not firstchar:
        return
    if firstchar not in ">@":
        raise ValueError("Input file is not a FastQ or FastA file")

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

    return "fasta" if firstchar == ">" else "fastq", genreads()
#init_sequence_file_read


def run_tssv_lite(infile, outfile, reportfile, library, seqformat,
                  threshold, minimum, aggregate_filtered,
                  missing_marker_action, dirname, indel_score, workers, no_deduplicate):
    tssv = TSSV(library, threshold, indel_score, dirname, workers, not no_deduplicate, infile)
    tssv.process_file()
    tssv.filter_sequences(aggregate_filtered, minimum, missing_marker_action)

    try:
        tssv.write_sequence_tables(outfile, seqformat)
        tssv.write_statistics_table(reportfile)
    except IOError as e:
        if e.errno != EPIPE:
            raise
#run_tssv_lite


def add_arguments(parser):
    add_sequence_format_args(parser, "raw", False, True)
    add_input_output_args(parser, True, False, True)
    parser.add_argument("-D", "--dir",
        help="output directory for verbose output; when given, a subdirectory "
             "will be created for each marker, each with a separate "
             "sequences.csv file and a number of FASTA/FASTQ files containing "
             "unrecognised reads (unknown.fa), recognised reads "
             "(Marker/paired.fa), and reads that lack one of the flanks of a "
             "marker (Marker/noend.fa and Marker/nostart.fa)")
    parser.add_argument("-T", "--num-threads", metavar="THREADS",
        type=pos_int_arg, default=1,
        help="number of worker threads to use (default: %(default)s)")
    parser.add_argument("-X", "--no-deduplicate", action="store_true",
        help="disable deduplication of reads; by setting this option, memory usage will be "
             "reduced in expense of longer running time")
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument("-m", "--mismatches", type=float,
        default=_DEF_MISMATCHES,
        help="number of mismatches per nucleotide to allow in flanking "
             "sequences (default: %(default)s)")
    filtergroup.add_argument("-n", "--indel-score", metavar="N",
        type=pos_int_arg, default=_DEF_INDEL_SCORE,
        help="insertions and deletions in the flanking sequences are "
             "penalised this number of times more heavily than mismatches "
             "(default: %(default)s)")
    filtergroup.add_argument("-a", "--minimum", metavar="N", type=pos_int_arg,
        default=_DEF_MINIMUM,
        help="report only sequences with this minimum number of reads "
             "(default: %(default)s)")
    filtergroup.add_argument("-A", "--aggregate-filtered", action="store_true",
        help="if specified, sequences that have been filtered (as per the "
             "-a/--minimum option, the expected_allele_length section in the "
             "library file, as well as all sequences with ambiguous bases) "
             "will be aggregated per marker and reported as 'Other sequences'")
    filtergroup.add_argument("-M", "--missing-marker-action", metavar="ACTION",
        choices=("include", "exclude", "halt"),
        default="include",
        help="action to take when no sequences are linked to a marker: one of "
             "%(choices)s (default: %(default)s)")
#add_arguments


def run(args):
    files = get_input_output_files(args, True, False)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    infile = sys.stdin if files[0] == "-" else open(files[0], "r")
    run_tssv_lite(infile, files[1], args.report, args.library,
                  args.sequence_format, args.mismatches, args.minimum,
                  args.aggregate_filtered, args.missing_marker_action,
                  args.dir, args.indel_score, args.num_threads, args.no_deduplicate)
    if infile != sys.stdin:
        infile.close()
#run
