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
import gzip
import sys

from contextlib import contextmanager, nullcontext
from functools import partial
from io import TextIOWrapper

from .seq import SEQ_SPECIAL_VALUES, ensure_sequence_format


@contextmanager
def std_or_open(name, mode="rt", *, zip_level=0):
    """
    Open the named file with the given mode.
    Use sys.stdin or sys.stdout if name is "-".

    Files are opened in text mode with UTF-8 encoding
    unless mode includes "b".

    If zip_level is nonzero, the file is opened through gzip.open()
    with the given compression level. The compression level only
    affects GZip output.
    """
    if name is None:
        yield None
        return

    if name == "-":
        # Read from sys.stdin, write to sys.stdout.
        name = sys.stdout if any(m in mode for m in "wxa") else sys.stdin

        if name is sys.stdin and sys.stdin.isatty():
            # Avoid attempting to read input from an interactive shell.
            raise ValueError(
                "Please specify an input file, or redirect output from another program")

        # Switch to binary I/O as needed.
        name = name.buffer if "b" in mode or zip_level else name

        # If we don't need a GZip wrapper, we're done here.
        # Else, we fall through to open a GZip wrapper below.
        if not zip_level:
            yield name
            return

    # Open the output file with GZip.open() if zip_level is nonzero.
    encoding = None if "b" in mode else "utf-8"
    open_func = partial(
        gzip.open, compresslevel=min(max(1, zip_level), 9)) if zip_level else open
    with open_func(name, mode, encoding=encoding) as f:
        yield f
#std_or_open


def get_column_ids(column_names, *names, optional=False):
    """Find all names in column_names and return their indices."""
    result = []
    for name in names:
        try:
            result.append(column_names.index(name))
        except ValueError:
            if optional:
                result.append(None)
            else:
                raise ValueError("Column not found in input file: %s" % name)
    if len(result) == 1:
        return result[0]
    return tuple(result)
#get_column_ids


def parse_flags(flags):
    """Convert comma-separated string of flags to list."""
    return [flag for flag in map(str.strip, flags.split(",")) if flag]


def parse_allelelist(allelefile, *, convert=None, library=None):
    """Read allele list from allelefile."""
    with std_or_open(allelefile) as allelelist:
        column_names = allelelist.readline().rstrip("\r\n").split("\t")
        if column_names == [""]:
            return {}  # Empty file.
        colid_sample, colid_marker, colid_allele = get_column_ids(column_names,
            "sample", "marker", "allele")
        alleles = {}
        for line in allelelist:
            line = line.rstrip("\r\n").split("\t")
            sample = line[colid_sample]
            marker = line[colid_marker]
            allele = line[colid_allele]
            if convert is not None:
                allele = ensure_sequence_format(allele, convert, library=library, marker=marker)
            alleles.setdefault(sample, {}).setdefault(marker, set()).add(allele)
    return alleles
#parse_allelelist


def prune_allelelist(allelelist, retained_tags, *, homozygotes_only=False, only_marker=None):
    """Delete unused and unwanted entries from the given allelelist."""
    for tag in tuple(allelelist):
        if tag not in retained_tags:
            del allelelist[tag]
            continue
        allelelist[tag] = {
            marker: alleles
            for marker, alleles in allelelist[tag].items()
            if (only_marker is None or marker == only_marker)
            and (not homozygotes_only or len(alleles) == 1)}
#prune_allelelist


def read_sample_data_file(infile, seqformat=None, library=None, *, drop_special_seq=False,
                          after_correction=False, combine_strands=False, extra_columns=None,
                          only_marker=None, alleles_out=None):
    """
    Read sample data from infile and return {(marker, sequence): [total]}
    if combine_strands is True, else {(marker, sequence): [forward, reverse]}.

    If seqformat is specified, convert sequences to this format.

    If drop_special_seq is True, ignore the SEQ_SPECIAL_VALUES.

    If after_correction is True, use '*_corrected' columns if present.

    The extra_columns can be a mapping of {column_name: is_optional},
    if not None append a dict to the output lists with this data.

    If only_marker is not None, retain only the specified marker.

    If alleles_out is given, require the flags column to be present and
    call alleles_out.setdefault(marker, set()).add(sequence) for all
    sequences that have the "allele" flag set.
    """
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return {}  # Empty file.
    colid_marker, colid_sequence = get_column_ids(column_names, "marker", "sequence")
    if combine_strands:
        colid_total = None
        numtype = int
        if after_correction:
            colid_total = get_column_ids(column_names, "total_corrected", optional=True)
        if colid_total is None:
            colid_total = get_column_ids(column_names, "total")
        else:
            numtype = float
    else:
        colid_forward = None
        colid_reverse = None
        numtype = int
        if after_correction:
            colid_forward, colid_reverse = get_column_ids(column_names,
                "forward_corrected", "reverse_corrected", optional=True)
        if colid_forward is None:
            colid_forward = get_column_ids(column_names, "forward")
        else:
            numtype = float
        if colid_reverse is None:
            colid_reverse = get_column_ids(column_names, "reverse")
        else:
            numtype = float
    if extra_columns is not None:
        extra_colids = {name: i for name, i in
            ((name, get_column_ids(column_names, name, optional=optional))
                for name, optional in extra_columns.items())
            if i is not None}

    # Require flags column if we are supposed to fill alleles_out.
    colid_flags = None if alleles_out is None else get_column_ids(column_names, "flags")

    data = {}
    for line in infile:
        if line.startswith("info_") or line.isspace():
            continue

        line = line.rstrip("\r\n").split("\t")
        if drop_special_seq and line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        marker = line[colid_marker]
        if only_marker is not None and marker != only_marker:
            continue
        sequence = line[colid_sequence] if seqformat is None else ensure_sequence_format(
            line[colid_sequence], seqformat, library=library, marker=marker)
        if colid_flags is not None and "allele" in parse_flags(line[colid_flags]):
            alleles_out.setdefault(marker, set()).add(sequence)
        if combine_strands:
            data[marker, sequence] = [numtype(line[colid_total])]
        else:
            data[marker, sequence] = list(map(numtype, (line[colid_forward], line[colid_reverse])))
        if extra_columns is not None:
            data[marker, sequence].append({name: line[i] for name, i in extra_colids.items()})
    return data
#read_sample_data_file


def read_fasta(infile):
    """Generate tuples of (header, sequence) from FastA input."""
    header = infile.readline().lstrip(">").strip()
    seq = infile.readline().strip()
    for line in infile:
        if line.startswith(">"):
            yield (header, seq)
            header = line[1:].strip()
            seq = infile.readline().strip()
        else:
            seq += line.strip()
    yield (header, seq)
#read_fasta


def read_fasta_as_fastq(infile):
    """Generate tuples of (header, seq, "#" * len(seq)) from FastA input."""
    header = infile.readline().lstrip(">").strip()
    seq = infile.readline().strip()
    for line in infile:
        if line.startswith(">"):
            yield (header, seq, "#" * len(seq))
            header = line[1:].strip()
            seq = infile.readline().strip()
        else:
            seq += line.strip()
    yield (header, seq, "#" * len(seq))
#read_fasta_as_fastq


def read_fastq(infile):
    """Generate tuples of (header, sequence, quality) from FastQ input."""
    header = infile.readline().lstrip("@").strip()
    seq = infile.readline().strip()
    reading_seq = True
    for line in infile:
        if reading_seq:
            if line.startswith("+"):
                reading_seq = False
                qual = infile.readline().strip()
            else:
                seq += line.strip()
        else:
            if len(qual) >= len(seq) or line.startswith("@"):
                yield (header, seq, qual)
                header = line[1:].strip()
                seq = infile.readline().strip()
                reading_seq = True
            else:
                qual += line.strip()
    yield (header, seq, qual)
#read_fastq


@contextmanager
def parse_reads(filename, as_fastq=False):
    """
    Context manager opening the given file for reading, then yielding
    a 2-tuple (input_format, reads_generator).

    The input_format is either ".fa", ".fa.gz", ".fq" or ".fq.gz".
    The generator produces tuples of (header, sequence) from a FastA
    stream if as_fastq is False or tuples of (header, sequence, quality)
    otherwise (with FastQ input or if as_fastq is True).
    If filename is "-", read from sys.stdin.
    """
    with std_or_open(filename, "rb") as rawfile:
        magic = rawfile.peek(2)[:2]
        if magic == b"\x1f\x8b":
            # GZipped input file.
            wrapper = gzip.open(rawfile, "rt", encoding="utf-8")
            filetype = ".gz"
        elif not magic or magic[0] in b">@":
            # Empty input file, or FastA/FastQ input file.
            wrapper = TextIOWrapper(rawfile, encoding="utf-8")
            filetype = ""
        else:
            raise ValueError("Input file is not a GZip, FastQ or FastA file")

        with wrapper as infile:
            firstchar = infile.read(1)
            if not firstchar:
                try:
                    filetype = (".fa" if ".fa" in rawfile.name else ".fq") + filetype
                except AttributeError:
                    filetype = ".fq" + filetype
                yield filetype, []
            elif firstchar == ">":
                filetype = ".fa" + filetype
                if as_fastq:
                    yield filetype, read_fasta_as_fastq(infile)
                yield filetype, read_fasta(infile)
            elif firstchar == "@":
                filetype = ".fq" + filetype
                yield filetype, read_fastq(infile)
            else:
                raise ValueError("Input GZip file does not contain a FastQ or FastA file")
#parse_reads


@contextmanager
def parse_read_pairs(filename1, filename2, as_fastq=False):
    """
    Context manager opening the given files for reading, then yielding a
    generator that produces two-tuples of read records, each as produced
    by the reads_generator created by parse_reads(filename, as_fastq).

    If filename2 is None, filename1 is treated as an interleaved file.
    If filename2 is not None, it must have the same format as filename1.
    """
    with parse_reads(filename1, as_fastq) as (read1_format, read1_generator):
        if filename2 is None:
            read2_manager = nullcontext((read1_format, read1_generator))
        else:
            read2_manager = parse_reads(filename2, as_fastq)
        with read2_manager as (read2_format, read2_generator):
            if read1_format.removesuffix(".gz") != read2_format.removesuffix(".gz"):
                raise ValueError("R1 and R2 files should be the same format (FastA/FastQ)!")

            # TODO: When we drop support for Py3.9, add strict=True.
            yield zip(read1_generator, read2_generator)
#parse_read_pairs


def write_sequence_record(outfile, record):
    """
    Write a tuple of (header, sequence) to a FastA stream or
    a tuple of (header, sequence, quality) to a FastQ stream.
    """
    outfile.write((">%s\n%s\n" if len(record) == 2 else "@%s\n%s\n+\n%s\n") % record)
#write_sequence_record
