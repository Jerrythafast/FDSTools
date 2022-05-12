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

import sys

from .seq import SEQ_SPECIAL_VALUES, ensure_sequence_format


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


def parse_allelelist(allelelist, *, convert=None, library=None):
    """Read allele list from open file handle."""
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
        if sample not in alleles:
            alleles[sample] = {}
        if marker not in alleles[sample]:
            alleles[sample][marker] = set()
        alleles[sample][marker].add(allele)
    return alleles
#parse_allelelist


def read_sample_data_file(infile, data, annotation_column=None, seqformat=None, library=None,
                          default_marker=None, drop_special_seq=False, after_correction=False,
                          combine_strands=False, extra_columns=None):
    """Add data from infile to data dict as [marker, sequence]=reads."""
    # TODO: require keywords
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    if column_names == [""]:
        return []  # Empty file.
    colid_sequence = get_column_ids(column_names, "sequence")
    if combine_strands:
        colid_total = None
        numtype = int
        if after_correction:
            colid_total = get_column_ids(column_names, "total_corrected",
                optional=(after_correction != "require"))
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
                "forward_corrected", "reverse_corrected",
                optional=(after_correction != "require"))
        if colid_forward is None:
            colid_forward = get_column_ids(column_names, "forward")
        else:
            numtype = float
        if colid_reverse is None:
            colid_reverse = get_column_ids(column_names, "reverse")
        else:
            numtype = float
    if extra_columns is not None:
        extra_colids = {c: i for c, i in
            ((c, get_column_ids(column_names, c, optional=extra_columns[c]))
                for c in extra_columns)
            if i is not None}

    # Get marker name column if it exists.
    colid_marker = get_column_ids(column_names, "marker", optional=True)

    # Also try to get annotation column if we have one.
    if annotation_column is not None:
        try:
            colid_annotation = get_column_ids(column_names, annotation_column)
        except:
            annotation_column = None

    found_alleles = []
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if drop_special_seq and line[colid_sequence] in SEQ_SPECIAL_VALUES:
            continue
        marker = line[colid_marker] if colid_marker is not None else default_marker
        sequence = line[colid_sequence] if seqformat is None else ensure_sequence_format(
            line[colid_sequence], seqformat, library=library, marker=marker)
        if annotation_column is not None and line[colid_annotation].startswith("ALLELE"):
            found_alleles.append((marker, sequence))
        if combine_strands:
            data[marker, sequence] = [numtype(line[colid_total])]
        else:
            data[marker, sequence] = list(map(numtype, (line[colid_forward], line[colid_reverse])))
        if extra_columns is not None:
            data[marker, sequence].append({c: line[extra_colids[c]] for c in extra_colids})

    return found_alleles
#read_sample_data_file


def get_sample_data(tags_to_files, callback, allelelist=None, annotation_column=None,
                    seqformat=None, library=None, marker=None, homozygotes=False,
                    drop_special_seq=False, after_correction=False, combine_strands=False,
                    extra_columns=None):
    """
    Parse sample data.
    TODO: write full doc, and require keywords...
    """
    for tag in tags_to_files:
        data = {}
        alleles = set()
        for infile in tags_to_files[tag]:
            infile = sys.stdin if infile == "-" else open(infile, "rt", encoding="UTF-8")
            alleles.update(read_sample_data_file(
                infile, data, annotation_column, seqformat, library, marker,
                drop_special_seq, after_correction, combine_strands, extra_columns))
            if infile != sys.stdin:
                infile.close()
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
                for markerx in tuple(allelelist[tag]):
                    if len(allelelist[tag][markerx]) > 1:
                        del allelelist[tag][markerx]
        callback(tag, data)
#get_sample_data


def try_write_pipe(stream, *args):
    """Call stream.write(*args), ignore EPIPE errors."""
    try:
        stream.write(*args)
    except IOError as e:
        if e.errno == EPIPE:
            return
        raise
#try_write_pipe


def print_db(text, debug):
    """Print text if debug is True."""
    if debug:
        print(text)
#print_db
