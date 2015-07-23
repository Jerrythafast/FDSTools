#!/usr/bin/env python
"""
Find true alleles in a single-person reference sample.
"""
import argparse
import sys

from ..lib import get_column_ids, pos_int_arg

__version__ = "0.1dev"


_DEF_MIN_ALLELE_PCT = 30.0
_DEF_MAX_NOISE_PCT = 10.0
_DEF_MAX_ALLELES = 2


def find_alleles(infile, outfile, min_allele_pct, max_noise_pct,
                 stuttermark_column, max_alleles):
    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_total, colid_allele, colid_name = get_column_ids(column_names,
        "total", "allele", "name")

    # Also get stuttermark column if we have one.
    if stuttermark_column is not None:
        colid_stuttermark = get_column_ids(column_names, stuttermark_column)

    highest_noise = {}
    highest_allele = {}
    alleles = {}
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if (stuttermark_column is not None and
                not line[colid_stuttermark].startswith("ALLELE")):
            continue

        marker = line[colid_name]
        allele = line[colid_allele]
        reads = int(line[colid_total])

        if marker in alleles:
            if reads > highest_allele[marker]:
                # New highest allele!
                highest_allele[marker] = reads
                for allele in alleles[marker]:
                    if (alleles[marker][allele] <
                            marker_max[marker] * (min_allele_pct/100.)):
                        if alleles[marker][allele] > highest_noise[marker]:
                            highest_noise[marker] = alleles[marker][allele]
                        del alleles[marker][allele]
            elif reads >= highest_allele[marker]*(min_allele_pct/100.):
                # New secundary allele!
                alleles[marker][allele] = reads
            elif reads >= highest_noise[marker]:
                # New highest noise!
                highest_noise[marker] = reads
        else:
            alleles[marker] = {allele: reads}
            highest_allele[marker] = reads
            highest_noise[marker] = 0

    outfile.write("\t".join(["marker", "allele"]) + "\n")
    for marker in alleles:
        if len(alleles[marker]) > max_alleles:
            allele_order = sorted(alleles[marker],
                                  key=lambda x: -alleles[marker][x])
            highest_noise[marker] = alleles[marker][allele_order[max_alleles]]
            alleles[marker] = {x: alleles[marker][x]
                               for x in allele_order[:max_alleles]}
        for allele in alleles[marker]:
            outfile.write("\t".join([marker, allele]) + "\n")
        if highest_noise[marker] > highest_allele[marker]*(max_noise_pct/100.):
            outfile.write("\t".join([marker, "NOISY"]) + "\n")
#find_alleles


def add_arguments(parser):
    parser.add_argument('infile', nargs='?', metavar="IN", default=sys.stdin,
        type=argparse.FileType('r'),
        help="the CSV data file to process (default: read from stdin)")
    parser.add_argument('outfile', nargs='?', metavar="OUT",
        default=sys.stdout, type=argparse.FileType('w'),
        help="the file to write the output to (default: write to stdout)")
    parser.add_argument('-m', '--min-allele-pct', metavar="N", type=float,
        default=_DEF_MIN_ALLELE_PCT,
        help="call heterozygous if the second allele is at least this "
             "percentage of the highest allele (default: %(default)s)")
    parser.add_argument('-M', '--max-noise-pct', metavar="N", type=float,
        default=_DEF_MAX_NOISE_PCT, help="output additional \"NOISY\" allele "
             "if the highest non-allelic sequence is at least this "
             "percentage of the highest allele (default: %(default)s)")
    parser.add_argument('-a', '--max-alleles', metavar="N", type=pos_int_arg,
        default=_DEF_MAX_ALLELES, help="allow no more than this number of "
             "alleles per marker (default: %(default)s)")
    parser.add_argument('-c', '--stuttermark-column', metavar="COLNAME",
        default=None,
        help="name of column with Stuttermark output; if specified, sequences "
             "for which the value in this column does not start with ALLELE "
             "are ignored")
#add_arguments


def run(args):
    if args.infile.isatty() and args.outfile.isatty():
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    find_alleles(args.infile, args.outfile, args.min_allele_pct,
                 args.max_noise_pct, args.stuttermark_column, args.max_alleles)
#run


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    try:
        add_arguments(parser)
        run(parser.parse_args())
    except OSError as error:
        parser.error(error)
#main


if __name__ == "__main__":
    main()
