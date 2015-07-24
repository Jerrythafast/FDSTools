#!/usr/bin/env python
"""
Find true alleles in reference samples and detect possible
contaminations.
"""
import argparse
import sys
import re

from ..lib import get_column_ids, pos_int_arg, get_tag

__version__ = "0.1dev"


# Default values for parameters are specified below.

# Default regular expression to capture sample tags in file names.
# This value can be overridden by the -e command line option.
_DEF_TAG_EXPR = "^(.+?)(?:\.[^.]+)?$"

# Default formatting template to write sample tags.
# This value can be overridden by the -f command line option.
_DEF_TAG_FORMAT = "\\1"

# Default minimum number of reads required for the highest allele.
# This value can be overridden by the -n command line option.
_DEF_MIN_READS = 50

# Default minimum number of reads required for an allele to be called,
# as a percentage of the number of reads of the highest allele.
# This value can be overridden by the -m command line option.
_DEF_MIN_ALLELE_PCT = 30.0

# Default maximum amount of noise to allow, as a percentage of the
# number of reads of the highest allele of each marker.  If any noise
# (i.e., non-allelic sequences) above this threshold are detected, the
# sample is considered 'noisy' for this marker.
# This value can be overridden by the -M command line option.
_DEF_MAX_NOISE_PCT = 10.0

# Default maximum number of alleles to expect for each marker.
# This value can be overridden by the -a command line option.
_DEF_MAX_ALLELES = 2

# Default maximum number of noisy markers allowed per sample.
# This value can be overridden by the -x command line option.
_DEF_MAX_NOISY = 2


def find_alleles(filelist, reportfile, tag_expr, tag_format, min_reads,
                 min_allele_pct, max_noise_pct, max_alleles, max_noisy,
                 stuttermark_column):

    print("\t".join(["sample", "marker", "total", "allele"]))
    for infile in filelist:
        tag = get_tag(infile.name, tag_expr, tag_format)
        find_alleles_sample(infile, reportfile, tag, min_reads, min_allele_pct,
            max_noise_pct, max_alleles, max_noisy, stuttermark_column)
#find_alleles


def find_alleles_sample(infile, reportfile, tag, min_reads, min_allele_pct,
                        max_noise_pct, max_alleles, max_noisy,
                        stuttermark_column):

    # Get column numbers.
    column_names = infile.readline().rstrip("\r\n").split("\t")
    colid_total, colid_allele, colid_name = get_column_ids(column_names,
        "total", "allele", "name")

    # Also get stuttermark column if we have one.
    if stuttermark_column is not None:
        colid_stuttermark = get_column_ids(column_names, stuttermark_column)

    top_noise = {}
    top_allele = {}
    alleles = {}
    for line in infile:
        line = line.rstrip("\r\n").split("\t")
        if (stuttermark_column is not None and
                not line[colid_stuttermark].startswith("ALLELE")):
            continue

        marker = line[colid_name]
        allele = line[colid_allele]
        reads = int(line[colid_total])

        if marker not in alleles:
            alleles[marker] = {allele: reads}
            top_allele[marker] = reads
            top_noise[marker] = ["-", 0]
        else:
            if reads > top_allele[marker]:
                # New highest allele!
                top_allele[marker] = reads
                for allelex in alleles[marker].keys():
                    if (alleles[marker][allelex] <
                            top_allele[marker] * (min_allele_pct/100.)):
                        if alleles[marker][allelex] > top_noise[marker][1]:
                            top_noise[marker] = [
                                allelex, alleles[marker][allelex]]
                        del alleles[marker][allelex]
                alleles[marker][allele] = reads
            elif reads >= top_allele[marker]*(min_allele_pct/100.):
                # New secundary allele!
                alleles[marker][allele] = reads
            elif reads >= top_noise[marker][1]:
                # New highest noise!
                top_noise[marker] = [allele, reads]

    # Find and eliminate noisy markers in this sample first.
    noisy_markers = 0
    for marker in alleles:
        if top_allele[marker] < min_reads:
            if reportfile:
                reportfile.write(
                    "Sample %s is not suitable for marker %s:\n"
                    "highest allele has only %i reads\n\n" %
                        (tag, marker, top_allele[marker]))
            alleles[marker] = {}
            continue
        if len(alleles[marker]) > max_alleles:
            allele_order = sorted(alleles[marker],
                                  key=lambda x: -alleles[marker][x])
            top_noise[marker] = [allele_order[max_alleles],
                alleles[marker][allele_order[max_alleles]]]
            alleles[marker] = {x: alleles[marker][x]
                               for x in allele_order[:max_alleles]}
        if top_noise[marker][1] > top_allele[marker]*(max_noise_pct/100.):
            if reportfile:
                reportfile.write(
                    "Sample %s is not suitable for marker %s:\n"
                    "highest non-allele is %.1f%% of the highest allele\n" %
                        (tag, marker,
                        100.*top_noise[marker][1]/top_allele[marker]))
                for allele in sorted(alleles[marker],
                                     key=lambda x: -alleles[marker][x]):
                    reportfile.write("%i\tALLELE\t%s\n" %
                        (alleles[marker][allele], allele))
                reportfile.write("%i\tNOISE\t%s\n\n" %
                    (top_noise[marker][1], top_noise[marker][0]))
            noisy_markers += 1
            alleles[marker] = {}

    # Drop this sample completely if it has too many noisy markers.
    if noisy_markers > max_noisy:
        if reportfile:
            reportfile.write("Sample %s appears to be contaminated!\n\n" % tag)
        return

    # The sample is OK, write out its alleles.
    for marker in alleles:
        for allele in sorted(alleles[marker],
                             key=lambda x: -alleles[marker][x]):
            print("\t".join(
                [tag, marker, str(alleles[marker][allele]), allele]))
#find_alleles_sample


def add_arguments(parser):
    parser.add_argument('filelist', nargs='*', metavar="FILE",
        default=[sys.stdin], type=argparse.FileType('r'),
        help="the data file(s) to process (default: read from stdin)")
    parser.add_argument('-r', '--report', metavar="OUTFILE",
        type=argparse.FileType("w"),
        help="write a report to the given file, detailing possibly "
             "contaminated or otherwise unsuitable samples")
    parser.add_argument('-e', '--tag-expr', metavar="REGEX", type=re.compile,
        default=_DEF_TAG_EXPR,
        help="regular expression that captures (using one or more capturing "
             "groups) the sample tags from the file names; by default, the "
             "entire file name except for its extension (if any) is captured")
    parser.add_argument('-f', '--tag-format', metavar="EXPR",
        default=_DEF_TAG_FORMAT,
        help="format of the sample tags produced; a capturing group reference "
             "like '\\n' refers to the n-th capturing group in the regular "
             "expression specified with -e/--tag-expr (the default of '\\1' "
             "simply uses the first capturing group); with a single sample, "
             "you can enter the samle tag here explicitly")
    parser.add_argument('-n', '--min-reads', metavar="N", type=pos_int_arg,
        default=_DEF_MIN_READS,
        help="require at least this number of reads for the highest allele "
             "(default: %(default)s)")
    parser.add_argument('-m', '--min-allele-pct', metavar="N", type=float,
        default=_DEF_MIN_ALLELE_PCT,
        help="call heterozygous if the second allele is at least this "
             "percentage of the highest allele (default: %(default)s)")
    parser.add_argument('-M', '--max-noise-pct', metavar="N", type=float,
        default=_DEF_MAX_NOISE_PCT,
        help="a sample is considered contaminated/unsuitable for a marker if "
             "the highest non-allelic sequence is at least this percentage of "
             "the highest allele (default: %(default)s)")
    parser.add_argument('-a', '--max-alleles', metavar="N", type=pos_int_arg,
        default=_DEF_MAX_ALLELES,
        help="allow no more than this number of alleles per marker (default: "
             "%(default)s)")
    parser.add_argument('-x', '--max-noisy', metavar="N", type=pos_int_arg,
        default=_DEF_MAX_NOISY,
        help="entirely reject a sample if more than this number of markers "
             "have a high non-allelic sequence (default: %(default)s)")
    parser.add_argument('-c', '--stuttermark-column', metavar="COLNAME",
        default=None,
        help="name of column with Stuttermark output; if specified, sequences "
             "for which the value in this column does not start with ALLELE "
             "are ignored")
#add_arguments


def run(args):
    if args.filelist == [sys.stdin] and sys.stdin.isatty():
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")

    find_alleles(args.filelist, args.report, args.tag_expr, args.tag_format,
                 args.min_reads, args.min_allele_pct, args.max_noise_pct,
                 args.max_alleles, args.max_noisy, args.stuttermark_column)
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
