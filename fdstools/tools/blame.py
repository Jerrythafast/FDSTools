#!/usr/bin/env python
"""
Find dirty reference samples or recurring contaminating alleles.

Match background noise profiles (as obtained from bgestimate) to the
database of reference samples from which they were obtained to detect
samples with unexpectedly high additional alleles.  The 'amount' column
in the output lists the number of noise/contaminant reads per true
allelic read.

(The tool is called 'blame' because with default settings it might
produce a DNA profile of an untidy laboratory worker.)
"""
import argparse
#import numpy as np  # Only imported when actually running this tool.

from ..lib import pos_int_arg, add_input_output_args, get_input_output_files,\
                  add_allele_detection_args, nnls, ensure_sequence_format,\
                  parse_allelelist, load_profiles, add_sequence_format_args,\
                  parse_library, get_sample_data

__version__ = "0.1dev"


# Default values for parameters are specified below.

# Default number of results per marker.
# This value can be overridden by the -n command line option.
_DEF_NUM = 5


def add_sample_data(data, sample_data, sample_tag, alleles):
    # Add this sample to the data.
    use_markers = set()
    for marker in alleles:
        genotype = [data[marker]["seqs"].index(x)
                    if x in data[marker]["seqs"] and
                        data[marker]["seqs"].index(x) < data[marker]["n"]
                    else -1
                    for x in alleles[marker]]
        if -1 in genotype:
            # Don't have a profile for all of this sample's alleles.
            continue
        data[marker]["samples_forward"].append([0] * data[marker]["m"])
        data[marker]["samples_reverse"].append([0] * data[marker]["m"])
        data[marker]["genotypes"].append(genotype)
        data[marker]["sample_tags"].append(sample_tag)
        use_markers.add(marker)

    for marker, allele in sample_data:
        if marker not in use_markers:
            continue
        try:
            i = data[marker]["seqs"].index(allele)
        except ValueError:
            continue
        data[marker]["samples_forward"][-1][i] = sample_data[marker, allele][0]
        data[marker]["samples_reverse"][-1][i] = sample_data[marker, allele][1]
#add_sample_data


def blame(samples_in, outfile, allelefile, annotation_column, mode,
          profilefile, num, seqformat, libfile, marker):
    import numpy as np
    library = parse_library(libfile) if libfile else None
    allelelist = {} if allelefile is None \
                    else parse_allelelist(allelefile, "raw", library)
    data = load_profiles(profilefile, library)
    if marker:
        data = {marker: data[marker]} if marker in data else {}
    for marker in data:
        data[marker]["samples_forward"] = []
        data[marker]["samples_reverse"] = []
        data[marker]["genotypes"] = []
        data[marker]["sample_tags"] = []

    # Read sample data.
    get_sample_data(
        samples_in,
        lambda tag, sample_data: add_sample_data(
            data, sample_data, tag,
            {m: allelelist[tag][m] for m in data if m in allelelist[tag]}),
        allelelist, annotation_column, "raw", library)

    outfile.write("\t".join(["marker",
                     "allele" if mode == "common" else "sample",
                     "amount"]) + "\n")
    for marker in data:
        if not data[marker]["sample_tags"]:
            continue

        # Estimate which alleles are present in the samples.
        P1 = np.matrix(data[marker]["forward"])
        P2 = np.matrix(data[marker]["reverse"])
        C1 = np.matrix(data[marker]["samples_forward"])
        C2 = np.matrix(data[marker]["samples_reverse"])
        A = nnls(np.hstack([P1, P2]).T, np.hstack([C1, C2]).T).T

        # Zero out the true alleles in each sample and scale the others.
        for i in range(len(data[marker]["genotypes"])):
            indices = data[marker]["genotypes"][i]
            scale = A[i, indices].sum()
            A[i, indices] = 0
            A[i, :] /= scale

        if mode == "common":
            # The columns with the highest means correspond to the most
            # common contaminants.
            A = A.mean(0).flat
            for i in np.argsort(A)[:-num-1:-1]:  # Print the top N.
                if A[i] == 0:
                    break
                outfile.write("\t".join([marker, ensure_sequence_format(
                    data[marker]["seqs"][i], seqformat, library=library,
                    marker=marker), str(A[i])]) + "\n")
        else:
            # The rows with the highest maxima/sums correspond to the
            # samples with the highest amounts of contaminant/s.
            A = A.max(1).flat if mode == "highest" else A.sum(1).flat
            for i in np.argsort(A)[:-num-1:-1]:  # Print the top N.
                if A[i] == 0:
                    break
                outfile.write("\t".join(
                    [marker, data[marker]["sample_tags"][i], str(A[i])])+"\n")
#blame


def add_arguments(parser):
    parser.add_argument('profiles', metavar="PROFILES",
        type=argparse.FileType('r'),
        help="file containing background noise profiles to match")
    parser.add_argument("-m", "--mode", metavar="MODE", default="common",
        choices=("common", "highest", "dirtiest"),
        help="controls what kind of information is printed; 'common' (the "
             "default) prints the top N contaminant alleles per marker, "
             "'highest' prints the top N samples with the highest single "
             "contaminant per marker, and 'dirtiest' prints the top N samples "
             "with the highest total amount of contaminants per marker")
    add_input_output_args(parser)
    add_allele_detection_args(parser)
    filtergroup = parser.add_argument_group("filtering options")
    filtergroup.add_argument('-n', '--num', metavar="N", type=pos_int_arg,
        default=_DEF_NUM,
        help="print the top N results per marker (default: %(default)s)")
    filtergroup.add_argument('-M', '--marker', metavar="MARKER",
        help="work only on MARKER")
    add_sequence_format_args(parser, "raw")
#add_arguments


def run(args):
    files = get_input_output_files(args)
    if not files:
        raise ValueError("please specify an input file, or pipe in the output "
                         "of another program")
    blame(files[0], files[1], args.allelelist, args.annotation_column,
          args.mode, args.profiles, args.num, args.sequence_format,
          args.library, args.marker)
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
