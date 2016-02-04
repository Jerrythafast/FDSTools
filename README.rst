Forensic DNA Sequencing Tools
=============================

Tools for filtering and interpretation of Next Generation Sequencing data of
forensic DNA samples. To obtain a list of included tools with a brief
description of each tool, run:

    ``fdstools -h``

For a complete description of the command line arguments of a specific tool,
run:

    ``fdstools -h TOOLNAME``


Installation
------------

The recommended way to install FDSTools is by using the ``pip`` package
installer. If you have ``pip`` installed, you can easily install FDSTools by
typing:

    ``pip install fdstools``

Alternatively, FDSTools can be installed by running:

    ``python setup.py install``


FDSTools Changelog
------------------
v0.0.3
    - Includes Stuttermark v1.5

v0.0.2
    - Added global -d/--debug switch
    - Includes Stuttermark v1.4

v0.0.1
    - Initial version
    - Includes Stuttermark v1.3


Stuttermark
-----------

Mark potential stutter products by assuming a fixed maximum percentage of
stutter product vs the parent allele.

Input
    Tab-seperated file with the following columns:
        - 'allele': the allele name, as a TSSV_-style sequence, e.g.,
          "``AGAT(12)TGAT(4)``" (required)
        - 'total': the total number of reads (required)
        - 'name': the name of the marker (optional)

    This format is compatible with 'knownalleles.csv' files created by TSSV_.
    If raw sequences or allele names are provided, Stuttermark can convert
    those to TSSV-style sequences automatically if a library file is given as
    well.

Output
    The same file, with an additional column (named 'annotation' by default).
    The new column contains '``STUTTER``' for possible stutter products, or
    '``ALLELE``' otherwise. Lines that were not evaluated are annotated as
    '``UNKNOWN``'. The ``STUTTER`` annotation contains additional information.
    For example,

        ``STUTTER:146.6x1(2-1):10.4x2(2-1x9-1)``

    This is a stutter product for which at most 146.6 reads have come from the
    first sequence in the output file ("``146.6x1``") and at most 10.4 reads
    have come from the second sequence in the output file ("``10.4x2``"). This
    sequence differs from the first sequence in the output file by a loss of
    one repeat of the second repeat block ("``2-1``") and it differs from the
    second sequence by the loss of one repeat in the second block *and* one
    repeat in the ninth block ("``2-1x9-1``").


Changelog
~~~~~~~~~

v1.5
    - Changed column names 'name' and 'allele' to 'marker' and 'sequence',
      respectively. WARNING: Stuttermark is now INCOMPATIBLE with output
      from TSSV but made compatible with TSSV-Lite instead.

v1.4
    - Stuttermark now accepts raw sequences and allele names as input, which
      are automatically rewritten as TSSV-style sequences using a specified
      library file
    - The 'name' column is now optional

v1.3
    - First version of Stuttermark to be included in ``fdstools``
    - Fixed crash that occurred when an empty allele (e.g., a primer dimer)
      was encountered
    - Stuttermark now prints a warning if an allele is encountered that is
      not a TSSV_-style sequence

v1.2
    - All settings are now available from the command line
    - Use 1-based indexing in ``STUTTER`` annotations

v1.1
    - Stuttermark now accepts file names and the minimum number of reads to
      evaluate as command line arguments

v1.0
    - Initial version


Libconvert
----------

Convert between TSSV (tab-separated) and FDSTools (ini-style) library formats.


Seqconvert
----------

Convert between raw sequences, TSSV-style sequences, and allele names.


Allelefinder
------------

Find true alleles in a single-person reference sample.


.. _TSSV: https://pypi.python.org/pypi/tssv/
