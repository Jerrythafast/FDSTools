Forensic DNA Sequencing Tools
=============================

Tools for filtering and interpretation of Next Generation Sequencing data of
forensic DNA samples. To obtain a list of included tools with a brief
description of each tool, run:

    ``fdstools -h``

For a complete description of a specific tool and its command line arguments,
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
    - Updated bundled JavaScript library Vega to v2.5.0
    - Updated bundled JavaScript library D3 to v3.5.12
    - Includes Allelefinder v1.0.0
    - Includes BGCorrect v1.0.0
    - Includes BGEstimate v1.0.0
    - Includes BGHomRaw v1.0.0
    - Includes BGHomStats v1.0.0
    - Includes BGMerge v1.0.0
    - Includes BGPredict v1.0.0
    - Includes Blame v1.0.0
    - Includes FindNewAlleles v1.0.0
    - Includes Libconvert v1.0.0
    - Includes Samplestats v1.0.0
    - Includes Seqconvert v1.0.0
    - Includes Stuttermark v1.5.0
    - Includes Stuttermodel v1.0.0
    - Includes TSSV v1.0.0
    - Includes Vis v1.0.0
    - Includes Allelevis v1.0.0beta1
    - Includes BGRawvis v1.0.0
    - Includes Profilevis v1.0.0
    - Includes Samplevis v2.0.0
    - Includes Stuttermodelvis v1.0.0beta1

v0.0.2
    - Added global -d/--debug switch
    - Includes Stuttermark v1.4

v0.0.1
    - Initial version
    - Includes Stuttermark v1.3


Allelefinder
~~~~~~~~~~~~
v1.0.0
    - Initial version


BGCorrect
~~~~~~~~~
v1.0.0
    - Initial version


BGEstimate
~~~~~~~~~~
v1.0.0
    - Initial version


BGHomRaw
~~~~~~~~
v1.0.0
    - Initial version


BGHomStats
~~~~~~~~~~
v1.0.0
    - Initial version


BGMerge
~~~~~~~
v1.0.0
    - Initial version


BGPredict
~~~~~~~~~
v1.0.0
    - Initial version


Blame
~~~~~
v1.0.0
    - Initial version


FindNewAlleles
~~~~~~~~~~~~~~
v1.0.0
    - Initial version


Libconvert
~~~~~~~~~~
v1.0.0
    - Initial version


Samplestats
~~~~~~~~~~~
v1.0.0
    - Initial version


Seqconvert
~~~~~~~~~~
v1.0.0
    - Initial version


Stuttermark
~~~~~~~~~~~
v1.5.0
    - Changed column names 'name' and 'allele' to 'marker' and 'sequence',
      respectively. WARNING: Stuttermark is now INCOMPATIBLE with output
      from TSSV_ but made compatible with TSSV-Lite and the new, bundled TSSV
      tool instead.

v1.4.0
    - Stuttermark now accepts raw sequences and allele names as input, which
      are automatically rewritten as TSSV-style sequences using a specified
      library file
    - The 'name' column is now optional

v1.3.0
    - First version of Stuttermark to be included in ``fdstools``
    - Fixed crash that occurred when an empty allele (e.g., a primer dimer)
      was encountered
    - Stuttermark now prints a warning if an allele is encountered that is
      not a TSSV_-style sequence

v1.2.0
    - All settings are now available from the command line
    - Use 1-based indexing in ``STUTTER`` annotations

v1.1.0
    - Stuttermark now accepts file names and the minimum number of reads to
      evaluate as command line arguments

v1.0.0
    - Initial version


Stuttermodel
~~~~~~~~~~~~
v1.0.0
    - Initial version


TSSV
~~~~
v1.0.0
    - Initial version


Vis
~~~
v1.0.0
    - Initial version


Allelevis
~~~~~~~~~
v1.0.0beta1
    - Initial version


BGRawvis
~~~~~~~~
v1.0.0
    - Initial version


Profilevis
~~~~~~~~~~
v1.0.0
    - Initial version


Samplevis
~~~~~~~~~
v2.0.0
    - Initial version


Stuttermodelvis
~~~~~~~~~~~~~~~
v1.0.0beta1
    - Initial version


.. _TSSV: https://pypi.python.org/pypi/tssv/
