Forensic DNA Sequencing Tools
=============================
Tools for filtering and interpretation of Massively Parallel Sequencing data of
forensic DNA samples. To obtain a list of included tools with a brief
description of each tool, run:

    fdstools --help

For a complete description of a specific tool and its command line arguments,
run:

    fdstools --help TOOLNAME


Installation
------------
FDSTools requires Python version 3.5 or later.

The recommended way to install FDSTools is by using the `pip` package
installer. If you have `pip` installed, you can easily install FDSTools by
running the following command:

    pip install -U fdstools

Alternatively, FDSTools can be installed by running:

    python setup.py install


Release Notes
-------------
### Version 2.0.4 (2023-03-02)
Microhaplotype alleles with SNPs adjacent to any of the microhaplotype
positions are now recognised. For example, an allele previously named
"MH_N_123NC>GA" will now be named "MH_G_124C>A" instead.

Greatly improved loading speed for large library files.


### Version 2.0.3 (2022-08-27)
Added NimaGen IDseek OmniSTR built-in library file (ID-OmniSTR).


### Version 2.0.2 (2022-07-25)
Fixed a compatibility issue with Python 3.10, that caused TSSV to crash.
Fixed some issues that could cause crashes in the Pipeline tool.

Corrected the reported range of AmelongeninY in the ForenSeq-A and ForenSeq-B
library files which were introduced in v2.0.1, by adding 1 to the starting
position. The range is now identical to the UAS Flanking Region Report range.
If still you need to work with data previously generated using the ForenSeq-A
or ForenSeq-B library as shipped with v2.0.1, use the Library tool to obtain a
copy of the built-in library file and make sure the range of AmelongeninY
starts at 6869870.


### Version 2.0.1 (2022-05-12)
Added support for microhaplotype targets. For markers configured as such,
the 'allelename' sequence format will look like "MH_AGGTC".

Added ForenSeq DNA Signature Prep Kit library files with ranges specifically
optimized for FDSTools (often longer than UAS Flanking Region Report).


### Version 2.0.0 (2021-07-15)
Transition to Python 3.5, along with a major upgrade to all tools to
provide an overall better experience.

Allele naming is now handled by STRNaming, eliminating the need for
complex library files. Library files of commonly-used kits are built in.

Kits that sequence a target on a single strand (such as the ForenSeq DNA
Signature Prep Kit by Verogen) are now supported.


### Version 1.2.1 (2021-07-15)
This release focuses on finishing support for Python2 before the transition
to Python3. FDSTools will now display the help page if no command is given.
A message about the transition to Python3 in the next version of FDSTools
is added to the command-line help pages. Furthermore, dependency version
numbers have been updated to ensure smooth installation on Python2 as well
as a smooth transition to Python3.


### Version 1.2.0 (2019-03-29)
Major improvements and fixes to the TSSV tool. Most notably, it no longer
relies on the external `tssvl` program because that is no longer
compatible with FDSTools. Furthermore, the new TSSV tool v2.0.0 comes with
a major performance upgrade and has some updated command-line arguments.

This release also fixes an issue in Samplestats and adds the ability to
apply graph filtering before noise correction in Samplevis, making the
effects of noise correction more apparent.


### Version 1.1.1 (2017-03-15)
Fixeds incorrect calculation of tLeft, fLeft, rLeft, tRight and fRight
columns in the report output file of TSSV, when -T/--num-threads was set to
2 or higher. The primary output was unaffected.


### Version 1.1.0 (2017-03-14)
In STR allele names for sequences that don't exactly match the description
given in the library file, no more insertions are produced at the end of
the prefix or the beginning of the suffix, in favour of extra STR blocks.

Empty input files and broken pipelines are now handled gracefully across
all tools. Specifically, an empty input file is now treated as if the
expected columns existed, but no lines of actual data were present. This
greatly helps in tracking down issues in pipelines involving multiple
tools, as tools will now shutdown gracefully if an upstream tool fails to
write output. Only the failing tool will output an error.

Furthermore, a new option has been added to the TSSV tool, enabling
multithreading support. This can greatly reduce analysis time by using
more (or all) cores of the system's processor simultaneously.

Finally, various small bugs and glitches were fixed.


### Version 1.0.1 (2016-12-21)
FDSTools library files may now contain IUPAC ambiguous bases in the prefix
prefix and suffix sequences of STR markers (except the first sequence, as
it is used as the reference). Additionally, optional bases may be
represented by lowercase letters.

An option was added to the Pipeline tool to skip running Allelefinder,
using a user-supplied allele list file instead. Multiple options have been
added to the Vis tool and some have been regrouped to more easily find the
option you are looking for.

It is now possible to save the a Samplevis HTML visualisation after having
made changes, preserving the changes made.

And various minor bug fixes and improvements throughout.


### Version 1.0.0 (2016-10-03)
Fixed an issue with variant descriptions in allele names of non-STR markers
that made it impossible to convert those back to raw sequences.

Added various useful options. Most notably, Samplevis now displays a
tooltip when the mouse pointer is over an allele, providing various details
about that allele.

And various minor bug fixes.


### Version 0.0.5 (2016-09-06)
Added the Library tool, for creating a template library file that includes
helpful commentary and examples to get new users started. Creating an empty
library file used to be a somewhat confusing option in the Libconvert tool.
Also, the Blame tool was replaced with the more advanced BGAnalyse tool.

Added the Pipeline tool, which implements some ready-made pipelines
involving most of the other tools in FDSTools. Three pipelines are
provided: one for noise reference sample analysis, one for case sample
analysis, and one for generating a background noise database from the
reference samples.

In Samplestats, the default allele calling option thresholds have changed:
    - Changed default value of -m/--min-pct-of-max from 5.0 to 2.0
    - Changed default value of -p/--min-pct-of-sum from 3.0 to 1.5

The TSSV tool was updated with an option to increase the penalty given to
insertions and deletions in the flanking sequences. It now requires TSSV
version 0.4.0 to be installed.

Various upgrades to visualisations, bringing a new responsive design to all
HTML visualisations and fixing various issues.


### Version 0.0.4 (2016-07-26)
Improved debugging: FDSTools will now print profiling information to stdout
when the -d/--debug option was specified. Also, all tools now correctly
interpret '-' as the output filename as 'write to standard out'.

BGEstimate has gained a new option to require a minimum number of unique
genotypes in which a specific allele must have been seen before it will be
considered for noise estimation. This is to avoid 'contamination' of the
noise profile of one allele with the noise of another. If homozygous
samples are available for an allele, this filter is not applied to that
allele.

Reduced the memory usage of BGPredict and BGMerge. Also, BGPredict will now
output nonzero values below the threshold set by -n/--min-pct if the
predicted noise ratio of the same stutter on the other strand is above the
threshold. Previously, values below the threshold were clipped to zero,
which may cause unnecessarily high strand bias in the predicted profile.
Similarly, by default Stuttermodel will no longer output a fit on one
strand if no fit could be optained on the other strand.

Changes have been made to rounding and column order in Samplestats.

Various minor fixes and enhancements have been made, mostly to the
visualisations.


### Version 0.0.3 (2016-02-02)
First version of FDSTools with all strings attached. Introduces 15 new tools
and five visualisations.

In Stuttermark, the column names 'name' and 'allele' have been changed to
'marker' and 'sequence', respectively, reflecting those of all the other
tools. WARNING: Stuttermark is now INCOMPATIBLE with output from TSSV, but
made compatible with TSSV-Lite and the new, bundled TSSV tool instead.


### Version 0.0.2 (2015-07-23)
Added a new global option: -d/--debug. This option disables the suppression
of technical details that would normally be visible when an error occurs.

Stuttermark now accepts raw sequences and allele names as input, which are
automatically rewritten as TSSV-style sequences using a specified library
file. Also, the 'name' column is now optional.


### Version 0.0.1 (2015-07-02)
Initial version of FDSTools, featuring a single tool: Stuttermark v1.3.
