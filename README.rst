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
v1.0.0rc1
    - Fixed bug that caused variant descriptions in allele names of non-STR
      markers to be prepended with plus signs similar to suffix variants
      in STR markers; when attempting to convert these allele names back to raw
      sequences, FDSTools would crash with an 'Invalid allele name' error
    - Tools that take a list of files as their argument (through the -i option
      or as positionals) now explicitly support '*' and '?' wildcards
    - Includes BGEstimate v1.1.1
    - Includes BGMerge v1.0.2
    - Includes Library v1.0.1
    - Includes Pipeline v1.0.1
    - Includes Stuttermodel v1.1.1
    - Includes Allelevis v2.0.1
    - Includes Samplevis v2.1.1
    - Includes Stuttermodelvis v2.0.1

v0.0.5
    - The Blame tool was removed in favour of BGAnalyse
    - Includes BGAnalyse v1.0.0
    - Includes Libconvert v1.1.0
    - Includes Library v1.0.0
    - Includes Pipeline v1.0.0
    - Includes Samplestats v1.1.0
    - Includes TSSV v1.0.2
    - Includes Vis v1.0.2
    - Includes Allelevis v2.0.0
    - Includes BGAnalysevis v1.0.0
    - Includes BGRawvis v2.0.0
    - Includes Profilevis v2.0.0
    - Includes Samplevis v2.1.0
    - Includes Stuttermodelvis v2.0.0

v0.0.4
    - FDSTools will now print profiling information to stdout when the
      -d/--debug option was specified
    - Fixed bug where specifying '-' as the output filename would be taken
      literally, while it should have been interpreted as 'write to standard
      out' (Affected tools: BGCorrect, Samplestats, Seqconvert, Stuttermark)
    - Added more detailed license information to FDSTools
    - Updated bundled JavaScript library Vega to v2.6.0
    - Updated bundled JavaScript library D3 to v3.5.17
    - Includes BGCorrect v1.0.1
    - Includes BGEstimate v1.1.0
    - Includes BGMerge v1.0.1
    - Includes BGPredict v1.0.1
    - Includes Libconvert v1.0.1
    - Includes Samplestats v1.0.1
    - Includes Seqconvert v1.0.1
    - Includes Stuttermodel v1.1.0
    - Includes TSSV v1.0.1
    - Includes Vis v1.0.1
    - Includes Allelevis v1.0.0beta2
    - Includes BGRawvis v1.0.1
    - Includes Profilevis v1.0.1
    - Includes Samplevis v2.0.1
    - Includes Stuttermodelvis v1.0.0beta2

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


BGAnalyse
~~~~~~~~~
v1.0.0
    - Initial version


BGCorrect
~~~~~~~~~
v1.0.1
    - Added new column 'weight' to the output. The value in this column
      expresses the number of times that the noise profile of that allele
      fitted in the sample.

v1.0.0
    - Initial version


BGEstimate
~~~~~~~~~~
v1.1.1
    - Added option -p/--profiles which can be used to provide a previously
      created background noise profiles file, from which starting values will
      be read instead of assuming zero noise

v1.1.0
    - Added a new option -g/--min-genotypes (default: 3). Only alleles that
      occur in at least this number of unique heterozygous genotypes will be
      considered. This is to avoid 'contamination' of the noise profile of one
      allele with the noise of another. If homozygous samples are available for
      an allele, this filter is not applied to that allele. Setting this option
      to 1 effectively disables it. This option has the same cascading effect
      as the -s/--min-samples option, that is, if one allele does not meet the
      threshold, the samples with this allele are excluded which may cause some
      of the other alleles of these samples to fall below the threshold as
      well.

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
v1.0.2
    - Minor changes to facilitate explicit filename wildcard support

v1.0.1
    - Reduced memory usage

v1.0.0
    - Initial version


BGPredict
~~~~~~~~~
v1.0.1
    - Greatly reduced memory usage.
    - BGPredict will now output nonzero values below the threshold set by
      -n/--min-pct if the predicted noise ratio of the same stutter on the
      other strand is above the threshold. Previously, values below the
      threshold were clipped to zero, which may cause unnecessarily high strand
      bias in the predicted profile.

v1.0.0
    - Initial version


FindNewAlleles
~~~~~~~~~~~~~~
v1.0.0
    - Initial version


Libconvert
~~~~~~~~~~
v1.1.0
    - When converting to FDSTools format, Libconvert automatically creates an
      empty FDSTools library file with the same contents as what would be
      obtained from the new Library tool without arguments.
    - The -a/--aliases option was modified such that it has the same effect as
      the -a/--aliases option of the new Library tool. This means that without
      this option specified, the [aliases] section will not be present in the
      output anymore.
    - The ability of the Libconvert tool to produce an empty FDSTools library
      file if no input file was given has been removed from the documentation
      (but not from the tool itself).

v1.0.1
    - Specifying '-' as the first positional argument to libconvert will now
      correctly interpret this as "read from stdin" instead of throwing a "file
      not found" error (or reading from a file named "-" if it exists)

v1.0.0
    - Initial version


Library
~~~~~~~
v1.0.1
    - Updated some of the comments describing the sections
    - Added proper examples for non-STR markers and aliases

v1.0.0
    - Initial version


Pipeline
~~~~~~~~
v1.0.1
    - Removed checking of the existence of the files specified for the
      -S/--in-samples option; instead, this is left to the downstream tools to
      find out, consistent with how this works with other input file options
    - Only output the running commands if the -d/--debug option was specified

v1.0.0
    - Initial version


Samplestats
~~~~~~~~~~~
v1.1.0
    - Changed default allele calling option thresholds:
        - Changed default value of -m/--min-pct-of-max from 5.0 to 2.0
        - Changed default value of -p/--min-pct-of-sum from 3.0 to 1.5
    - Mentioned allele calling in the tool descriptions

v1.0.1
    - Samplestats will now round to 4 or 5 significant digits if a value is
      above 1000 or 10000, respectively. Previously, this was only done for the
      combined 'Other sequences' values
    - The 'Other sequences' lines will now also include values for
      total_recovery, forward_recovery, and reverse_recovery
    - The total_recovery, forward_recovery, and reverse_recovery columns are no
      longer placed to the left of all the other columns generated by
      Samplestats
    - The help text for Samplestats erroneously listed the X_recovery_pct
      instead of X_recovery
    - Added support for the new 'weight' column produced by BGCorrect when the
      -a/--filter-action option is set to 'combine'

v1.0.0
    - Initial version


Seqconvert
~~~~~~~~~~
v1.0.1
    - Internal naming of the first positional argument was changed from
      'format' to 'sequence-format'. This was done for consistency with the
      -F/--sequence-format option in other tools, giving it the same name in
      Pipeline configuration files.

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
v1.1.1
    - Minor change to internal variant representation

v1.1.0
    - Stuttermodel will now only output a fit for one strand if it could also
      obtain a fit for the other strand (for the same marker, unit, and stutter
      depth). This new behaviour can be disabled with a new -O/--orphans
      option.
    - Fixed bug that caused Stuttermodel to output only the raw data points for
      -1 and +1 stutter when normal output was supressed

v1.0.0
    - Initial version


TSSV
~~~~
v1.0.2
    - Added new option -n/--indel-score which can be used to increase the
      penalty given to insertions and deletions in the flanking sequences
      w.r.t. the penalty given to mismatches.
    - NOTE: Requires TSSV v0.4.0 or newer to be installed.

v1.0.1
    - Renamed the '--is_fastq' option to '--is-fastq', which was the only
      option with an underscore instead of a hyphen in FDSTools
    - Fixed crash that would occur if -F/--sequence-format was set to anything
      other than 'raw'

v1.0.0
    - Initial version


Vis
~~~
v1.0.2
    - Changed default value of -n/--min-abs from 15 to 5
    - Added -I/--input2 option, which allows for specifying a file with raw
      data points for Stuttermodelvis and Profilevis
    - Added support for creating BGAnalysevis visualisations

v1.0.1
    - Added -j/--jitter option for Stuttermodelvis (default: 0.25)
    - Fixed bug where Vis would not allow the -n/--min-abs and the
      -s/--min-per-strand options to be set to 0

v1.0.0
    - Initial version


Allelevis
~~~~~~~~~
v2.0.1
    - Added tooltip support to HTML visualisations

v2.0.0
    - Replaced the simple Options overlay with responsive design options panels
      in HTML visualisations
    - Reduced Vega graph spec complexity by using the new Rank transform to
      position the subgraphs
    - Fixed glitch that caused unnecessary padding around the graph

v1.0.0beta2
    - Fixed potential crash/corruption that could occur with very unfortunate
      combinations of sample names and marker names
    - HTML visualisations made with the -O/--online option of the Vis tool will
      now contain https URLs instead of http
    - Added two more colours to the legend, such that a maximum of 22 markers
      is now supported without re-using colours

v1.0.0beta1
    - Initial version


BGAnalysevis
~~~~~~~~~~~~
v1.0.0
    - Initial version


BGRawvis
~~~~~~~~
v2.0.0
    - Replaced the simple Options overlay with responsive design options panels
      in HTML visualisations
    - Sequences are now sorted by CE allele length when applicable
    - Changed default minimum number of reads from 15 to 5
    - Added marker selection menu for easier filtering

v1.0.1
    - Fixed a JavaScript crash that would occur in HTML visualisations if the
      Marker name filter resulted in an invalid regular expression (e.g., when
      the entered value ends with a backslash)
    - Reduced Vega graph spec complexity by using the new Rank transform to
      position the subgraphs
    - HTML visualisations made with the -O/--online option of the Vis tool will
      now contain https URLs instead of http

v1.0.0
    - Initial version


Profilevis
~~~~~~~~~~
v2.0.0
    - Replaced the simple Options overlay with responsive design options panels
      in HTML visualisations
    - Alleles and sequences are now sorted by CE allele length when applicable
    - Added option to plot BGHomRaw data on top of the profiles
    - Added marker selection menu for easier filtering

v1.0.1
    - Fixed a JavaScript crash that would occur in HTML visualisations if the
      Marker name filter resulted in an invalid regular expression (e.g., when
      the entered value ends with a backslash)
    - Reduced Vega graph spec complexity by using the new Rank transform to
      position the subgraphs.
    - HTML visualisations made with the -O/--online option of the Vis tool will
      now contain https URLs instead of http

v1.0.0
    - Initial version


Samplevis
~~~~~~~~~
v2.1.1
    - Added tooltip support to HTML visualisations
    - The tooltip may include a 'new allele' note if the input sample was
      analysed with FindNewAlleles
    - The allele tables in HTML visualisations will now grow much wider than
      before if the screen (or window) is very narrow
    - Improved line breaking behaviour in the tables in HTML visualisations
    - Improved determination of column widths of the allele tables when
      printing an HTML visualisation
    - When printing an HTML visualisation, the graph and the corresponding
      table of a marker will be kept on the same page in all browsers now
    - Fixed glitch that caused 'Infinity%' or 'NaN%' to be written in some
      cells in the allele tables in HTML visualisations

v2.1.0
    - Changed default minimum number of reads for graph filtering from 15 to 5
    - Changed default table filtering options:
        - Percentage of highest allele per marker changed from 5% to 2%
        - Percentage of the marker's total reads changed from 3% to 1.5%
        - Minimum number of reads in both orientations changed from 0 to 1

v2.0.1
    - Fixed a JavaScript crash that would occur in HTML visualisations if the
      Marker name filter resulted in an invalid regular expression (e.g., when
      the entered value ends with a backslash)
    - Reduced Vega graph spec complexity by using the new Rank transform to
      position the subgraphs
    - Fixed a glitch in HTML visualisations where clicking the 'Truncate
      sequences to' label would select the marker spacing input
    - In HTML visualisations, the 'Notes' table cells with 'BGPredict' in them
      now get a light orange background to warn the user that their background
      profile was computed. If a sequence was explicitly 'not corrected', 'not
      in ref db', or 'corrected as background only', the same colour is used.
    - The message bar at the bottom of Samplevis HTML visualisations will now
      grow no larger than 3 lines. A scroll bar will appear as needed.
    - HTML visualisations made with the -O/--online option of the Vis tool will
      now contain https URLs instead of http

v2.0.0
    - Initial version


Stuttermodelvis
~~~~~~~~~~~~~~~
v2.0.1
    - Changed the unit in the horizontal axis title from 'bp' to 'nt'

v2.0.0
    - Replaced the simple Options overlay with responsive design options panels
      in HTML visualisations
    - Fixed glitch that caused the graphs to be re-rendered twice when loading
      a file by drag-and-drop in HTML visualisations
    - Fixed glitch that made it possible to replace the data that was embedded
      in an HTML visualisation through drag-and-drop
    - Added repeat unit selection menu for easier filtering

v1.0.0beta2
    - HTML visualisations now support drawing raw data points on top of the fit
      functions. The points can be drawn with an adjustable jitter to reduce
      overlap.
    - Fixed a JavaScript crash that would occur in HTML visualisations if the
      Repeat unit or Marker name filter resulted in an invalid regular
      expression (e.g., when the entered value ends with a backslash)
    - Reduced Vega graph spec complexity by using the new Rank transform to
      position the subgraphs.
    - HTML visualisations made with the -O/--online option of the Vis tool will
      now contain https URLs instead of http

v1.0.0beta1
    - Initial version


.. _TSSV: https://pypi.python.org/pypi/tssv/
