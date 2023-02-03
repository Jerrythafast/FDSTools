FDSTools Changelog
==================
### Version 2.0.4
- Microhaplotype alleles with SNPs adjacent to any of the microhaplotype
  positions are now recognised.
- Greatly improved loading speed for large library files.
- Updated dependency requirement (strnaming~=1.1.4).


### Version 2.0.3
- Added NimaGen IDseek OmniSTR built-in library file (ID-OmniSTR).
- Updated dependency requirement (strnaming~=1.1.3).


### Version 2.0.2
- Fixed a compatibility issue with Python 3.10, that caused TSSV to crash.
- Corrected AmelogeninY range (starting position +1) in ForenSeqA and ForenSeqB
  built-in libraries. Range is now identical to the UAS Flanking Region Report.
- Includes [Pipeline v1.1.1](#Pipeline-111)


### Version 2.0.1
- Added support for microhaplotype targets. For markers configured as such,
  the 'allelename' sequence format will look like "MH_AGGTC".
- Added the ForenSeqA and ForenSeqB built-in libraries, which use ranges
  optimized for FDSTools (often longer than UAS Flaking Region Report).
- Updated dependency requirement (strnaming~=1.1.2).
- Includes [FindNewAlleles v1.1.1](#FindNewAlleles-111)
- Includes [Library v1.1.1](#Library-111)
- Includes [Samplestats v1.3.1](#Samplestats-131)
- Includes [Stuttermark v1.6.1](#Stuttermark-161)
- Includes [TSSV v2.1.1](#TSSV-211)


### Version 2.0.0
- FDSTools can now be run using `python -m fdstools` too, which is helpful
  when the regular `fdstools` command is not available on the PATH.
- FDSTools now contains built-in libraries for commonly-used kits, which can be
  accessed by specifying a predefined library name instead of a path to a file.
- Directly supplying a legacy TSSV library (tab-separated format) is no longer
  supported; instead, the Libconvert tool should explicitly be used to obtain
  a library file in the native format of FDSTools.
- Specifying multiple prefix and suffix sequences for one marker is no longer
  supported; users migrating their library files from v1.x to v2.0 will have
  to make sure only the first prefix or suffix sequence is retained.
- In library files, a semicolon (;) can now be used to add a comment to the end
  of the line. As a side-effect, it is no longer possible to use semicolons as
  value separators.
- Aliases, custom allele names for specific sequences, have been removed.
- FDSTools will now always use UTF-8 encoding for file I/O.
- FDSTools will now display a 'Failed to load X' message if importing tool X
  failed. The other tools will remain available.
- FDSTools will now display a 'Failed to configure X' message if configuring
  tool X failed. The other tools will remain available, and more information
  about the error can be obtained by running the troubled tool in debug mode.
- FDSTools will now display detailed error information when the -d/--debug
  argument is given, if an unexpected error occurs during argument parsing.
- Updated dependency requirements (python>=3.5, numpy>=1.17 and strnaming==1.1.*).
- Includes [Allelefinder v1.1.0](#Allelefinder-110)
- Includes [BGAnalyse v1.1.0](#BGAnalyse-110)
- Includes [BGCorrect v1.1.0](#BGCorrect-110)
- Includes [BGEstimate v1.2.0](#BGEstimate-120)
- Includes [BGHomRaw v1.1.0](#BGHomRaw-110)
- Includes [BGHomStats v1.1.0](#BGHomStats-110)
- Includes [BGMerge v1.1.0](#BGMerge-110)
- Includes [BGPredict v1.1.0](#BGPredict-110)
- Includes [FindNewAlleles v1.1.0](#FindNewAlleles-110)
- Includes [Libconvert v1.2.0](#Libconvert-120)
- Includes [Library v1.1.0](#Library-110)
- Includes [Pipeline v1.1.0](#Pipeline-110)
- Includes [Samplestats v1.3.0](#Samplestats-130)
- Includes [Seqconvert v1.1.0](#Seqconvert-110)
- Includes [Stuttermark v1.6.0](#Stuttermark-160)
- Includes [Stuttermodel v1.2.0](#Stuttermodel-120)
- Includes [TSSV v2.1.0](#TSSV-210)
- Includes [Vis v1.1.0](#Vis-110)
- Includes [Profilevis v2.0.2](#Profilevis-202)
- Includes [Samplevis v2.3.0](#Samplevis-230)
- Includes [Stuttermodelvis v2.0.4](#Stuttermodelvis-204)


### Version 1.2.1
- FDSTools will now display the help page if no command is given.
- The command-line help pages and debug output will now display a message
  that FDSTools v2 is available and how to get it.
- When running setup.py on python3, a message is displayed that this
  version is only compatible with python2 but a newer is available.
- Updated dependency requirements (python>=2.7.9,<3 and numpy<1.17).
- Includes [Samplevis v2.2.2](#Samplevis-222)


### Version 1.2.0
- Includes [Pipeline v1.0.4](#Pipeline-104)
- Includes [Samplestats v1.2.0](#Samplestats-120)
- Includes [TSSV v2.0.0](#TSSV-200)
- Includes [Samplevis v2.2.1](#Samplevis-221)


### Version 1.1.1
- Includes [TSSV v1.1.1](#TSSV-111)


### Version 1.1.0
- Allele name heuristics: don't produce insertions at the end of the prefix
  or at the beginning of the suffix; just include extra STR blocks.
- FDSTools will no longer crash with a 'column not found' error when
  an input file is empty. This situation is now treated as if the
  expected columns existed, but no lines of actual data were present.
  This greatly helps in tracking down issues in pipelines involving
  multiple tools, as tools will now shutdown gracefully if an upstream
  tool fails to write output.
- Includes [Allelefinder v1.0.1](#Allelefinder-101)
- Includes [BGAnalyse v1.0.1](#BGAnalyse-101)
- Includes [BGCorrect v1.0.2](#BGCorrect-102)
- Includes [BGEstimate v1.1.2](#BGEstimate-112)
- Includes [BGHomRaw v1.0.1](#BGHomRaw-101)
- Includes [BGHomStats v1.0.1](#BGHomStats-101)
- Includes [BGMerge v1.0.3](#BGMerge-103)
- Includes [BGPredict v1.0.2](#BGPredict-102)
- Includes [FindNewAlleles v1.0.1](#FindNewAlleles-101)
- Includes [Libconvert v1.1.2](#Libconvert-112)
- Includes [Library v1.0.3](#Library-103)
- Includes [Pipeline v1.0.3](#Pipeline-103)
- Includes [Samplestats v1.1.1](#Samplestats-111)
- Includes [Seqconvert v1.0.2](#Seqconvert-102)
- Includes [Stuttermark v1.5.1](#Stuttermark-151)
- Includes [Stuttermodel v1.1.2](#Stuttermodel-112)
- Includes [TSSV v1.1.0](#TSSV-110)
- Includes [Vis v1.0.4](#Vis-104)
- Includes [BGRawvis v2.0.1](#BGRawvis-201)
- Includes [Profilevis v2.0.1](#Profilevis-201)
- Includes [Samplevis v2.2.0](#Samplevis-220)
- Includes [Stuttermodelvis v2.0.3](#Stuttermodelvis-203)


### Version 1.0.1
- Fixed crash that occurred when using the -i option to run the same
  command on multiple input files.
- The 'usage' line now always starts with 'fdstools', even if FDSTools was
  invoked through some other command (e.g. on Windows, FDSTools gets
  invoked through a file called 'fdstools-script.py').
- Fixed bug with the -d/--debug option being ignored if placed before the
  tool name on systems running Python 2.7.9 or later.
- FDSTools library files may now contain IUPAC ambiguous bases in the
  prefix and suffix sequences of STR markers (except the first sequence,
  as it is used as the reference). Additionally, optional bases may be
  represented by lowercase letters.
- If no explicit prefix/suffix is given for an alias, the prefix/suffix of
  the corresponding marker is assumed instead. This situation was not
  handled correctly when converting from raw sequences to TSSV or
  allelename format, which resulted in the alias remaining unused.
- Includes [Libconvert v1.1.1](#Libconvert-111)
- Includes [Library v1.0.2](#Library-102)
- Includes [Pipeline v1.0.2](#Pipeline-102)
- Includes [Vis v1.0.3](#Vis-103)
- Includes [Samplevis v2.1.2](#Samplevis-212)
- Includes [Stuttermodelvis v2.0.2](#Stuttermodelvis-202)


### Version 1.0.0
- Fixed bug that caused variant descriptions in allele names of non-STR
  markers to be prepended with plus signs similar to suffix variants
  in STR markers; when attempting to convert these allele names back to raw
  sequences, FDSTools would crash with an 'Invalid allele name' error.
- Tools that take a list of files as their argument (through the -i option
  or as positionals) now explicitly support '*' and '?' wildcards.
- Includes [BGEstimate v1.1.1](#BGEstimate-111)
- Includes [BGMerge v1.0.2](#BGMerge-102)
- Includes [Library v1.0.1](#Library-101)
- Includes [Pipeline v1.0.1](#Pipeline-101)
- Includes [Stuttermodel v1.1.1](#Stuttermodel-111)
- Includes [Allelevis v2.0.1](#Allelevis-201)
- Includes [Samplevis v2.1.1](#Samplevis-211)
- Includes [Stuttermodelvis v2.0.1](#Stuttermodelvis-201)


### Version 0.0.5
- The Blame tool was removed in favour of BGAnalyse.
- Includes [BGAnalyse v1.0.0](#BGAnalyse-100)
- Includes [Libconvert v1.1.0](#Libconvert-110)
- Includes [Library v1.0.0](#Library-100)
- Includes [Pipeline v1.0.0](#Pipeline-100)
- Includes [Samplestats v1.1.0](#Samplestats-110)
- Includes [TSSV v1.0.2](#TSSV-102)
- Includes [Vis v1.0.2](#Vis-102)
- Includes [Allelevis v2.0.0](#Allelevis-200)
- Includes [BGAnalysevis v1.0.0](#BGAnalysevis-100)
- Includes [BGRawvis v2.0.0](#BGRawvis-200)
- Includes [Profilevis v2.0.0](#Profilevis-200)
- Includes [Samplevis v2.1.0](#Samplevis-210)
- Includes [Stuttermodelvis v2.0.0](#Stuttermodelvis-200)


### Version 0.0.4
- FDSTools will now print profiling information to stdout when the
  -d/--debug option was specified.
- Fixed bug where specifying '-' as the output filename would be taken
  literally, while it should have been interpreted as 'write to standard
  out'. (Affected tools: BGCorrect, Samplestats, Seqconvert, Stuttermark)
- Added more detailed license information to FDSTools.
- Updated bundled JavaScript library Vega to v2.6.0
- Updated bundled JavaScript library D3 to v3.5.17
- Includes [BGCorrect v1.0.1](#BGCorrect-101)
- Includes [BGEstimate v1.1.0](#BGEstimate-110)
- Includes [BGMerge v1.0.1](#BGMerge-101)
- Includes [BGPredict v1.0.1](#BGPredict-101)
- Includes [Libconvert v1.0.1](#Libconvert-101)
- Includes [Samplestats v1.0.1](#Samplestats-101)
- Includes [Seqconvert v1.0.1](#Seqconvert-101)
- Includes [Stuttermodel v1.1.0](#Stuttermodel-110)
- Includes [TSSV v1.0.1](#TSSV-101)
- Includes [Vis v1.0.1](#Vis-101)
- Includes [Allelevis v1.0.0beta2](#Allelevis-100beta2)
- Includes [BGRawvis v1.0.1](#BGRawvis-101)
- Includes [Profilevis v1.0.1](#Profilevis-101)
- Includes [Samplevis v2.0.1](#Samplevis-201)
- Includes [Stuttermodelvis v1.0.0beta2](#Stuttermodelvis-100beta2)


### Version 0.0.3
- Updated bundled JavaScript library Vega to v2.5.0
- Updated bundled JavaScript library D3 to v3.5.12
- Includes [Allelefinder v1.0.0](#Allelefinder-100)
- Includes [BGCorrect v1.0.0](#BGCorrect-100)
- Includes [BGEstimate v1.0.0](#BGEstimate-100)
- Includes [BGHomRaw v1.0.0](#BGHomRaw-100)
- Includes [BGHomStats v1.0.0](#BGHomStats-100)
- Includes [BGMerge v1.0.0](#BGMerge-100)
- Includes [BGPredict v1.0.0](#BGPredict-100)
- Includes [Blame v1.0.0](#Blame-100)
- Includes [FindNewAlleles v1.0.0](#FindNewAlleles-100)
- Includes [Libconvert v1.0.0](#Libconvert-100)
- Includes [Samplestats v1.0.0](#Samplestats-100)
- Includes [Seqconvert v1.0.0](#Seqconvert-100)
- Includes [Stuttermark v1.5.0](#Stuttermark-150)
- Includes [Stuttermodel v1.0.0](#Stuttermodel-100)
- Includes [TSSV v1.0.0](#TSSV-100)
- Includes [Vis v1.0.0](#Vis-100)
- Includes [Allelevis v1.0.0beta1](#Allelevis-100beta1)
- Includes [BGRawvis v1.0.0](#BGRawvis-100)
- Includes [Profilevis v1.0.0](#Profilevis-100)
- Includes [Samplevis v2.0.0](#Samplevis-200)
- Includes [Stuttermodelvis v1.0.0beta1](#Stuttermodelvis-100beta1)


### Version 0.0.2
- Added global -d/--debug switch.
- Includes [Stuttermark v1.4.0](#Stuttermark-140)


### Version 0.0.1
- Initial version.
- Includes [Stuttermark v1.3.0](#Stuttermark-130)



Allelefinder
------------
### Allelefinder 1.1.0
- Removed the -c/--stuttermark-column argument. Allelefinder will look for stutter
  annotations in the 'flags' column instead, if present.
- The default value of the -x/--max-noisy argument has changed from 2 to 0.1; the
  value will now be interpreted as a fraction of markers if it is below 1.
- Allelefinder will now use the total_corrected column if present, else it will use
  the total column. The forward and reverse columns are no longer used by this tool.
- The -a/--max-alleles option now defaults to 1 for markers on the mitochondrial
  genome or the Y chromosome.
- Added the -N/--min-reads-lowest argument (default: 15). Allelefinder will not call
  any alleles on a marker if the lowest allele appears to have less than this number
  of reads (analogous to the -n/--min-reads option for the highest allele).


### Allelefinder 1.0.1
- Fixed crash that occurred when converting sequences to allele name format
  when no library file was provided.
- Shut down cleanly when the output pipe is closed.


### Allelefinder 1.0.0
- Initial version.



BGAnalyse
---------
### BGAnalyse 1.1.0
- Fixed bug that caused non-integer percentiles to appear rounded down in the output.
- BGAnalyse will now use the total_corrected column instead of the forward_corrected
  and reverse_corrected columns. As a result, samples with 0 reads on one strand for
  a genuine allele are now accepted.


### BGAnalyse 1.0.1
- Shut down cleanly when the output pipe is closed.


### BGAnalyse 1.0.0
- Initial version.



BGCorrect
---------
### BGCorrect 1.1.0
- Added the -C/--combine-strands option, to apply noise correction only for the total
  number of reads instead of separately for either strand.
- The correction_flags column may now contain a comma-separated list of 'corrected_X'
  flags, one for each tool used in the noise profile of that particular allele.
- Greatly reduced memory usage.


### BGCorrect 1.0.2
- Don't crash on empty input files.
- Shut down cleanly when the output pipe is closed.


### BGCorrect 1.0.1
- Added new column 'weight' to the output. The value in this column
  expresses the number of times that the noise profile of that allele
  fitted in the sample.


### BGCorrect 1.0.0
- Initial version.



BGEstimate
----------
### BGEstimate 1.2.0
- Added the -C/--combine-strands option, to estimate noise profiles for the total number
  of reads instead of estimating separate noise profiles for either strand.
- Removed the random subsampling arguments.
- The 'tool' output column has been renamed to 'tools'.


### BGEstimate 1.1.2
- Shut down cleanly when the output pipe is closed.


### BGEstimate 1.1.1
- Added option -p/--profiles which can be used to provide a previously
  created background noise profiles file, from which starting values will
  be read instead of assuming zero noise.


### BGEstimate 1.1.0
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


### BGEstimate 1.0.0
- Initial version.



BGHomRaw
--------
### BGHomRaw 1.1.0
- Added the -C/--combine-strands option, to calculate noise ratios for the total number
  of reads only, instead of calculating separate noise ratios for either strand.


### BGHomRaw 1.0.1
- Clarified the 'Allele x of marker y has 0 reads' error message with the
  name of the sample that triggered the error.
- Shut down cleanly when the output pipe is closed.


### BGHomRaw 1.0.0
- Initial version.



BGHomStats
----------
### BGHomStats 1.1.0
- Added the -C/--combine-strands option, to calculate statistics for the total number
  of reads only, instead of calculating separate statistics for either strand.
- Removed the random subsampling arguments.
- The 'tool' output column has been renamed to 'tools'.


### BGHomStats 1.0.1
- Error messages about the input data now contain the name of the sample
  that triggered the error.
- Shut down cleanly when the output pipe is closed.


### BGHomStats 1.0.0
- Initial version.



BGMerge
-------
### BGMerge 1.1.0
- Noise profiles for forward, reverse and total noise are now merged independently.
- The 'tool' output column has been renamed to 'tools' and now contains all tools used
  for the merged noise profile.


### BGMerge 1.0.3
- Shut down cleanly when the output pipe is closed.


### BGMerge 1.0.2
- Minor changes to facilitate explicit filename wildcard support.


### BGMerge 1.0.1
- Reduced memory usage.


### BGMerge 1.0.0
- Initial version.



BGPredict
---------
### BGPredict 1.1.0
- The default value of the -t/--min-r2 option has been changed to 0, effectively
  disabling the filter by default.
- Added the -C/--combine-strands option, to use stutter models trained on the total
  number of reads instead of the separate models per strand.
- Flanks are now ignored even if the library file is given; repeats may
  be interpreted as being slightly shorter if they continue into the flanks.
- The 'tool' output column has been renamed to 'tools'.


### BGPredict 1.0.2
- Don't crash on empty input files.
- Shut down cleanly when the output pipe is closed.


### BGPredict 1.0.1
- Greatly reduced memory usage.
- BGPredict will now output nonzero values below the threshold set by
  -n/--min-pct if the predicted noise ratio of the same stutter on the
  other strand is above the threshold. Previously, values below the
  threshold were clipped to zero, which may cause unnecessarily high strand
  bias in the predicted profile.


### BGPredict 1.0.0
- Initial version.



FindNewAlleles
--------------
### FindNewAlleles 1.1.1
- Added the -r/--remove-allele-flags option. When specified, the 'allele' flag is
  removed from those alleles that are flagged as 'novel' by this tool.
- When a marker is not present in the file with known alleles, no alleles will be
  marked 'novel' for that marker anymore.


### FindNewAlleles 1.1.0
- New alleles are now flagged as 'novel' in the 'flags' column.
- Removed the -m/--marker option.


### FindNewAlleles 1.0.1
- Don't crash on empty input files.
- Shut down cleanly when the output pipe is closed.


### FindNewAlleles 1.0.0
- Initial version.



Libconvert
----------
### Libconvert 1.2.0
- Removed conversion from FDSTools to TSSV format.
- Removed the -a/--aliases option.


### Libconvert 1.1.2
- Shut down cleanly when the output pipe is closed.


### Libconvert 1.1.1
- Adjustments for supporting IUPAC notation in prefix and suffix sequences
  when converting from FDSTools to TSSV library format.


### Libconvert 1.1.0
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


### Libconvert 1.0.1
- Specifying '-' as the first positional argument to libconvert will now
  correctly interpret this as "read from stdin" instead of throwing a "file
  not found" error (or reading from a file named "-" if it exists).


### Libconvert 1.0.0
- Initial version.



Library
-------
### Library 1.1.1
- Added the -m/--microhaplotypes option, to include the new [microhaplotype_positions]
  section in the library file.


### Library 1.1.0
- Introduce the new 'smart' library file and make it the default. All explanatory comments
  have been updated.
- Added the -b/--builtin option to pre-fill the new library file with built-in info.
- Removed the -a/--aliases option.


### Library 1.0.3
- Shut down cleanly when the output pipe is closed.


### Library 1.0.2
- Added documentation for IUPAC support to the descriptive comment of the
  [prefix] section.


### Library 1.0.1
- Updated some of the comments describing the sections.
- Added proper examples for non-STR markers and aliases.


### Library 1.0.0
- Initial version.



Pipeline
--------
### Pipeline 1.1.1
- Fixed issues with handling of TSSV's aggregate-filtered option.
- Fixed issues with handling of the max-alleles option of Samplestats.
- Removed access to the uncall-alleles option of Samplestats.


### Pipeline 1.1.0
- Added the -C/--combine-strands option, to analyse noise for the total number of
  reads instead of separately for either strand.
- When the -d/--debug option is used, it is now applied to all tools ran in the pipeline.
- Fixed bug that required specifying the -r/--store-predictions option to run a pipeline
  file that contained store-predictions=True.


### Pipeline 1.0.4
- Removed reference to the 'is-fastq' option of TSSV.


### Pipeline 1.0.3
- Fixed glitch that caused the 'bgprofiles.html' output file of the
  reference-database analysis to lack a proper title.


### Pipeline 1.0.2
- Added -A/--in-allelelist option, with which an existing allele list file
  can be provided when running the reference-database analysis pipeline,
  bypassing Allelefinder.


### Pipeline 1.0.1
- Removed checking of the existence of the files specified for the
  -S/--in-samples option; instead, this is left to the downstream tools to
  find out, consistent with how this works with other input file options.
- Only output the running commands if the -d/--debug option was specified.


### Pipeline 1.0.0
- Initial version.



Samplestats
-----------
### Samplestats 1.3.1
- Fixed glitch that caused addition of an empty flag in the flags column if it was
  present but empty in the input.
- Added -U/--uncall-alleles option, which lets Samplestats remove allele calls from
  sequences that do not meet the provided criteria.


### Samplestats 1.3.0
- Added the -E/--min-allele-reads and -D/--max-nonallele-pct options, which can be used
  to suppress allele calls for low-coverage and high-noise markers, respectively.
- Added the -G/--max-alleles option, which can be used to suppress allele calls if more
  alleles than expected are called at a marker.
- Added the -l/--library and -F/--sequence-format options.
- Changed the default value of the -b/--min-per-strand option from 1 to 0.
- Changed the default value of the -B/--min-per-strand-filt option from 1 to 0.


### Samplestats 1.2.0
- Fixed bug where the 'Other sequences' could be treated as the maximum
  against which the '*_mp_max' columns are calculated.


### Samplestats 1.1.1
- Don't crash on empty input files.
- Shut down cleanly when the output pipe is closed.


### Samplestats 1.1.0
- Changed default allele calling option thresholds:
    - Changed default value of -m/--min-pct-of-max from 5.0 to 2.0.
    - Changed default value of -p/--min-pct-of-sum from 3.0 to 1.5.
- Mentioned allele calling in the tool descriptions.


### Samplestats 1.0.1
- Samplestats will now round to 4 or 5 significant digits if a value is
  above 1000 or 10000, respectively. Previously, this was only done for the
  combined 'Other sequences' values.
- The 'Other sequences' lines will now also include values for
  total_recovery, forward_recovery, and reverse_recovery.
- The total_recovery, forward_recovery, and reverse_recovery columns are no
  longer placed to the left of all the other columns generated by
  Samplestats.
- The help text for Samplestats erroneously listed the X_recovery_pct
  instead of X_recovery.
- Added support for the new 'weight' column produced by BGCorrect when the
  -a/--filter-action option is set to 'combine'.


### Samplestats 1.0.0
- Initial version.



Seqconvert
----------
### Seqconvert 1.1.0
- STRNaming is now used as the back-end for generating allele names.


### Seqconvert 1.0.2
- Shut down cleanly when the output pipe is closed.


### Seqconvert 1.0.1
- Internal naming of the first positional argument was changed from
  'format' to 'sequence-format'. This was done for consistency with the
  -F/--sequence-format option in other tools, giving it the same name in
  Pipeline configuration files.


### Seqconvert 1.0.0
- Initial version.



Stuttermark
-----------
### Stuttermark 1.6.1
- Fixed glitch that caused addition of an empty flag in the flags column if it was
  present but empty in the input.


### Stuttermark 1.6.0
- Stuttermark now writes its annotations to the 'flags' column.
- Removed the -c/--column-name option to change the output column name.
- Removed the 'ALLELE' and 'UNKNOWN' annotations.


### Stuttermark 1.5.1
- Don't crash on empty input files.
- Shut down cleanly when the output pipe is closed.


### Stuttermark 1.5.0
- Changed column names 'name' and 'allele' to 'marker' and 'sequence',
  respectively. WARNING: Stuttermark is now INCOMPATIBLE with output
  from TSSV, but made compatible with TSSV-Lite and the new, bundled TSSV
  tool instead.


### Stuttermark 1.4.0
- Stuttermark now accepts raw sequences and allele names as input, which
  are automatically rewritten as TSSV-style sequences using a specified
  library file.
- The 'name' column is now optional.


### Stuttermark 1.3.0
- First version of Stuttermark to be included in `fdstools`.
- Fixed crash that occurred when an empty allele (e.g., a primer dimer)
  was encountered.
- Stuttermark now prints a warning if an allele is encountered that is
  not a TSSV-style sequence.


### Stuttermark 1.2.0
- All settings are now available from the command line.
- Use 1-based indexing in 'STUTTER' annotations.


### Stuttermark 1.1.0
- Stuttermark now accepts file names and the minimum number of reads to
  evaluate as command line arguments.


### Stuttermark 1.0.0
- Initial version.



Stuttermodel
------------
### Stuttermodel 1.2.0
- Changed the default value for -t/--min-r2 from 0.8 to 0.5. Stutter correction
  based on a model with an r-squared score of 0.5 is usually better than no correction.
- Added the -C/--combine-strands option, to model stutter for the total number of
  reads instead of fitting a separate stutter prediction model for either strand.
- Added the -T/--num-threads option (default: 1), which controls the number worker
  threads to use when scanning for potential stutter events.
- Removed the random subsampling arguments.
- Flanks are now ignored even if the library file is given; repeats may
  be interpreted as being slightly shorter if they continue into the flanks.
- Duplicate models are no longer produced for palindromic repeat units such as AT.


### Stuttermodel 1.1.2
- Shut down cleanly when the output pipe is closed.


### Stuttermodel 1.1.1
- Minor change to internal variant representation.


### Stuttermodel 1.1.0
- Stuttermodel will now only output a fit for one strand if it could also
  obtain a fit for the other strand (for the same marker, unit, and stutter
  depth). This new behaviour can be disabled with a new -O/--orphans
  option.
- Fixed bug that caused Stuttermodel to output only the raw data points for
  -1 and +1 stutter when normal output was supressed.


### Stuttermodel 1.0.0
- Initial version.



TSSV
----
### TSSV 2.1.1
- Fixed "error: too many values to unpack (expected 2)" when using the -D/--dir option.
  Thanks to @agynna for reporting!


### TSSV 2.1.0
- Changed the default value for -a/--minimum from 1 to 2.
- Replaced the -A/--aggregate-filtered option with the -B/--no-aggregate-filtered option.
- TSSV will no longer assign one read to multiple, seemingly overlapping markers.
  This makes it possible to configure copies of multi-copy markers as separate markers
  and have TSSV figure out from which copy a read originates. It also solves situations
  like the DYS389I/II ambiguity.
- Added the -L/--flank-length option, to specify the length of the anchor sequences to
  use, unless explicitly specified per marker in the library file.
- The -m/--mismatches option value will now be taken as an absolute number of allowed
  mismatches if it is larger than 1.
- The flanks may now contain IUPAC ambiguity codes, which can be used to match against
  degenerate bases in the primers or bisulfite-converted Cytosine bases in methylation-based
  essays (C/T=Y).
- Added support for GZipped FastA/FastQ files.
- Removed the error that occurred when the -D/--dir option is used and
  the output directory already exists.
- An informative error message is now displayed when an empty library file is given.
- Don't abort if the report stream is closed; log a message to the
  report stream if the main output stream is closed.


### TSSV 2.0.0
- Removed dependency on external `tssv` package (it is no longer compatible).
- Greatly increased performance by deduplicating the input reads.
- Removed the -q/--is-fastq option in favour of automatic detection.
- Changed the default value for -m/--mismatches from 0.08 to 0.1.
- Changed the default value for -n/--indel-score from 1 to 2.
- Added the -X/--no-deduplicate option to disable deduplication.
- The -D/--dir option can now be used together with -T/--num-threads.
- Fixed potential crash that could occur under very specific circumstances.


### TSSV 1.1.1
- Fixed incorrect calculation of tLeft, fLeft, rLeft, tRight and fRight
  columns in the report output file when -T/--num-threads was set to 2 or
  higher. The primary output was unaffected.


### TSSV 1.1.0
- Added option '-T/--num-threads' (default: 1), which controls the number
  of worker threads TSSV may spawn to run the analysis in parallel.
- Shut down cleanly when the output pipe is closed.


### TSSV 1.0.2
- Added new option -n/--indel-score which can be used to increase the
  penalty given to insertions and deletions in the flanking sequences
  w.r.t. the penalty given to mismatches.
- NOTE: Requires TSSV v0.4.0 or newer to be installed.


### TSSV 1.0.1
- Renamed the '--is_fastq' option to '--is-fastq', which was the only
  option with an underscore instead of a hyphen in FDSTools.
- Fixed crash that would occur if -F/--sequence-format was set to anything
  other than 'raw'.


### TSSV 1.0.0
- Initial version.



Vis
---
### Vis 1.1.0
- Changed the default value of the -B/--bias-threshold option from 25 to 0.
- Changed the default value of the -Z/--allele-min-per-strand option from 1 to 0.
- Removed the -O/--online option.


### Vis 1.0.4
- Shut down cleanly when the output pipe is closed.


### Vis 1.0.3
- The -n/--min-abs and -s/--min-per-strand options now accept non-integer
  values as well.
- Added six options to control the Table Filtering Options of Samplevis.
- Grouped some options as 'Display Options' in the command line help.


### Vis 1.0.2
- Changed default value of -n/--min-abs from 15 to 5.
- Added -I/--input2 option, which allows for specifying a file with raw
  data points for Stuttermodelvis and Profilevis.
- Added support for creating BGAnalysevis visualisations.


### Vis 1.0.1
- Added -j/--jitter option for Stuttermodelvis (default: 0.25).
- Fixed bug where Vis would not allow the -n/--min-abs and the
  -s/--min-per-strand options to be set to 0.


### Vis 1.0.0
- Initial version.



Allelevis
---------
### Allelevis 2.0.1
- Added tooltip support to HTML visualisations.


### Allelevis 2.0.0
- Replaced the simple Options overlay with responsive design options panels
  in HTML visualisations.
- Reduced Vega graph spec complexity by using the new Rank transform to
  position the subgraphs.
- Fixed glitch that caused unnecessary padding around the graph.


### Allelevis 1.0.0beta2
- Fixed potential crash/corruption that could occur with very unfortunate
  combinations of sample names and marker names.
- HTML visualisations made with the -O/--online option of the Vis tool will
  now contain https URLs instead of http.
- Added two more colours to the legend, such that a maximum of 22 markers
  is now supported without re-using colours.


### Allelevis 1.0.0beta1
- Initial version.



BGAnalysevis
------------
### BGAnalysevis 1.0.0
- Initial version.



BGRawvis
--------
### BGRawvis 2.0.1
- Changed default save filename in HTML visualisations to 'bgprofiles-raw'.
- Fixed glitch where, in HTML visualisations with embedded data and a
  custom title, the custom title was truncated to the last '.' as if it
  were a file name.


### BGRawvis 2.0.0
- Replaced the simple Options overlay with responsive design options panels
  in HTML visualisations.
- Sequences are now sorted by CE allele length when applicable.
- Changed default minimum number of reads from 15 to 5.
- Added marker selection menu for easier filtering.


### BGRawvis 1.0.1
- Fixed a JavaScript crash that would occur in HTML visualisations if the
  Marker name filter resulted in an invalid regular expression (e.g., when
  the entered value ends with a backslash).
- Reduced Vega graph spec complexity by using the new Rank transform to
  position the subgraphs.
- HTML visualisations made with the -O/--online option of the Vis tool will
  now contain https URLs instead of http.


### BGRawvis 1.0.0
- Initial version.



Profilevis
----------
### Profilevis 2.0.2
- Added ability to display background noise profiles operating on total read counts.


### Profilevis 2.0.1
- Changed default save filename in HTML visualisations to 'bgprofiles'.
- Fixed glitch where, in HTML visualisations with embedded data and a
  custom title, the custom title was truncated to the last '.' as if it
  were a file name.


### Profilevis 2.0.0
- Replaced the simple Options overlay with responsive design options panels
  in HTML visualisations.
- Alleles and sequences are now sorted by CE allele length when applicable.
- Added option to plot BGHomRaw data on top of the profiles.
- Added marker selection menu for easier filtering.


### Profilevis 1.0.1
- Fixed a JavaScript crash that would occur in HTML visualisations if the
  Marker name filter resulted in an invalid regular expression (e.g., when
  the entered value ends with a backslash).
- Reduced Vega graph spec complexity by using the new Rank transform to
  position the subgraphs.
- HTML visualisations made with the -O/--online option of the Vis tool will
  now contain https URLs instead of http.


### Profilevis 1.0.0
- Initial version.



Samplevis
---------
### Samplevis 2.3.0
- Fixed an issue that prevented manually marking additional alleles.
- Fixed an issue with restoring allele calls after the visualisation has been saved.
- Added ability to display background noise profiles operating on total read counts.
- Added ability to display multiple noise correction notes for a single sequence.
- Added ability to display and save notes based on the optional 'flags' column, such
  as Novel allele calls from FindNewAlleles.
- Changed default minimum number of reads per strand for table filtering from 1 to 0.
- Changed the default bias threshold from 25% to 0% (disabled).


### Samplevis 2.2.2
- Minor change to the calculation of percentage of forward reads to prevent
  roundoff effects.


### Samplevis 2.2.1
- Added an option to apply graph filtering before noise correction (on by
  default).


### Samplevis 2.2.0
- Fixed incorrect calculation of 'percentage of highest' if the 'sequence'
  with the highest read count within a marker is the aggregated 'Other
  sequences' data. In exceptional cases, this could have resulted in the
  erroneous omission of an allele in the visualisation (graphs and/or
  tables).


### Samplevis 2.1.2
- Added 'Save page' link to HTML visualisations, which offers for download
  a copy of the entire HTML visualisation including the user's changes.
- Added automatic allele calling to static visualisations.
- The net effect of the allele calling thresholds (table filtering options)
  is now visualised in the graphs as a dashed vertical red line.


### Samplevis 2.1.1
- Added tooltip support to HTML visualisations.
- The tooltip may include a 'new allele' note if the input sample was
  analysed with FindNewAlleles.
- The allele tables in HTML visualisations will now grow much wider than
  before if the screen (or window) is very narrow.
- Improved line breaking behaviour in the tables in HTML visualisations.
- Improved determination of column widths of the allele tables when
  printing an HTML visualisation.
- When printing an HTML visualisation, the graph and the corresponding
  table of a marker will be kept on the same page in all browsers now.
- Fixed glitch that caused 'Infinity%' or 'NaN%' to be written in some
  cells in the allele tables in HTML visualisations.


### Samplevis 2.1.0
- Changed default minimum number of reads for graph filtering from 15 to 5.
- Changed default table filtering options:
    - Percentage of highest allele per marker changed from 5% to 2%.
    - Percentage of the marker's total reads changed from 3% to 1.5%.
    - Minimum number of reads in both orientations changed from 0 to 1.


### Samplevis 2.0.1
- Fixed a JavaScript crash that would occur in HTML visualisations if the
  Marker name filter resulted in an invalid regular expression (e.g., when
  the entered value ends with a backslash).
- Reduced Vega graph spec complexity by using the new Rank transform to
  position the subgraphs.
- Fixed a glitch in HTML visualisations where clicking the 'Truncate
  sequences to' label would select the marker spacing input.
- In HTML visualisations, the 'Notes' table cells with 'BGPredict' in them
  now get a light orange background to warn the user that their background
  profile was computed. If a sequence was explicitly 'not corrected', 'not
  in ref db', or 'corrected as background only', the same colour is used.
- The message bar at the bottom of Samplevis HTML visualisations will now
  grow no larger than 3 lines. A scroll bar will appear as needed.
- HTML visualisations made with the -O/--online option of the Vis tool will
  now contain https URLs instead of http.


### Samplevis 2.0.0
- Initial version.



Stuttermodelvis
---------------
### Stuttermodelvis 2.0.4
- Added ability to display stutter models fitted on total read counts.


### Stuttermodelvis 2.0.3
- Fixed bug that caused HTML visualisations with embedded data to fail
  while loading.
- Fixed glitch where, in HTML visualisations with embedded data and a
  custom title, the custom title was truncated to the last '.' as if it
  were a file name.


### Stuttermodelvis 2.0.2
- Added filtering option for the stutter amount (-1, +1, -2, etc.).
- Added filtering option for the coefficient of determination (r squared
  value) of the fit functions.


### Stuttermodelvis 2.0.1
- Changed the unit in the horizontal axis title from 'bp' to 'nt'.


### Stuttermodelvis 2.0.0
- Replaced the simple Options overlay with responsive design options panels
  in HTML visualisations.
- Fixed glitch that caused the graphs to be re-rendered twice when loading
  a file by drag-and-drop in HTML visualisations.
- Fixed glitch that made it possible to replace the data that was embedded
  in an HTML visualisation through drag-and-drop.
- Added repeat unit selection menu for easier filtering.


### Stuttermodelvis 1.0.0beta2
- HTML visualisations now support drawing raw data points on top of the fit
  functions. The points can be drawn with an adjustable jitter to reduce
  overlap.
- Fixed a JavaScript crash that would occur in HTML visualisations if the
  Repeat unit or Marker name filter resulted in an invalid regular
  expression (e.g., when the entered value ends with a backslash).
- Reduced Vega graph spec complexity by using the new Rank transform to
  position the subgraphs.
- HTML visualisations made with the -O/--online option of the Vis tool will
  now contain https URLs instead of http.


### Stuttermodelvis 1.0.0beta1
- Initial version
