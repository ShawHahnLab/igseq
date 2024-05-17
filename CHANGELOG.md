# Changelog

## dev

### Added

 * `tree` now writes seq set groupings in NEXUS output files ([#67])

### Changed

 * `igblast` has its version pinned at 1.21.0 rather than left unspecified
   ([#70])

### Fixed

 * `igblast` now handles crashes during file input more clearly ([#75])
 * `phix` can accept a custom path for the read counts CSV file ([#74])

[#75]: https://github.com/ShawHahnLab/igseq/pull/75
[#74]: https://github.com/ShawHahnLab/igseq/pull/74
[#70]: https://github.com/ShawHahnLab/igseq/pull/70
[#67]: https://github.com/ShawHahnLab/igseq/pull/67

## 0.5.1 - 2023-03-31

### Changed

 * `convert` now handles edge cases for sequence input and tabular output by
   always including a sequence description column in the output ([#63])

### Fixed

 * `msa` will now bypass calling MUSCLE when called with just a single input
   sequence, avoiding a crash ([#62])
 * `convert` will now obey a custom sequence description column name if one is
   given with `--col-seq-desc` ([#60])
 * `getreads` command is now compatible with the latest available version of
   bcl2fastq, v2.20.0.422 ([#58])
 * `tree` command can now handle assigning a color code when exactly one
   sequence set is defined ([#57])

[#63]: https://github.com/ShawHahnLab/igseq/pull/63
[#62]: https://github.com/ShawHahnLab/igseq/pull/62
[#60]: https://github.com/ShawHahnLab/igseq/pull/60
[#58]: https://github.com/ShawHahnLab/igseq/pull/58
[#57]: https://github.com/ShawHahnLab/igseq/pull/57

## 0.5.0 - 2023-01-04

### Added

 * `summarize` command will automatically use all available references for a
   given species if a species is given but no references ([#50])
 * `tree` command for creating and formatting phylogenetic trees ([#44])
 * support for additional arguments for `getreads` command passed through to
   bcl2fastq ([#43])
 * `msa` command for building multiple sequence alignments with
   [MUSCLE](https://drive5.com/muscle5/) ([#41])

### Fixed

 * `convert` command and underlying input/output features now handles sequence
   descriptions ([#51])
 * `identity` command now uses a custom sequence ID column if one is given
   ([#49])

[#51]: https://github.com/ShawHahnLab/igseq/pull/51
[#50]: https://github.com/ShawHahnLab/igseq/pull/50
[#49]: https://github.com/ShawHahnLab/igseq/pull/49
[#44]: https://github.com/ShawHahnLab/igseq/pull/44
[#43]: https://github.com/ShawHahnLab/igseq/pull/43
[#41]: https://github.com/ShawHahnLab/igseq/pull/41

## 0.4.0 - 2022-09-17

### Added

 * Automatic usage of all available references for a given species in `igblast`
   command ([#39])
 * `identity` command for calculating pairwise identity between arbitrary
   queries and references ([#31], [#37])
 * Support for showing basic tree topology for Newick-format files in `show`
   command ([#33])

### Fixed

 * Uppercase file extensions are now supported by the `show` command ([#35])
 * broken pipes (such as from `igseq something | something else`) are now
   handled gracefully ([#30])

[#39]: https://github.com/ShawHahnLab/igseq/pull/39
[#37]: https://github.com/ShawHahnLab/igseq/pull/37
[#35]: https://github.com/ShawHahnLab/igseq/pull/35
[#33]: https://github.com/ShawHahnLab/igseq/pull/33
[#31]: https://github.com/ShawHahnLab/igseq/pull/31
[#30]: https://github.com/ShawHahnLab/igseq/pull/30

## 0.3.0 - 2022-07-14

### Added

 * Support for .afa (as FASTA) and .fq (as FASTQ) in input/output handling
   ([#26])

### Changed

 * Replaced Rhesus germline reference "bernat2021" with "kimdb" version 1.1,
   containing additions and corrections since the original publication, and
   added a CSV of germline reference information ([#28])

### Fixed

 * vdj-match now skips references that are missing gene segments but were only
   included implicitly rather than directly named ([#25])

[#28]: https://github.com/ShawHahnLab/igseq/pull/28
[#26]: https://github.com/ShawHahnLab/igseq/pull/26
[#25]: https://github.com/ShawHahnLab/igseq/pull/25

## 0.2.0 - 2022-02-15

### Added

 * Human germline FASTAs from IMGT ([#21])
 * support for FASTA/FASTQ/CSV/TSV query inputs for the igblast and related
   commands ([#18], [#19])
 * convert command for FASTA/FASTQ/CSV/TSV file conversion, in place of the
   more limited tab2seq command ([#14], [#16])
 * Rhesus germline HV and HJ allele FASTAs from
   [10.4049/jimmunol.1800342](https://doi.org/10.4049/jimmunol.1800342) ([#13])

[#21]: https://github.com/ShawHahnLab/igseq/pull/21
[#19]: https://github.com/ShawHahnLab/igseq/pull/19
[#18]: https://github.com/ShawHahnLab/igseq/pull/18
[#16]: https://github.com/ShawHahnLab/igseq/pull/16
[#14]: https://github.com/ShawHahnLab/igseq/pull/14
[#13]: https://github.com/ShawHahnLab/igseq/pull/13

## 0.1.1 - 2021-12-07

### Changed

 * Show a warning for zero file matches in list and show commands ([#8])

### Fixed

 * IgBLAST output is now shown even if it crashes ([#12])
 * Header-only CSV/TSV files no longer crash the show command ([#8])
 * summarize command now works with multiple references ([#7])
 * vdj-match and summarize commands now work as intended via igblast, and all
   have basic automated tests ([#5])
 * igblastn arguments can now be given with a two-dash prefix to ensure they
   aren't interpreted as igseq arguments ([#4])
 * Duplicate FASTA paths found in vdj-gather will no longer result in
   duplicated output sequences ([#2])

[#12]: https://github.com/ShawHahnLab/igseq/pull/12
[#8]: https://github.com/ShawHahnLab/igseq/pull/8
[#7]: https://github.com/ShawHahnLab/igseq/pull/7
[#5]: https://github.com/ShawHahnLab/igseq/pull/5
[#4]: https://github.com/ShawHahnLab/igseq/pull/4
[#2]: https://github.com/ShawHahnLab/igseq/pull/2

## 0.1.0 - 2021-11-19

First beta release
