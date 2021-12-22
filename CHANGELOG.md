# Changelog

## dev

### Added

 * convert command for FASTA/FASTQ/CSV/TSV file conversion, in place of the
   more limited tab2seq command ([#14])
 * Rhesus germline HV and HJ allele FASTAs from
   [10.4049/jimmunol.1800342](https://doi.org/10.4049/jimmunol.1800342) ([#13])

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
