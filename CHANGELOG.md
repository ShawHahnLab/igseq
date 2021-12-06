# Changelog

## 0.1.1 - 2021-12-06

### Changed

 * Show a warning for zero file matches in list and show commands ([#8])

### Fixed

 * Header-only CSV/TSV files no longer crash the show command ([#8])
 * summarize command now works with multiple references ([#7])
 * vdj-match and summarize commands now work as intended via igblast, and all
   have basic automated tests ([#5])
 * igblastn arguments can now be given with a two-dash prefix to ensure they
   aren't interpreted as igseq arguments ([#4])
 * Duplicate FASTA paths found in vdj-gather will no longer result in
   duplicated output sequences ([#2])

[#8]: https://github.com/ShawHahnLab/igseq/pull/8
[#7]: https://github.com/ShawHahnLab/igseq/pull/7
[#5]: https://github.com/ShawHahnLab/igseq/pull/5
[#4]: https://github.com/ShawHahnLab/igseq/pull/4
[#2]: https://github.com/ShawHahnLab/igseq/pull/2

## 0.1.0 - 2021-11-19

First beta release
