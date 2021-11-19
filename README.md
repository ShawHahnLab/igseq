# IgSeq Utilities

**Work in Progress**

A command-line tool for various common Ig-Seq tasks.  These are heavily biased
toward the peculiarities of our protocol and for rhesus macaque antibody
sequences.  Your mileage may vary.

## Install

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

Install from the latest version via <https://anaconda.org/ShawHahnLab/igseq>:

    conda create --name igseq -c conda-forge -c bioconda -c ShawHahnLab igseq

Or, install from the latest source here:

    git clone https://github.com/ShawHahnLab/igseq.git
    cd igseq
    conda env update --file igseq/data/environment.yml
    conda activate igseq
    pip install .

## Some Instructions

The `igseq` command is organized into subcommands, grouped into two
categories: early read processing tasks (demultiplex, trim, merge, etc.), and
various convenience tools (IgBLAST this against that, what database has the
closest V gene, etc.).

Read processing subcommands:

 * getreads: Run [bcl2fastq] with an autogenerated sample sheet and some
   customized settings to write Undetermined I1/R1/R2 fastq.gz files.
 * demux: Demultiplex the I1/R1/R2 files according to per-sample barcodes
 * phix: Map reads left unassigned post-demux to the PhiX genome for
   troubleshooting
 * trim: Run [cutadapt] to remove adapter/primer/barcode and low-quality
   sequences on a per-sample basis
 * merge: Merge R1/R2 for each sample with [PEAR].

Each step in the read processing produces a read counts summary CSV table
`<step>.counts.csv` and has default output paths derived from the inputs, and
most can work per-sample or per-directory, so it's easy to chain together the
steps for a given run:

    igseq getreads /seq/runs/211105_M05588_0469_000000000-JWV49
    igseq demux -s samples.csv analysis/reads/211105_M05588_0469_000000000-JWV49
    igseq phix analysis/demux/211105_M05588_0469_000000000-JWV49
    igseq trim -s samples.csv -S rhesus analysis/demux/211105_M05588_0469_000000000-JWV49
    igseq merge analysis/trim/211105_M05588_0469_000000000-JWV49

The `samples.csv` file is a table matching sample names to run IDs, barcode
IDs, and antibody chain types, like:

    Sample,Run,BarcodeFwd,BarcodeRev,Type
    wk12H,211105_M05588_0469_000000000-JWV49,1,1,gamma
    wk12K,211105_M05588_0469_000000000-JWV49,2,2,kappa
    wk24H,211105_M05588_0469_000000000-JWV49,3,3,gamma
    wk24K,211105_M05588_0469_000000000-JWV49,4,4,kappa

The barcode IDs refer to the numbered barcodes for the protocol, with the
varying-length randomized prefix for the forward barcodes:

    $ igseq show barcodes
     Direction BC              Seq
             F  1     NNNNAACCACTA
             F  2    NNNNNAACTCTAA
             F  3   NNNNNNAAGGCCCT
             F  4  NNNNNNNAATATGTC
             F  5 NNNNNNNNAATCGTCA
    ...
             R  1         TAGTGGTT
             R  2         TTAGAGTT
             R  3         AGGGCCTT
             R  4         GACATATT
             R  5         TGACGATT
    ...

The chain type is used to select the appropriate constant region primer:

    $ igseq show primers
     Species    Type                               Seq
       human   gamma  GCCAGGGGGAAGACCGATGGGCCCTTGGTGGA
       human   alpha GAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGG
       human      mu AGGAGACGAGGGGGAAAAGGGTTGGGGCGGATG
       human epsilon GCGGGTCAAGGGGAAGACGGATGGGCTCTGTGT
       human   delta CTGATATGATGGGGAACACATCCGGAGCCTTGG
       human   kappa GCGGGAAGATGAAGACAGATGGTGCAGCCACAG
       human  lambda GGCCTTGTTGGCTTGAAGCTCCTCAGAGGAGGG
      rhesus   gamma  GCCAGGGGGAAGACCGATGGGCCCTTGGTGGA
      rhesus   alpha GAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGG
      rhesus      mu GAGACGAGGGGGAAAAGGGTTGGGGCGGATGCA
      rhesus epsilon CGGGTCAAGGGGAAGACGGATGGGCTCTGTGTG
      rhesus   delta CTGATATGATGGGGAACACATCCGGAGCCTTGG
      rhesus   kappa GCGGGAAGATGAAGACAGATGGTGCAGCCACAG
      rhesus  lambda GGCCTTGTTGGCTTGAAGCTCCTCAGAGGAGGG

See `igseq/data/examples/readproc.sh` for an example read processing workflow
from start to finish with a small set of reads.

[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
[cutadapt]: https://cutadapt.readthedocs.io
[PEAR]: https://cme.h-its.org/exelixis/web/software/pear/
