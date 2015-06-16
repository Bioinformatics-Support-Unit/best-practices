# Exon Level Quantification

## Introduction

After alignment has been done, a set of sam<sup>1</sup> (sequence alignment map) or bam (binary sam) files, that contain the sequence alignment data, are available for further processing. Usually sam files are very large and it is more efficient to convert the sam files to bam which is a binary format and much for space efficient.


## Recommended tool

The tool reccommended for quanitifcation are two Python scripts that are provided as part of DEXSeq [1]. A requirement for running the scripts is that the HTSeq [2] package be installed.


## How to use

The first of the two Python scripts, dexseq_prepare_annotation.py, is for preparing annotation. The file requires a GTF format file with gene models of the species that is appropriate to the analysis. GTF files can be downloaded from different providers but the authors of DEXSeq suggests using Ensembl.

```python /path/to/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py
Drosophila_melanogaster.BDGP5.72.gtf Dmel.BDGP5.25.62.DEXSeq.chr.gff```

The actual counting of reads is done with the dexseq_count.py script.

```python /path/to/library/DEXSeq/python_scripts/dexseq_count.py
Dmel.BDGP5.25.62.DEXSeq.chr.gff untreated1.sam untreated1fb.txt```

## References
[1]Anders, Simon, Alejandro Reyes, and Wolfgang Huber. *Detecting differential usage of exons from RNA-seq data.* Genome research 22.10 (2012): 2008-2017.

[2]Simon Anders, Paul Theodor Pyl, Wolfgang Huber. HTSeq — *A Python framework to work with high-throughput sequencing data.* Bioinformatics (2014), in print, online at doi:10.1093/bioinformatics/btu638

<a name="1">1</a>: https://samtools.github.io/hts-specs/SAMv1.pdf
