# RNA-Seq Quality Control Read Quality

## Introduction

Raw sequences from any sequence technology will be of variable quality, with a number of problems that may be detectable in the raw data. Some of these problems may be systematic (such as degrading quality scores with increasing read length – a notable problem with Illumina sequencers), or limited to a particular subset of reads (such as adapter contamination). QC can be performed by any program, or combination of programs that recognise these types of issue with sequencing data, and provide methods for correcting it, such as filtering out bad reads, or trimming low quality or contaminating bases.

## Recommended tools

**FastQC** (<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>) is a Java program from the Babraham Institute in Cambridge, that produces a number of graphical reports on the quality of a particular FASTQ file. These reports can be valuable in detecting the types of problem that can be exhibited by RNA-Seq data. They must, however, be interpreted with care - as what can be seen by FastQC as 'aberrant' behaviour can be entirely expected.

**MultiQC** (<http://multiqc.info>) is an excellent Python package for aggregating QC reports from a number of tools, including FastQC. As RNA-Seq experiments become increasingly large, the number of FastQC reports becomes problematic. MultiQC provides an excellent means of bringing together multiple reports in a single browser.

**RSeQC** (<http://rseqc.sourceforge.net>) is another Python package that provides a range of QC tools specifically designed for RNA-Seq data. Many of these tools are useful _post-alignment_. 

## How to use
