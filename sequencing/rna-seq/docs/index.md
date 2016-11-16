# Recommended Best Practices for RNA-Seq analysis

## Motivation
RNA Sequencing is a relatively new technology for assaying the RNA content of samples. In the last few years it has gradually replaced microarrays as the technology of choice in this area. During this time, the bioinformatics protocols for analysing RNA Seq data have evolved rapidly. As such, keeping up with the 'state of the art' in transcriptomics analysis presents a significant challenge -- especially for bench-based researchers for whom bioinformatics is an occasional engagement.

The purpose of these documents is to set out an evidence-based (where possible) set of recommendations for how to best analyse an RNA Seq experiment from a variety of perspectives.

  * [Quality Control](QC) - good QC is universal to all RNA Seq experiments, whatever your downstream analysis considerations. This section contains recommendations for QC assessment and remedial action.
  * [Gene Level Analysis](gene-level) - probably the most common way of considering an RNA Seq analysis. Most closely mirrors the analysis of a 3'-IVT expression microarray. The expression of each gene is considered _in toto_.
  * [Transcript Level Analysis](transcript-level) - a more 'biologically realistic' interpretation of RNA Seq data. The expression of each gene is considered, and analysed, in terms of its component transcripts.
  * [Exon Level Analysis](exon-level) - consideration of exon usage in RNA Seq data allows for the exploration of alternative splicing, including discovery of novel, non-annotated transcripts. 

## Things you should know
Stuff here about strandedness, library prep, etc.

## Understanding your
