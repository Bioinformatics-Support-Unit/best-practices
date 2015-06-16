# Transcript Level Quantification

## Introduction

Much current RNA-Seq analysis essentially attempts to replicate the function of
3' expression microarrays such as AffyMetrix GeneChips or Illumina BeadArrays. Distilling the complexity of gene expression down to single observations for each gene, and summary statistics to determine which of these genes is differentially expressed. The primary problem with this kind of analysis is that it is not immediately reflective of biology. It is easy to envisage the case where a gene, which is not differentially expressed, is composed of one or more transcripts which are (see figure 1 for an example of this phenomenon from a real data set). Other confounding factors can occur, and indeed do so frequently.

![Figure 1][fig1]

The complexity of the information available in even a small RNA-Seq dataset allows us to go far beyond this standard analysis. Examining individual exons for expression and dysregulation is one way of working with this complexity - and there are [recommended tools and workflows for this type of analysis][exon] on this site. However, since most RNA-Seq experiments are designed to have fragments larger than the average exon,* we have the flexibility to do more than simply quantify exons. It is also possible to determine how the exons are joined to form cogent transcripts, and these transcripts can then be assessed for *their* differential expression. This is a biologically more meaningful analysis scenario.

It should be noted that despite the feasibility of this analysis approach with current RNA-Seq data, transcript-level quantification remains an extremely challenging problem, and one which many tools have attempted to address. It is only with vey recent developments that it seems like we are getting close to an optimal solution for this problem.

\* Fragment size is usually 250 bases, the average human exon is 248.08 bases, the median is 129 bases - determined from Ensembl release 79 genes, see figure 2.

![Figure 2][fig2]

**Figure2** - distribution of human exon sizes from GRCh38 (Ensembl release 79).

### Available Tools

#### Cufflinks
[Link to paper][cufflinks]
#### RSEM
[Link to paper][rsem]
#### eXpress
[Link to paper][express]
#### Sailfish
[Link to paper][sailfish]
#### Salmon
[Link to paper][salmon]
#### Kallisto
[Link to paper][kallisto]

## Recommended tool

Salmon/Kallisto

Salmon has - gene level quantification, bias correction

Kallisto has - bootstrapping, deterministic

[A note on TPM vs FPKM][tpm]

## How to use



 [tpm]: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
 [cufflinks]: http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html
 [rsem]: http://www.biomedcentral.com/1471-2105/12/323/abstract
 [express]: http://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2251.html
 [sailfish]: http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html
 [salmon]: http://thereisntone.com/
 [kallisto]: http://arxiv.org/abs/1505.02710
 [exon]: ../exon-level
 [fig1]: figure-1.png
 [fig2]: figure-2.png
