# Gene Level Differential Expression

## Introduction

By now it is assumed that you have mapped and quantified your reads and you want to evaluate differential expression (DE) at the gene level. At an intuitive level, one could imagine, in each sample, summing the counts mapping to a particular gene and doing a t-test between the samples in the control and intervention groups of the experiment, rejecting the null hypothesis (i.e. accepting that the gene is differentially expressed) if the  p-value is small enough. And indeed, this is the simple idea that underlies most of the packages that carry out the task. However there are several difficulties that make it inadequate by itself:

1. The distribution of read counts in an RNA-Seq experiment is non-normal. The [negative binomial distribution](http://en.wikipedia.org/wiki/Negative_binomial_distribution), which has a heavy upper tail, is a better model.

1. There is often only a small number of biological replicate samples for each condition (group) of the experiment. 3 is typical, though more is desirable. This means that the estimate of the within-sample variance is unstable: in particular, a simple variance estimate will very often be much too small if the counts for all 3 samples are about the same, leading to false positive "DE" genes. Hence a moderated estimate of the variance should be used, which is a weighted sum of the variance estimate for a particular gene and the average variance of all other genes that have a similar expression level. This approach was first implemented (in the microarray context) by the `ebayes()` function of the `limma` package [ref]. 

1. The design of the experiment may well be more complex than assays on independent samples falling into two groups. It is deriable to have the ability to handle  more complex experimental designs, for example allowing more than two groups (ANOVA rather than t-test), measurements made repeatedly on the same subjects (paired t-test rather than simple), time course measurments etc. 

1. Because many genes will be assayed (typically about 20000 protein coding genes in a mammal), it will be necessary to make a multiple testing correction to the single-test p-values if the results are not to be affected by many false positives. It is convenient to have the ability to make this multiple testing correction within the analysis tool. 

There several bioconductor pacakges that allow all these considerations to be taken into account: `edgeR` [ref], `limma+voom` [ref(s)] and `DESeq2` [ref]. 


## Recommended tool

Use the bioconductor package DESeq2. We recommend it because it addresses the issues raised in the introduction and performed well in a methods comparison trial [rapaport]. In addition it was produced by the same authors as htseq-count and works seamlessly with it.

### Installation

```
source("http://bioconductor.org/biocLite.R")  
biocLite("DESeq2")
```

### Reference 

[Love MI, Huber W and Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, pp. 550.](http://dx.doi.org/10.1186/s13059-014-0550-8)

A long [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) and 
[reference manual](http://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)
 are available from the bioconductor web site.

## How to use

### simple usage:

1. Read data 
1. set up model 
1. (define contrasts if necessary) 
1. call DESeq2
```
dds<-DESeq(dds)
out<-results(dds)
(how to get the top n?). The column names of the output data frame are [TODO]


```

### Analysis of output: Quality Control (some of this can be before: you dont' need to run DESeq2 itself)  
  1. MA-plots
  1. heatmaps and sample-sample differene matrix
  1. PCA plots 

but some has to come after
  1.  goodness of fit (dispersion plot)(cooks distance) 
  1.  histogram of p-values 

### Analysis of output: Genes 
  1. volcano plot . This is not built in. 

### Caveats:

be careful of versioning if you are not using an up-to-date version of R. 
With R 3.0 and installing through biocLite you will end up with DESeq 1.x.x
which is over-optimistic in calling DE. Later versions (1.y.y.) are OK. 

DESeq2 allows the analysis of datasets without replicates, although it prints a warning message. As always, having replicates should be viewed as in practice
essential to getting reliable results. Estimation of dispersions from the different groups themselves will be poor. 

DESeq2 should only be used with raw counts as input. 


### to-do

introduce "dispersion"; variance stabilising transformation; 
