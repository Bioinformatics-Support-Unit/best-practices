# Gene Level Differential Expression

## Introduction

It is assumed that you have mapped and quantified your reads and you now want to evaluate differential expression (DE) of genes. At a simple level, one could imagine considering the counts mapping to a particular gene in each sample and doing a t-test between the samples in the control and intervention groups of the experiment, rejecting the null hypothesis (i.e. accepting that the gene is differentially expressed) if the  p-value is small enough. And indeed, this idea underlies most of the packages that carry out the task. However there are several difficulties that make it inadequate by itself:

1. The distribution of read counts in an RNA-Seq experiment is non-normal. The [negative binomial distribution](http://en.wikipedia.org/wiki/Negative_binomial_distribution), which is skew with a long upper tail, is a better model.

1. There is often only a small number of biological replicate samples for each condition (group) of the experiment. 3 is typical, though more is desirable. This means that the estimate of the within-sample variance is unstable: in particular, a simple variance estimate will very often be much too small if the counts for all 3 samples are about the same, leading to false positive "DE" genes. Hence a moderated estimate of the variance should be used, which is a weighted sum of the variance estimate for a particular gene and the average variance of all other genes that have a similar expression level. This approach was first implemented (in the microarray context) by the `ebayes()` function of the `limma` package [ref]. 

1. The design of the experiment may well be more complex than assays on independent samples falling into two groups. It is deriable to have the ability to handle  more complex experimental designs, for example allowing more than two groups (ANOVA rather than t-test), measurements made repeatedly on the same subjects (paired t-test rather than simple), time course measurments etc. 

1. Because many genes will be assayed (typically about 20000 protein coding genes in a mammal), it will be necessary to make a multiple testing correction to the single-test p-values if the results are not to be affected by many false positives. It is convenient to have the ability to make this multiple testing correction within the analysis tool. 

1. It is desirable that the package should provide functions to assist with data quality control at the level both of whole samples and of particular genes. 

There several bioconductor pacakges that allow all these considerations to be taken into account: `edgeR` [ref], `limma+voom` [ref(s)] and `DESeq2` [ref]. 


## Recommended tool

We recommend `DESeq2`, which addresses all the issues raised in the introduction. In addition it was produced by the same authors as htseq-count and works seamlessly with it, but can be also used with count data in any tabular format.

[//]: # (COMMENT: and performed well in a methods comparison trial - rapaport)

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

1. Describe the experimental samples. This should be a dataframe (called `colData` for example) describing the group or groups each experimental sample belongs to. In the simplest example, experimental groups correspond to the levels of a single factor,`condition` (say), which is the name of the single column of the dataframe. The rownames() functions is used to assign the sample names:
`colData`
|    |  condition|
|--- |---|
|KO1 |KO |
|KO2 |KO |
|KO3 |KO |
|WT1 |WT |
|WT2 |WT |
In more complex experiments there will be additional columns, one for each additional factor or for the experimental source (animal), if each produced more than one sample. 
1. Read data and set up DESeqDataSet (called `ddsDS` for example). Two functions are provided for this, `DESeqDataSetFromMatrix()` and `DESeqDataSetFromHTSeqCount()`. The DESeqDataSet contains the count data, the colData dataframe and the model. If using ...`FromMatrix`, the count data will already have been read in, is column names corresponding to the rownames of the `colData` object; if ...`FromHTSeqCount` it is read at this stage by passing file names to the function. For the simple example above, the model can be specified by the argument `design= ~ condition`.

1. call DESeq2, e.g. `dds <- DSEeq(ddsDS)`
1. Obtain results, `res <- results(dds)` . In the example above, only one contast is possible, but if the `condition` factor had more than 2 levels, or if there are additional factors, result sets for different contrasts can also be obtained by giving the `contrasts` argument to `results`. 

The result set, `res`, is a data frame with column names baseMean (average expression across samples), log2FoldChange, 2lfcSE (std. error of previous), stat (test statistic), pvalue, padj (multiple test corrected p value) 


### Analysis of output: Quality Control

Some of this can be run before the call to `DESeq()`, but in most cases they are run subsequently (in the examples below, either the `res` or `dds` objects are used, or an `rld` object, which is a modified log-transform of the count data, obtained by `rld <- rlog(dds)` ). 

  1. MA-plots  (`plotMA(res)`)
  1. heatmaps of sample-sample difference matrix. The difference matrix is calculated by `dists <- dist(t(assay(rld)))` and plotted with `heatmap.2` in the `gplots()` library.
  1. PCA plots (`plotPCA(rld)`). Groups should be separated, but samples of poor quality can often be detected because they are a long way away from other samples in the same group.  The dataset should be reanalysed with these groups omitted. The basic `plotPCA` can be customised by `ggplot()` (`ggplot2` library).
  1. library sizes - `sizeFactors(dds)` . This is an indication of the total read depth in each sample, (the total number of reads mapped to all genes in the sample is another). Because it is based on a median read depth, it is less sensitive to skewing by a few genes associated with a very large number of reads. 
  1.  goodness of fit (dispersion plot)(cooks distance)
  1.  histogram of p-values : can be made by applying base graphics `hist()` function to `res$pvalues`. Should be approximately flat, hopefully with a peak near zero where genuinely DE genes are. If there are peaks at higher values of p, especially near 1, it may be an indication of a large number of genes in the sample with a low overall level of expression (see also "library sizes" above), which can be filtered out prior to the call to `DESeq()`.

### Analysis of output: Genes 
  1. Results table: There is no analogue of the `toptable` function in limma. Sort the output data frame by pvalue and cut off where pdaj < the desired FDR. 
  1. volcano plot . This is not built in, but can be obtained by making a scatter plot of `res$log2FoldChange` against `-log10(res$adj)`. 

### Caveats:

be careful of versioning if you are not using an up-to-date version of R. 
With R 3.0 (still, as of 30/6/15, the default on Ubuntu linux) and installing through biocLite you will end up with a relatively
early version, DESeq 1.2.4, which is over-optimistic in the number of genes it calls differentially expressed. Later versions (1.6.3 and after) are OK.

DESeq2 allows the analysis of datasets without replicates, although it prints a warning message. As always, having replicates should be viewed as essential in practice to getting reliable results. 

DESeq2 should only be used with raw counts as input. 

[//]: # (COMMENT: to-do introduce "dispersion"; variance stabilising transformation; )
