#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : NOF and OA Microarray Study                                                |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Illumina HT-12v3 Microarray of NOF and OA Samples. Grouped analysis and    |
#                Regression by Age.
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Documents/Bioinformatics/Customers/Young/NOF_OA/Microarray/Microarray/v3/")

source("http://bioconductor.org/biocLite.R")
biocLite()

library(sva)
library(lumi)
library(gplots)
library(ggplot2)
library(annotate)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(limma)
library(reshape2)
library(scales)
library(devtools)
library(lumidat)
##'-----------------------------------------------------------------------------------------#



##'Read in Raw Data and attach pheno data
##'-----------------------------------------------------------------------------------------#
raw_data            <- lumiR("Raw_data.txt")
pheno               <- read.table("Pheno.txt",
                                  sep="\t",
                                  stringsAsFactors=F,
                                  header=T,
                                  row.names=1)

pheno$sample_match  <- colnames(raw_data)
rownames(pheno)     <- pheno$sample_match
pData(raw_data)     <- pheno
##'-----------------------------------------------------------------------------------------#



##'Subset Arrays by phenotype
##'-----------------------------------------------------------------------------------------#
arrays_in           <- pheno[pheno$Exclusions != "Exclude" &
                             pheno$Age != "ND" &
                             pheno$Tissue_Type == "NOF",]
raw_data_in         <- raw_data[, match(arrays_in$sample_match, pheno$sample_match)]

# raw_data_in         <- raw_data
##'-----------------------------------------------------------------------------------------#



##'Normalisation - VST and RSN
##'-----------------------------------------------------------------------------------------#
vst_data            <- lumiT(raw_data_in, method = "vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
lumi.Q              <- lumiQ(rsn_data)
##'-----------------------------------------------------------------------------------------#



##'Probe level filtering - Detection P Values
##'-----------------------------------------------------------------------------------------#
exprs_data          <- exprs(lumi.Q)
present_count       <- detectionCall(lumi.Q)
normalised_data     <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#



##'Probe level filtering - Detection P Values
##'-----------------------------------------------------------------------------------------#
probe_list          <- rownames(normalised_data)
nuIDs               <- probeID2nuID(probe_list)[, "nuID"]
symbol              <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name                <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df             <- data.frame(ID=nuIDs, probe_list, symbol, name)
##'-----------------------------------------------------------------------------------------#



##'Limma Model Design - Linear Variable - Age
##'-----------------------------------------------------------------------------------------#
design    <- model.matrix(~as.numeric(pData(raw_data_in)$Age))
fit       <- lmFit(normalised_data, design)
fit       <- eBayes(fit)
fit$genes <- anno_df
gene_list <- topTable(fit,
                      coef          = 2,
                      adjust.method ="none",
                      number        = Inf,
                      sort.by       = "P",
                      p.value       = 0.01)
##'-----------------------------------------------------------------------------------------#
