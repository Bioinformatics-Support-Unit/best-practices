#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Interferon Study of a TET2 Mutation in pediatric Patients - Rare Samples   |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
#  Description : Illumina HT-12v4 Microarray Interferon Study - Utilitses a grouped model   |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Hambleton/Raw_Data/")

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
##'-----------------------------------------------------------------------------------------#



##'Read in Raw Data and attach pheno data
##'-----------------------------------------------------------------------------------------#
raw_data              <- lumiR("Hambleton_Sample_Probe_Profile.txt", verbose=T)
pheno_table           <- read.table("pheno.txt", header=T, sep="\t", stringsAsFactors=F)
pheno_table$batches   <- as.numeric(as.factor(pheno_table$Line_ID))
rownames(pheno_table) <- pheno_table$ArrayID
pData(raw_data)       <- pheno_table
##'-----------------------------------------------------------------------------------------#



##'Normalisation - VST and RSN
##'-----------------------------------------------------------------------------------------#
vst_data               <- lumiT(raw_data, method='vst')
rsn_data               <- lumiN(vst_data, method="rsn")
lumi.Q                 <- lumiQ(rsn_data)
##'-----------------------------------------------------------------------------------------#



##'Probe level filtering - Detection P Values
##'-----------------------------------------------------------------------------------------#
exprs_data             <- exprs(lumi.Q)
present_count          <- detectionCall(lumi.Q)
normalised_data        <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#



##'Probe level filtering - Detection P Values
##'-----------------------------------------------------------------------------------------#
probe_list             <- rownames(normalised_data)
nuIDs                  <- probeID2nuID(probe_list)[, "nuID"]
symbol                 <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name                   <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df                <- data.frame(ID=nuIDs, probe_list, symbol, name)
##'-----------------------------------------------------------------------------------------#



##'Limma Model Design - Multi-Level Analysis
##'-----------------------------------------------------------------------------------------#
foo              <- data.frame(ArrayName=colnames(normalised_data))
foo$CellLine     <- pData(lumi.Q)$Line_ID
foo$Type         <- "Stock"
foo$Type[grep("_1|_2|_3", pData(lumi.Q)$Sample_Name)] <- "Patient"
foo$Treatment    <- pData(lumi.Q)$Treatment

Treat            <- factor(paste(foo$Type,foo$Treatment, sep="."))
design           <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

corfit           <- duplicateCorrelation(normalised_data, design, block=foo$CellLine)
fit              <- lmFit(normalised_data,
                          design,
                          block=foo$CellLine,
                          correlation=corfit$consensus)

cm               <- makeContrasts(ControlStockPatient      = Patient.C-Stock.C,
                                  AlphaStockPatient        = Patient.A-Stock.A,
                                  BetaStockPatient         = Patient.B-Stock.B,
                                  GammaStockPatient        = Patient.G-Stock.G,

                                  ZeroAlphaStock           = Stock.A-Stock.C,
                                  ZeroGammaStock           = Stock.G-Stock.C,
                                  ZeroAlphaPatient         = Patient.A-Patient.C,
                                  ZeroBetaPatient          = Patient.B-Patient.C,
                                  ZeroGammaPatient         = Patient.G-Patient.C,
                                  ZeroBetaStock            = Stock.B-Stock.C,

                                  AlphaPatientStockControl = Patient.A-Stock.C,
                                  BetaPatientStockControl  = Patient.B-Stock.C,

                                  levels=design)
fit2            <- contrasts.fit(fit, cm)
fit2            <- eBayes(fit2)
fit2$genes      <- anno_df

# pData(raw_data_det)$ggCont <- Treat
# pData(lumi.Q)$ggCont       <- Treat
##'-----------------------------------------------------------------------------------------#



##'Model Interogation - Multi Level Analysis
##'-----------------------------------------------------------------------------------------#
comparisons <- c("ControlStockPatient",      "AlphaStockPatient",
                 "BetaStockPatient",         "GammaStockPatient",
                 "ZeroAlphaStock",           "ZeroGammaStock",
                 "ZeroAlphaPatient",         "ZeroBetaPatient",
                 "ZeroGammaPatient",         "ZeroBetaStock",
                 "AlphaPatientStockControl", "BetaPatientStockControl")

p_cut_off   <- 0.01
fold_change <- 2
mtc         <- 'BH'
i           <- 1

gene_list            <- topTable(fit2,
                                 coef=comparisons[i],
                                 p.value=p_cut_off,
                                 lfc=log2(fold_change),
                                 number=Inf,
                                 adjust.method="BH")

gene_list_unfiltered <- topTable(fit2,
                                 coef=comparisons[i],
                                 number=Inf,
                                 adjust.method="BH")
##'-----------------------------------------------------------------------------------------#
