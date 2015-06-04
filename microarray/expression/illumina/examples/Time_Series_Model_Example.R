#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : si Knockout and Interlukin Treatment Time Series                           |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Illumina HT-12v4 Microarray of si Knockout and Interlukin time series      |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("/Volumes/MRG_DATA/Bioinformatics/project_backups/Dimitra_Dec_2014/raw_data/")

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
raw.lumi               <- lumiR.idat(files=list.files("R/"),
                                     path=file.path(getwd(), "R/"),
                                     probeID="ProbeID",
                                     manifestfile="HumanHT-12_V4_0_R2_15002873_B.txt",
                                     controls=TRUE,
                                     detectionTh = 0.01,
                                     backgroundCorrect=T,
                                     collapseMode="none",
                                     QC=TRUE,
                                     memory="-Xmx4080m",
                                     verbose=T)
# raw.lumi.bgc           <- lumiB(raw.lumi, method="bgAdjust", verbose=T)
pheno_table            <- read.table("pheno.txt", header=T, row.names=1,
                                     sep="\t", stringsAsFactors=F)
raw.lumi.bgc           <- raw.lumi
pData(raw.lumi.bgc)    <- pheno_table
##'-----------------------------------------------------------------------------------------#



##'Normalisation - VST and RSN
##'-----------------------------------------------------------------------------------------#
vst_data               <- lumiT(raw.lumi.bgc, method='vst')
rsn_data               <- lumiN(vst_data, method = "rsn")
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
anno_df                <- data.frame(ID=nuIDs,
                                     probe_list,
                                     symbol,
                                     name)
##'-----------------------------------------------------------------------------------------#



##'Limma Model Design - Time Series
##'-----------------------------------------------------------------------------------------#
pData(lumi.Q)$TimeSeries      <- paste(pData(lumi.Q)$Sample_ID,
                                       pData(lumi.Q)$Treatment,
                                       sep=".")
lev                           <- c("Apicidin.no_IL1","Apicidin.IL1_1h","Apicidin.IL1_6h",
                                   "DMSO.no_IL1",    "DMSO.IL1_1h",    "DMSO.IL1_6h",
                                   "siHDAC3.no_IL1", "siHDAC3.IL1_1h", "siHDAC3.IL1_6h",
                                   "siCon.no_IL1",   "siCon.IL1_1h",   "siCon.IL1_6h")
f                             <- factor(pData(lumi.Q)$TimeSeries,
                                        levels=lev)
design                        <- model.matrix(~0+f)
colnames(design)              <- lev
fit                           <- lmFit(normalised_data, design)
cont.wt                       <- makeContrasts(Apicidin.IL1_1h-Apicidin.no_IL1,
                                               Apicidin.IL1_6h-Apicidin.IL1_1h,
                                               Apicidin.IL1_6h-Apicidin.no_IL1,

                                               DMSO.IL1_1h-DMSO.no_IL1,
                                               DMSO.IL1_6h-DMSO.IL1_1h,
                                               DMSO.IL1_6h-DMSO.no_IL1,

                                               siHDAC3.IL1_1h-siHDAC3.no_IL1,
                                               siHDAC3.IL1_6h-siHDAC3.IL1_1h,
                                               siHDAC3.IL1_6h-siHDAC3.no_IL1,

                                               siCon.IL1_1h-siCon.no_IL1,
                                               siCon.IL1_6h-siCon.IL1_1h,
                                               siCon.IL1_6h-siCon.no_IL1,
                                               levels=design)
fit2                          <- contrasts.fit(fit, cont.wt)
fit2                          <- eBayes(fit2)
fit2$genes                    <- anno_df
##'-----------------------------------------------------------------------------------------#



##'Model Interogation - Time Series
##'-----------------------------------------------------------------------------------------#
comparisons <- c("Apicidin.IL1_1h - Apicidin.no_IL1", "Apicidin.IL1_6h - Apicidin.IL1_1h",
                 "Apicidin.IL1_6h - Apicidin.no_IL1", "DMSO.IL1_1h - DMSO.no_IL1",
                 "DMSO.IL1_6h - DMSO.IL1_1h",         "DMSO.IL1_6h - DMSO.no_IL1",
                 "siHDAC3.IL1_1h - siHDAC3.no_IL1",   "siHDAC3.IL1_6h - siHDAC3.IL1_1h",
                 "siHDAC3.IL1_6h - siHDAC3.no_IL1",   "siCon.IL1_1h - siCon.no_IL1",
                 "siCon.IL1_6h - siCon.IL1_1h",       "siCon.IL1_6h - siCon.no_IL1")

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



##'Volcano Plots - Time Series
##'-----------------------------------------------------------------------------------------#
for(j in 1:length(comparisons)) {
  png(paste0("VPlot_", gsub(" ", "", comparisons[j], fixed=T), ".png"),
      width=4096,
      height=3096,
      units="px",
      res=300)
    gene_list_unfiltered <- topTable(fit2,
                                     coef=comparisons[j],
                                     number=Inf,
                                     adjust.method="BH")
    print(ggvolcano(gene_list_unfiltered,
                    contrast=comparisons[j]))
  dev.off()
}
##'-----------------------------------------------------------------------------------------#
