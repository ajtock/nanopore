#!/usr/bin/env Rscript

# Analysis:
# Order genes by decreasing Fleiss' kappa or stochasticity and make
# heatmaps of gene classifiers using this ordering

# Usage:
# conda activate R-4.0.3
# ./feature_among_read_variation_scoring_func_genes_heatmap.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions'
# conda deactivate

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "regions"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]

options(stringsAsFactors = F)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(circlize)
library(data.table)
library(stringr)
library(doParallel)
library(doFuture)
registerDoFuture()
plan(multicore)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

#library(parallel)
#library(dplyr)
#library(GenomicRanges)
##BiocManager::install("clusterProfiler")
##library(clusterProfiler)
##BiocManager::install("org.At.tair.db", character.only = TRUE)
##library("org.At.tair.db", character.only = TRUE)
##BiocManager::install("pathview") # installation of package ‘Rgraphviz’ had non-zero exit status; installation of package ‘pathview’ had non-zero exit status
##library(pathview)
##BiocManager::install("enrichplot")
##library(enrichplot)
##BiocManager::install("DOSE")
##library(DOSE)
#library(ggplot2)
#
#library(methods)
#library(plotrix)
#library(ggbeeswarm)
#library(ggthemes)
#library(grid)
#library(gridExtra)
#library(extrafont)
#
##library(parallel)
##library(GenomicRanges)
##library(irr)
##library(dplyr)
##library(tidyr)
##library(cluster)
##library(fpc)
###library(data.table)
###library(segmentSeq)
##library(ComplexHeatmap)
###library(RColorBrewer)
##library(scales)
###library(circlize)
## 
##library(ggplot2)
##library(cowplot)
###library(ggcorrplot)
##library(viridis)
##library(ggthemes)
##library(tidyquant)
###library(grid)

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
plotDir_kappa_mC <- paste0(outDir, "plots/heatmaps_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/heatmaps_", context, "_stocha_mC/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/heatmaps_", context, "_kappa_stocha/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_kappa_mC, " ] || mkdir -p ", plotDir_kappa_mC))
system(paste0("[ -d ", plotDir_stocha_mC, " ] || mkdir -p ", plotDir_stocha_mC))
#system(paste0("[ -d ", plotDir_kappa_stocha, " ] || mkdir -p ", plotDir_kappa_stocha))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Read in Lloyd et al. (2015) Plant Cell gene classifiers
ds1 <- read.csv(paste0("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/",
                       "TPC2015-00051-LSBR3_Supplemental_Data_set_1_Sheet1.csv"),
                header = T)
nrow(ds1)

lethal <- ds1[ds1[,2] == "Lethal" |
              ds1[,3] == "Yes",]
nonlethal <- ds1[ds1[,2] == "Non-Lethal" |
                 ds1[,3] == "No",]

ds3 <- read.csv(paste0("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/",
                       "TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv"),
                header = T, na.strings = c("NA", "?"))
nrow(ds3)

tandemdup <- ds3[which(ds3$Tandem.duplicate == 1),]
nrow(tandemdup)



# Load among-read and within-read mC data for featName featRegion
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)
featDF$parent <- sub(pattern = "\\.\\d+", replacement = "", x = featDF$name)
featDF$parent <- sub(pattern = "_\\d+", replacement = "", x = featDF$parent)

# Append intron retention ratio (calculated with IRFinder)
Col_Rep1_IRFinder <- fread(paste0("/home/ajt200/analysis/RNAseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_RNAseq_IRFinder_TAIR10_chr_all/REF/TAIR10_chr_all/",
                                  "Col_0_RNAseq_pooled_ERR96615/IRFinder-IR-dir.txt"),
                           sep = "\t", data.table = F)
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[grep("clean", Col_Rep1_IRFinder$Name),]
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[-which(Col_Rep1_IRFinder$Warnings %in% c("LowCover")),]
#nrow(Col_Rep1_IRFinder[which(Col_Rep1_IRFinder$Warnings == "-"),])
#[1] 45907
#[1] 22136
Col_Rep1_IRFinder$Name <- str_extract(Col_Rep1_IRFinder$Name, "AT\\wG\\d+")
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[-which(is.na(Col_Rep1_IRFinder$Name)),]

parentIDs <- unique(Col_Rep1_IRFinder$Name)

Col_Rep1_IRratio <- foreach(i = iter(parentIDs),
                            .combine = "rbind",
                            .multicombine = T,
                            .maxcombine = length(parentIDs)+1e1,
                            .inorder = F,
                            .errorhandling = "pass") %dopar% {
  tmpDF <- Col_Rep1_IRFinder[which(Col_Rep1_IRFinder$Name == i),]
  data.frame(chr = paste0("Chr", tmpDF$Chr[1]),
             start = min(tmpDF$Start, na.rm = T),
             end = max(tmpDF$End, na.rm = T),
             parent = i,
             strand = tmpDF$Strand[1],
             intronWidth_sum = sum(tmpDF$End - tmpDF$Start + 1, na.rm = T),
             excludedBases_sum = sum(tmpDF$ExcludedBases, na.rm = T),
             coverage_sum = sum(tmpDF$Coverage, na.rm = T),
             intronDepth_sum = sum(tmpDF$IntronDepth, na.rm = T),
             IRratio_mean = mean(tmpDF$IRratio, na.rm = T),
             IRratio_median = median(tmpDF$IRratio, na.rm = T),
             IRratio_sd = sd(tmpDF$IRratio, na.rm = T),
             IRratio_min = min(tmpDF$IRratio, na.rm = T),
             IRratio_max = max(tmpDF$IRratio, na.rm = T))
}

featDF_tab <- base::merge(x = featDF, y = Col_Rep1_IRratio,
                          by.x = "parent", by.y = "parent")


