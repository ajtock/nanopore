#!/usr/bin/env Rscript

# Analysis:
# Order TEs by decreasing Fleiss' kappa or stochasticity and make
# heatmaps of TE classifiers using this ordering

# Usage:
# conda activate R-4.0.3
# ./feature_among_read_variation_scoring_func_TEs_unfilt_heatmap.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies' 1
# conda deactivate

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CHG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "TE"
#featRegion <- "bodies"
#featMinLen <- 1

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]
featMinLen <- as.numeric(args[8])

options(stringsAsFactors = F)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(pals)
library(circlize)
library(data.table)
library(stringr)
library(GenomicRanges)
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
plotDir_kappa <- paste0(outDir, "plots/heatmaps_", context, "_kappa/")
plotDir_stocha <- paste0(outDir, "plots/heatmaps_", context, "_stocha/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/heatmaps_", context, "_kappa_stocha/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_kappa, " ] || mkdir -p ", plotDir_kappa))
system(paste0("[ -d ", plotDir_stocha, " ] || mkdir -p ", plotDir_stocha))
#system(paste0("[ -d ", plotDir_kappa_stocha, " ] || mkdir -p ", plotDir_kappa_stocha))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Load among-read and within-read mC data for featName featRegion
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)
featDF <- featDF[which(featDF$feature_width >= featMinLen),]
colnames(featDF)[which(colnames(featDF) == "fk_kappa_all")] <- "Kappa"
colnames(featDF)[which(colnames(featDF) == "mean_stocha_all")] <- "Stocha"
colnames(featDF)[which(colnames(featDF) == "mean_min_acf_all")] <- "Min_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mean_acf_all")] <- "Mean_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mC_all")] <- paste0("Mean_m", context)
featDF$Kappa_C_density <- featDF$fk_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$Stocha_C_density <- featDF$stocha_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$parent <- featDF$name

##
#gypsy <- featDF[which(featDF$score == "Gypsy_LTR"),]
#cor.test(gypsy$Kappa, gypsy$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
#cor.test(gypsy$Stocha, gypsy$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
#copia <- featDF[which(featDF$score == "Copia_LTR"),]
#cor.test(copia$Kappa, copia$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
#cor.test(copia$Stocha, copia$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
#unltr <- featDF[which(featDF$score == "Unclassified_LTR"),]
#cor.test(unltr$Kappa, unltr$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
#cor.test(unltr$Stocha, unltr$ltr_identity, use = "pairwise.complete.obs", method = "spearman")
##

featGR <- GRanges(seqnames = featDF$chr,
                  ranges = IRanges(start = featDF$start, end = featDF$end),
                  strand = "*",
                  featID = featDF$name)

superfamNames <- sort(unique(featDF$score))
superfamNames <- c(superfamNames[3], superfamNames[1], superfamNames[14],
                   superfamNames[10], superfamNames[7], superfamNames[13],
                   superfamNames[6], superfamNames[2], superfamNames[4],
                   superfamNames[5], superfamNames[8], superfamNames[9],
                   superfamNames[12], superfamNames[11])
superfamNamesPlot <- gsub("Pogo_Tc1_Mariner", "Pogo/Tc1/Mar", superfamNames)
superfamNamesPlot <- gsub("_", " ", superfamNamesPlot)
superfamNamesPlot <- gsub("classified", ".", superfamNamesPlot)

featDF_merge <- featDF
for(x in 1:length(superfamNames)) {
  featDF_merge[, superfamNamesPlot[x]] <- NA
  featDF_merge[which(featDF_merge$score == superfamNames[x]), superfamNamesPlot[x]] <- 1
  featDF_merge[which(featDF_merge$score != superfamNames[x]), superfamNamesPlot[x]] <- 0
}


# Load coordinates for mitochondrial insertion on Chr2, in BED format,
# for removal of overlapping epimutation hotspots and coldspots
mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
                       header = F)
colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
mito_ins_chr <- unique(mito_ins$chr)
mito_ins_start <- min(mito_ins$start)+1
mito_ins_end <- max(mito_ins$end)
#mito_ins_GR <- GRanges(seqnames = "Chr2",
#                       ranges = IRanges(start = min(mito_ins$start)+1,
#                                        end = max(mito_ins$end)),
#                       strand = "*")

# Read in methylation divergence (mD) files
mDdir <- "/home/ajt200/analysis/BSseq_leaf_Hazarika_Johannes_2022_NatPlants/snakemake_BSseq_t2t-col.20210610/MA1_2/coverage/report/alphabeta/"
MA1_2_mD_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0(mDdir,
                    "mD_at_dt62_genomeBinSize10kb_genomeStepSize1kb",
                    "_MA1_2_MappedOn_", refbase, "_", chrName[x], "_CpG.tsv"),
             header = T)
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  MA1_2_mD <- dplyr::bind_rows(MA1_2_mD_list)
} else {
  MA1_2_mD <- MA1_2_mD_list[[1]]
}
rm(MA1_2_mD_list); gc()

tab_mD <- MA1_2_mD
rm(MA1_2_mD); gc()

tab_mD <- tab_mD[
  with(tab_mD, order(chr, start, end)),
]

# Remove genomic bins in tab_mD overlapping the mitochondrial insertion on Chr2,
# as we cannot be sure that these BS-seq reads come from the nuclear genome
# Note that long reads wholly contained within the boundaries of the Chr2 mito insertion
# have already been removed for the among-read and within-read mC variation analysis,
# by among_read_variation_scoring.R
tab_mD_mito <- tab_mD %>%
  dplyr::filter(chr == mito_ins_chr) %>%
  dplyr::filter( (start >= mito_ins_start &
                  start <= mito_ins_end) |
                 (end >= mito_ins_start &
                  end <= mito_ins_end) )

stopifnot(unique(tab_mD_mito$chr) == mito_ins_chr)
stopifnot(tab_mD_mito$end[1] >= mito_ins_start)
stopifnot(tab_mD_mito$start[nrow(tab_mD_mito)] <= mito_ins_end)
print(head(tab_mD_mito[,1:6]))
print(tail(tab_mD_mito[,1:6]))
print(paste0("Mitochondrial insetion on ", mito_ins_chr, ":", mito_ins_start, "-", mito_ins_end))

tab_mD <- tab_mD %>%
  dplyr::anti_join(tab_mD_mito)

tab_mD <- tab_mD[!is.na(tab_mD$MA1_2_mean.D),]

tab_mD$rank <- rank(tab_mD$MA1_2_mean.D)/nrow(tab_mD)

# Define epimutation hotspots (mD_hs) and coldspots (mD_cs)
mD_hs <- tab_mD[tab_mD$rank >= 0.9,]
mD_cs <- tab_mD[tab_mD$rank <= 0.1,]

mD_hs_GR <- GRanges(seqnames = mD_hs$chr,
                    ranges = IRanges(start = mD_hs$start, end = mD_hs$end),
                    strand = "*",
                    mD = mD_hs$MA1_2_mean.D)
mD_cs_GR <- GRanges(seqnames = mD_cs$chr,
                    ranges = IRanges(start = mD_cs$start, end = mD_cs$end),
                    strand = "*",
                    mD = mD_cs$MA1_2_mean.D)

fOverlaps_feat_mD_hs <- findOverlaps(query = featGR,
                                     subject = mD_hs_GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = T)
featGR_mD_hs <- featGR[unique(queryHits(fOverlaps_feat_mD_hs))]
featID_mD_hs <- unique(featGR_mD_hs$featID)

fOverlaps_feat_mD_cs <- findOverlaps(query = featGR,
                                     subject = mD_cs_GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = T)
featGR_mD_cs <- featGR[unique(queryHits(fOverlaps_feat_mD_cs))]
featID_mD_cs <- unique(featGR_mD_cs$featID)

featDF_merge$mD_hotspot <- NA
featDF_merge[which(featDF_merge$name %in% featID_mD_hs),]$mD_hotspot <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_mD_hs)),]$mD_hotspot <- 0
featDF_merge$mD_coldspot <- NA
featDF_merge[which(featDF_merge$name %in% featID_mD_cs),]$mD_coldspot <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_mD_cs)),]$mD_coldspot <- 0

# DMRs in different DNA methylation mutants
# kss_hypoCHG
kss_hypoCHG <- read.table(paste0("/home/ajt200/analysis/BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/",
                                 "snakemake_BSseq_", refbase, "/coverage/report/DMRs/hypoDMRs/",
                                 paste0(chrName, collapse = "_"),
                                 "/features_6quantiles_by_change_in_kss_BSseq_Rep1_hypoCHG_DMRs_vs3reps_",
                                 "mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200_in_",
                                 refbase, "_", paste0(chrName, collapse = "_"), "_genomewide.bed"),
                          header = F)
kss_hypoCHG_GR <- GRanges(seqnames = kss_hypoCHG[,1],
                          ranges = IRanges(start = kss_hypoCHG[,2]+1,
                                           end = kss_hypoCHG[,3]),
                          strand = "*",
                          l2fc = kss_hypoCHG[,5])
fOverlaps_feat_kss_hypoCHG <- findOverlaps(query = featGR,
                                           subject = kss_hypoCHG_GR,
                                           type = "any",
                                           select = "all",
                                           ignore.strand = T)
featGR_kss_hypoCHG <- featGR[unique(queryHits(fOverlaps_feat_kss_hypoCHG))]
featID_kss_hypoCHG <- unique(featGR_kss_hypoCHG$featID)

featDF_merge$kss_hypoCHG <- NA
featDF_merge[which(featDF_merge$name %in% featID_kss_hypoCHG),]$kss_hypoCHG <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_kss_hypoCHG)),]$kss_hypoCHG <- 0

# cmt3_hypoCHG
cmt3_hypoCHG <- read.table(paste0("/home/ajt200/analysis/BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/",
                                  "snakemake_BSseq_", refbase, "/coverage/report/DMRs/hypoDMRs/",
                                  paste0(chrName, collapse = "_"),
                                  "/features_6quantiles_by_change_in_cmt3_BSseq_Rep1_hypoCHG_DMRs_vs3reps_",
                                  "mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200_in_",
                                  refbase, "_", paste0(chrName, collapse = "_"), "_genomewide.bed"),
                           header = F)
cmt3_hypoCHG_GR <- GRanges(seqnames = cmt3_hypoCHG[,1],
                           ranges = IRanges(start = cmt3_hypoCHG[,2]+1,
                                            end = cmt3_hypoCHG[,3]),
                           strand = "*",
                           l2fc = cmt3_hypoCHG[,5])
fOverlaps_feat_cmt3_hypoCHG <- findOverlaps(query = featGR,
                                            subject = cmt3_hypoCHG_GR,
                                            type = "any",
                                            select = "all",
                                            ignore.strand = T)
featGR_cmt3_hypoCHG <- featGR[unique(queryHits(fOverlaps_feat_cmt3_hypoCHG))]
featID_cmt3_hypoCHG <- unique(featGR_cmt3_hypoCHG$featID)

featDF_merge$cmt3_hypoCHG <- NA
featDF_merge[which(featDF_merge$name %in% featID_cmt3_hypoCHG),]$cmt3_hypoCHG <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_cmt3_hypoCHG)),]$cmt3_hypoCHG <- 0

# kss_hypoCHH
kss_hypoCHH <- read.table(paste0("/home/ajt200/analysis/BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/",
                                 "snakemake_BSseq_", refbase, "/coverage/report/DMRs/hypoDMRs/",
                                 paste0(chrName, collapse = "_"),
                                 "/features_6quantiles_by_change_in_kss_BSseq_Rep1_hypoCHH_DMRs_vs3reps_",
                                 "mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200_in_",
                                 refbase, "_", paste0(chrName, collapse = "_"), "_genomewide.bed"),
                          header = F)
kss_hypoCHH_GR <- GRanges(seqnames = kss_hypoCHH[,1],
                          ranges = IRanges(start = kss_hypoCHH[,2]+1,
                                           end = kss_hypoCHH[,3]),
                          strand = "*",
                          l2fc = kss_hypoCHH[,5])
fOverlaps_feat_kss_hypoCHH <- findOverlaps(query = featGR,
                                           subject = kss_hypoCHH_GR,
                                           type = "any",
                                           select = "all",
                                           ignore.strand = T)
featGR_kss_hypoCHH <- featGR[unique(queryHits(fOverlaps_feat_kss_hypoCHH))]
featID_kss_hypoCHH <- unique(featGR_kss_hypoCHH$featID)

featDF_merge$kss_hypoCHH <- NA
featDF_merge[which(featDF_merge$name %in% featID_kss_hypoCHH),]$kss_hypoCHH <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_kss_hypoCHH)),]$kss_hypoCHH <- 0

# cmt3_hypoCHH
cmt3_hypoCHH <- read.table(paste0("/home/ajt200/analysis/BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/",
                                  "snakemake_BSseq_", refbase, "/coverage/report/DMRs/hypoDMRs/",
                                  paste0(chrName, collapse = "_"),
                                  "/features_6quantiles_by_change_in_cmt3_BSseq_Rep1_hypoCHH_DMRs_vs3reps_",
                                  "mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200_in_",
                                  refbase, "_", paste0(chrName, collapse = "_"), "_genomewide.bed"),
                           header = F)
cmt3_hypoCHH_GR <- GRanges(seqnames = cmt3_hypoCHH[,1],
                           ranges = IRanges(start = cmt3_hypoCHH[,2]+1,
                                            end = cmt3_hypoCHH[,3]),
                           strand = "*",
                           l2fc = cmt3_hypoCHH[,5])
fOverlaps_feat_cmt3_hypoCHH <- findOverlaps(query = featGR,
                                            subject = cmt3_hypoCHH_GR,
                                            type = "any",
                                            select = "all",
                                            ignore.strand = T)
featGR_cmt3_hypoCHH <- featGR[unique(queryHits(fOverlaps_feat_cmt3_hypoCHH))]
featID_cmt3_hypoCHH <- unique(featGR_cmt3_hypoCHH$featID)

featDF_merge$cmt3_hypoCHH <- NA
featDF_merge[which(featDF_merge$name %in% featID_cmt3_hypoCHH),]$cmt3_hypoCHH <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_cmt3_hypoCHH)),]$cmt3_hypoCHH <- 0

DMRnames <- colnames(featDF_merge)[(ncol(featDF_merge)-3):ncol(featDF_merge)]
DMRnamesPlot <- gsub("_", " ", DMRnames)

featDF_kappa <- featDF_merge[
  with(featDF_merge, 
       order(Kappa, decreasing = T)),
]

featDF_stocha <- featDF_merge[
  with(featDF_merge, 
       order(Stocha, decreasing = T)),
]

# Define heatmap colours
rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
rich8 <- function() {manual_pal(values = c("#000041","#0000CB","#0081FF","#02DA81","#80FE1A","#FDEE02","#FFAB00","#FF3300"))}
rich6 <- function() {manual_pal(values = c("#000043","#0033FF","#01CCA4","#BAFF12","#FFCC00","#FF3300"))}
rich8to6equal <- c("#0000CB", "#0081FF", "#87CEFA", "#FDEE02", "#FFAB00", "#FF3300")
revSpectralScale11 <- rev(brewer.pal(11, "Spectral"))
viridisScale6 <- viridis_pal()(6)

# Heatmap plotting function
featureHeatmap <- function(mat,
                           colFun,
                           datName,
                           rowOrder) {
  Heatmap(matrix = mat,
          col = colFun,
          row_order = rowOrder,
#          column_title = datName,
#          column_title_rot = 45,
#          column_title_gp = gpar(fontsize = 13),
          column_labels = datName,
          column_names_gp = gpar(fontsize = 13),
          column_names_side = "top",
          column_names_rot = 45,
          column_names_centered = TRUE,
          cluster_columns = FALSE,
          cluster_column_slices = FALSE,
          cluster_rows = FALSE,
          cluster_row_slices = FALSE,
          heatmap_legend_param = list(title = datName,
                                      title_position = "topcenter",
                                      title_gp = gpar(font = 1, fontsize = 11),
                                      legend_direction = "horizontal",
                                      labels_gp = gpar(fontsize = 10)),
          heatmap_width = unit(2, "npc"),
          heatmap_height = unit(4, "npc"),
          column_gap = unit(0, "mm"),
          row_gap = unit(1.0, "mm"),
#          row_split = factor(tab_extend$phylo,
#                             levels = sort(unique(as.character(tab_extend$phylo)))),
          row_title = NULL,
          show_row_names = FALSE,
          border = FALSE,
          # If converting into png with pdfTotiffTopng.sh,
          # set use_raster to FALSE
          use_raster = FALSE)
          #use_raster = TRUE, raster_device = "png", raster_quality = 4)
}

Kappa_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Kappa")])
Stocha_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Stocha")])
Kappa_C_density_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Kappa_C_density")])
Stocha_C_density_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Stocha_C_density")])
Mean_mC_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == paste0("Mean_m", context))])
feature_width_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "feature_width")])
ltr_identity_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "ltr_identity")])
score_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "score")])
superfam_mat_list <- lapply(superfamNamesPlot, function(x) {
  as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == x)])
})
DNA_RNA_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "DNA_RNA")])
DMR_mat_list <- lapply(DMRnames, function(x) {
  as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == x)])
})
mD_hotspot_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "mD_hotspot")])
mD_coldspot_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "mD_coldspot")])

Kappa_colFun <- colorRamp2(quantile(
    Kappa_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
Stocha_colFun <- colorRamp2(quantile(
    Stocha_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  plasma(6))
Kappa_C_density_colFun <- colorRamp2(quantile(
    Kappa_C_density_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
Stocha_C_density_colFun <- colorRamp2(quantile(
    Stocha_C_density_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  plasma(6))
Mean_mC_colFun <- colorRamp2(quantile(
    Mean_mC_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
feature_width_colFun <- colorRamp2(quantile(
    feature_width_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  plasma(6))
ltr_identity_colFun <- colorRamp2(quantile(
    ltr_identity_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  rev(viridis(6)))
score_colFun <- cols25(n = 25)[-c(7:16, 25)]
names(score_colFun) <- superfamNames
superfam_colFun_list <- lapply(1:length(superfam_mat_list), function(x) {
  c("0" = "white", "1" = cols25(n = 25)[-c(7:16, 25)][x])
})
DNA_RNA_colFun <- c("RNA" = "darkorange1", "DNA" = "dodgerblue4", "Unclassified" = "grey80")
DMR_colFun_list <- lapply(1:length(DMR_mat_list), function(x) {
 c("0" = "white", "1" = cols25(n = 4)[x])
})
mD_hotspot_colFun <- c("0" = "white", "1" = "firebrick3") 
mD_coldspot_colFun <- c("0" = "white", "1" = "dodgerblue4") 


Kappa_htmp <- featureHeatmap(mat = Kappa_mat,
  colFun = Kappa_colFun,
  datName = "Agreement",
  rowOrder = c(1:nrow(Kappa_mat)))
Stocha_htmp <- featureHeatmap(mat = Stocha_mat,
  colFun = Stocha_colFun,
  datName = "Stochasticity",
  rowOrder = c(1:nrow(Stocha_mat)))
Kappa_C_density_htmp <- featureHeatmap(mat = Kappa_C_density_mat,
  colFun = Kappa_C_density_colFun,
  datName = paste0(context, " density kappa"),
  rowOrder = c(1:nrow(Kappa_C_density_mat)))
Stocha_C_density_htmp <- featureHeatmap(mat = Stocha_C_density_mat,
  colFun = Stocha_C_density_colFun,
  datName = paste0(context, " density"),
  rowOrder = c(1:nrow(Stocha_C_density_mat)))
Mean_mC_htmp <- featureHeatmap(mat = Mean_mC_mat,
  colFun = Mean_mC_colFun,
  datName = paste0("m", context, " avg."),
  rowOrder = c(1:nrow(Mean_mC_mat)))
feature_width_htmp <- featureHeatmap(mat = feature_width_mat,
  colFun = feature_width_colFun,
  datName = "TE length",
  rowOrder = c(1:nrow(feature_width_mat)))
ltr_identity_htmp <- featureHeatmap(mat = ltr_identity_mat,
  colFun = ltr_identity_colFun,
  datName = "LTR identity",
  rowOrder = c(1:nrow(ltr_identity_mat)))
score_htmp <- featureHeatmap(mat = score_mat,
  colFun = score_colFun,
  datName = "Superfamily",
  rowOrder = c(1:nrow(score_mat)))
superfam_htmp_list <- lapply(1:length(superfam_mat_list), function(x) {
  featureHeatmap(mat = superfam_mat_list[[x]],
                 colFun = superfam_colFun_list[[x]],
                 datName = superfamNamesPlot[x],
                 rowOrder = c(1:nrow(superfam_mat_list[[x]])))
})
DNA_RNA_htmp <- featureHeatmap(mat = DNA_RNA_mat,
  colFun = DNA_RNA_colFun,
  datName = "DNA/RNA",
  rowOrder = c(1:nrow(DNA_RNA_mat)))
DMR_htmp_list <- lapply(1:length(DMR_mat_list), function(x) {
  featureHeatmap(mat = DMR_mat_list[[x]],
                 colFun = DMR_colFun_list[[x]],
                 datName = DMRnamesPlot[x],
                 rowOrder = c(1:nrow(DMR_mat_list[[x]])))
})
mD_hotspot_htmp <- featureHeatmap(mat = mD_hotspot_mat,
  colFun = mD_hotspot_colFun,
  datName = expression(italic(D) ~ "hotspot"),
  rowOrder = c(1:nrow(mD_hotspot_mat)))
mD_coldspot_htmp <- featureHeatmap(mat = mD_coldspot_mat,
  colFun = mD_coldspot_colFun,
  datName = expression(italic(D) ~ "coldspot"),
  rowOrder = c(1:nrow(mD_coldspot_mat)))

htmps <- Kappa_htmp + Stocha_htmp +
         Kappa_C_density_htmp + Stocha_C_density_htmp +
         Mean_mC_htmp + feature_width_htmp +
#         ltr_identity_htmp + score_htmp +
         score_htmp +
         superfam_htmp_list[[1]] + superfam_htmp_list[[2]] + superfam_htmp_list[[3]] + superfam_htmp_list[[4]] + superfam_htmp_list[[5]] +
         superfam_htmp_list[[6]] + superfam_htmp_list[[7]] + superfam_htmp_list[[8]] + superfam_htmp_list[[9]] + superfam_htmp_list[[10]] +
         superfam_htmp_list[[11]] + superfam_htmp_list[[12]] + superfam_htmp_list[[13]] + superfam_htmp_list[[14]] +
         DNA_RNA_htmp +
         DMR_htmp_list[[1]] + DMR_htmp_list[[2]] + DMR_htmp_list[[3]] + DMR_htmp_list[[4]]
#         mD_hotspot_htmp + mD_coldspot_htmp

legendGap <- unit(20, "mm")

pdf(paste0(plotDir_kappa,
           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase, "_", context,
           "_NAmax", NAmax, "_all_unfilt_heatmap_kappa_", paste0(chrName, collapse = "_"), ".pdf"),
    width = 1.5*length(htmps), height = 10)
draw(htmps,
     gap = unit(1, "mm"),
     column_title = paste0("TEs ordered by decreasing among-read agreement (m", context, ") in TE ", featRegion), 
     column_title_gp = gpar(font = 2, fontsize = 16),
     heatmap_legend_side = "bottom",
     legend_gap = legendGap)
dev.off()

