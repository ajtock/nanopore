#!/usr/bin/env Rscript

# Analysis:
# Order genes by decreasing Fleiss' kappa or stochasticity and make
# heatmaps of gene classifiers using this ordering

# Usage:
# conda activate R-4.0.3
# ./feature_among_read_variation_scoring_func_genes_filt_heatmap.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions'
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

# Load coordinates for mitochondrial insertion on Chr2, in BED format
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

# Read in Lloyd et al. (2015) Plant Cell gene classifiers
ds1 <- read.csv(paste0("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/",
                       "TPC2015-00051-LSBR3_Supplemental_Data_set_1_Sheet1.csv"),
                header = T)
nrow(ds1)
colnames(ds1)[1] <- "parent"

lethal <- ds1[ds1[,2] == "Lethal" |
              ds1[,3] == "Yes",]
nonlethal <- ds1[ds1[,2] == "Non-Lethal" |
                 ds1[,3] == "No",]

ds3 <- read.csv(paste0("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/",
                       "TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv"),
                header = T, na.strings = c("NA", "?"))
nrow(ds3)
ds3 <- ds3[,c(1, 2, 4, 5, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 31, 34, 35, 41, 49, 50)]
ncol(ds3)
ds3 <- ds3[,c(1, 21, 2, 16, 12, 3, 7, 13, 5, 6, 4, 8, 14, 9, 15, 10, 11, 17, 18, 20, 22, 23, 19, 24)] 
ncol(ds3)
colnames(ds3) <- c("parent",
                   "gbM",
                   "Median_expression", "Expression_breadth", "Expression_variation", "Expression_correlation", "Expression_correlation_Ks_lt2",
                   "Coexpression_module_size", "AraNet_edges", "PPIs",
                   "Protein_domains", "Amino_acids",
                   "PI_with_paralog", "Ks_with_paralog", "KaKs_with_paralog",
                   "Nucleotide_diversity", "Sequence_conservation", "OrthoMCL_paralog_cluster_size", "Core_eukaryotic_gene",
                   "beta_gamma_WGD_paralog_retained", "alpha_WGD_paralog_retained",
                   "Tandem_duplicate", "Pseudogene_homolog_present", "No_homolog_in_rice") 

ds3$Lethal <- NA
ds3[which(ds3$parent %in% lethal$parent),]$Lethal <- 1
ds3[which(ds3$parent %in% nonlethal$parent),]$Lethal <- 0


# Read in methylation divergence (mD) files
mDdir <- "/home/ajt200/analysis/BSseq_leaf_Hazarika_Johannes_2022_NatPlants/snakemake_BSseq_t2t-col.20210610/MA1_2/coverage/report/alphabeta/"
MA1_2_mD_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0(mDdir,
                    "mD_at_dt62_genomeBinSize10kb_genomeStepSize1kb",
                    "_MA1_2_MappedOn_", refbase, "_", chrName[x], "_", context, ".tsv"),
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

# Load among-read and within-read mC data for featName featRegion
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)
colnames(featDF)[which(colnames(featDF) == "fk_kappa_all")] <- "Kappa"
colnames(featDF)[which(colnames(featDF) == "mean_stocha_all")] <- "Stocha"
colnames(featDF)[which(colnames(featDF) == "mean_min_acf_all")] <- "Min_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mean_acf_all")] <- "Mean_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mC_all")] <- paste0("Mean_m", context)
featDF$Kappa_C_density <- featDF$fk_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$Stocha_C_density <- featDF$stocha_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$parent <- sub(pattern = "\\.\\d+", replacement = "", x = featDF$name)
featDF$parent <- sub(pattern = "_\\d+", replacement = "", x = featDF$parent)

featGR <- GRanges(seqnames = featDF$chr,
                  ranges = IRanges(start = featDF$start, end = featDF$end),
                  strand = "*",
                  featID = featDF$name)

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

featDF_merge <- base::merge(x = featDF, y = Col_Rep1_IRratio,
                            by.x = "parent", by.y = "parent",
                            all.x = T)
featDF_merge <- base::merge(x = featDF_merge, y = ds3,
                            by.x = "parent", by.y = "parent",
                            all.x = T)

featDF_merge$mD_hotspot <- NA
featDF_merge[which(featDF_merge$name %in% featID_mD_hs),]$mD_hotspot <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_mD_hs)),]$mD_hotspot <- 0
featDF_merge$mD_coldspot <- NA
featDF_merge[which(featDF_merge$name %in% featID_mD_cs),]$mD_coldspot <- 1
featDF_merge[which(!(featDF_merge$name %in% featID_mD_cs)),]$mD_coldspot <- 0

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
gbM_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "gbM")])
Expression_breadth_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Expression_breadth")])
Expression_variation_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Expression_variation")])
Median_expression_mat <- log2(as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Median_expression")])+1)
Coexpression_module_size_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Coexpression_module_size")])
feature_width_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "feature_width")])
Amino_acids_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Amino_acids")])
Protein_domains_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Protein_domains")])
exons_count_per_kb_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "exons_count_per_kb")])
exons_width_prop_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "exons_width_prop")])
introns_count_per_kb_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "introns_count_per_kb")])
introns_width_prop_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "introns_width_prop")])
IRratio_median_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "IRratio_median")])
Lethal_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Lethal")])
Core_eukaryotic_gene_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Core_eukaryotic_gene")])
Tandem_duplicate_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Tandem_duplicate")])
OrthoMCL_paralog_cluster_size_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "OrthoMCL_paralog_cluster_size")])
Sequence_conservation_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Sequence_conservation")])
Nucleotide_diversity_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Nucleotide_diversity")])
Ks_with_paralog_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Ks_with_paralog")])
KaKs_with_paralog_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "KaKs_with_paralog")])
mD_hotspot_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "mD_hotspot")])
mD_coldspot_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "mD_coldspot")])

#Min_ACF_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Min_ACF")])
#AraNet_edges_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "AraNet_edges")])
#PPIs_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "PPIs")])
#Expression_correlation_Ks_lt2_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Expression_correlation_Ks_lt2")])
#Expression_correlation_mat <- as.matrix(featDF_kappa[,which(colnames(featDF_kappa) == "Expression_correlation")])

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
gbM_colFun <- c("0" = "dodgerblue4", "1" = "darkorange1") 
Expression_breadth_colFun <- colorRamp2(quantile(
    Expression_breadth_mat,
    c(0.05, 0.95),
    na.rm = T),
  c("blue", "red"))
Expression_variation_colFun <- colorRamp2(quantile(
    Expression_variation_mat,
    c(0.05, 0.95),
    na.rm = T),
  c("blue", "red"))
Median_expression_colFun <- colorRamp2(quantile(
    Median_expression_mat,
    c(0.05, 0.95),
    na.rm = T),
  c("blue", "red"))
Coexpression_module_size_colFun <- colorRamp2(quantile(
    Coexpression_module_size_mat,
    c(0.05, 0.95),
    na.rm = T),
  c("blue", "red"))
feature_width_colFun <- colorRamp2(quantile(
    feature_width_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
Amino_acids_colFun <- colorRamp2(quantile(
    Amino_acids_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
Protein_domains_colFun <- colorRamp2(quantile(
    Protein_domains_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
exons_count_per_kb_colFun <- colorRamp2(quantile(
    exons_count_per_kb_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
exons_width_prop_colFun <- colorRamp2(quantile(
    exons_width_prop_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
introns_count_per_kb_colFun <- colorRamp2(quantile(
    introns_count_per_kb_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
introns_width_prop_colFun <- colorRamp2(quantile(
    introns_width_prop_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
IRratio_median_colFun <- colorRamp2(quantile(
    IRratio_median_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  viridis(6))
Lethal_colFun <- c("0" = "white", "1" = "black") 
Core_eukaryotic_gene_colFun <- c("0" = "white", "1" = "firebrick3") 
Tandem_duplicate_colFun <- c("0" = "dodgerblue", "1" = "darkgreen") 
OrthoMCL_paralog_cluster_size_colFun <- colorRamp2(quantile(
    OrthoMCL_paralog_cluster_size_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
Sequence_conservation_colFun <- colorRamp2(quantile(
    Sequence_conservation_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
Nucleotide_diversity_colFun <- colorRamp2(quantile(
    Nucleotide_diversity_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
Ks_with_paralog_colFun <- colorRamp2(quantile(
    Ks_with_paralog_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
KaKs_with_paralog_colFun <- colorRamp2(quantile(
    KaKs_with_paralog_mat,
    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
    na.rm = T),
  heat.colors(6))
mD_hotspot_colFun <- c("0" = "white", "1" = "firebrick3") 
mD_coldspot_colFun <- c("0" = "white", "1" = "dodgerblue4") 


#Min_ACF_colFun <- colorRamp2(quantile(
#    Min_ACF_mat,
#    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
#    na.rm = T),
#  plasma(6))

#AraNet_edges_colFun <- colorRamp2(quantile(
#    AraNet_edges_mat,
#    c(0.05, 0.95),
#    na.rm = T),
#  c("blue", "red"))
#PPIs_colFun <- colorRamp2(quantile(
#    PPIs_mat,
#    c(0.05, 0.95),
#    na.rm = T),
#  c("blue", "red"))

#Expression_correlation_Ks_lt2_colFun <- colorRamp2(quantile(
#    Expression_correlation_Ks_lt2_mat,
#    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
#    na.rm = T),
#  magma(6))
#Expression_correlation_colFun <- colorRamp2(quantile(
#    Expression_correlation_mat,
#    c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
#    na.rm = T),
#  magma(6))


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
gbM_htmp <- featureHeatmap(mat = gbM_mat,
  colFun = gbM_colFun,
  datName = "gbM",
  rowOrder = c(1:nrow(gbM_mat)))
Expression_breadth_htmp <- featureHeatmap(mat = Expression_breadth_mat,
  colFun = Expression_breadth_colFun,
  datName = "Exp. breadth",
  rowOrder = c(1:nrow(Expression_breadth_mat)))
Expression_variation_htmp <- featureHeatmap(mat = Expression_variation_mat,
  colFun = Expression_variation_colFun,
  datName = "Exp. variation",
  rowOrder = c(1:nrow(Expression_variation_mat)))
Median_expression_htmp <- featureHeatmap(mat = Median_expression_mat,
  colFun = Median_expression_colFun,
  datName = "Exp. avg.",
  rowOrder = c(1:nrow(Median_expression_mat)))
Coexpression_module_size_htmp <- featureHeatmap(mat = Coexpression_module_size_mat,
  colFun = Coexpression_module_size_colFun,
  datName = "Exp. module size",
  rowOrder = c(1:nrow(Coexpression_module_size_mat)))
feature_width_htmp <- featureHeatmap(mat = feature_width_mat,
  colFun = feature_width_colFun,
  datName = "Gene length",
  rowOrder = c(1:nrow(feature_width_mat)))
Amino_acids_htmp <- featureHeatmap(mat = Amino_acids_mat,
  colFun = Amino_acids_colFun,
  datName = "Protein length",
  rowOrder = c(1:nrow(Amino_acids_mat)))
Protein_domains_htmp <- featureHeatmap(mat = Protein_domains_mat,
  colFun = Protein_domains_colFun,
  datName = "Protein domains",
  rowOrder = c(1:nrow(Protein_domains_mat)))
exons_count_per_kb_htmp <- featureHeatmap(mat = exons_count_per_kb_mat,
  colFun = exons_count_per_kb_colFun,
  datName = "Exons per kb",
  rowOrder = c(1:nrow(exons_count_per_kb_mat)))
exons_width_prop_htmp <- featureHeatmap(mat = exons_width_prop_mat,
  colFun = exons_width_prop_colFun,
  datName = "Proportion exonic",
  rowOrder = c(1:nrow(exons_width_prop_mat)))
introns_count_per_kb_htmp <- featureHeatmap(mat = introns_count_per_kb_mat,
  colFun = introns_count_per_kb_colFun,
  datName = "Introns per kb",
  rowOrder = c(1:nrow(introns_count_per_kb_mat)))
introns_width_prop_htmp <- featureHeatmap(mat = introns_width_prop_mat,
  colFun = introns_width_prop_colFun,
  datName = "Proportion intronic",
  rowOrder = c(1:nrow(introns_width_prop_mat)))
IRratio_median_htmp <- featureHeatmap(mat = IRratio_median_mat,
  colFun = IRratio_median_colFun,
  datName = "Intron retention",
  rowOrder = c(1:nrow(IRratio_median_mat)))
Lethal_htmp <- featureHeatmap(mat = Lethal_mat,
  colFun = Lethal_colFun,
  datName = "Lethal",
  rowOrder = c(1:nrow(Lethal_mat)))
Core_eukaryotic_gene_htmp <- featureHeatmap(mat = Core_eukaryotic_gene_mat,
  colFun = Core_eukaryotic_gene_colFun,
  datName = "Core eukaryotic gene",
  rowOrder = c(1:nrow(Core_eukaryotic_gene_mat)))
Tandem_duplicate_htmp <- featureHeatmap(mat = Tandem_duplicate_mat,
  colFun = Tandem_duplicate_colFun,
  datName = "Tandem duplicate",
  rowOrder = c(1:nrow(Tandem_duplicate_mat)))
OrthoMCL_paralog_cluster_size_htmp <- featureHeatmap(mat = OrthoMCL_paralog_cluster_size_mat,
  colFun = OrthoMCL_paralog_cluster_size_colFun,
  datName = "Ortho. cluster size",
  rowOrder = c(1:nrow(OrthoMCL_paralog_cluster_size_mat)))
Sequence_conservation_htmp <- featureHeatmap(mat = Sequence_conservation_mat,
  colFun = Sequence_conservation_colFun,
  datName = "Sequence conservation",
  rowOrder = c(1:nrow(Sequence_conservation_mat)))
Nucleotide_diversity_htmp <- featureHeatmap(mat = Nucleotide_diversity_mat,
  colFun = Nucleotide_diversity_colFun,
  datName = "Nucleotide diversity",
  rowOrder = c(1:nrow(Nucleotide_diversity_mat)))
Ks_with_paralog_htmp <- featureHeatmap(mat = Ks_with_paralog_mat,
  colFun = Ks_with_paralog_colFun,
  datName = "Ks with paralog",
  rowOrder = c(1:nrow(Ks_with_paralog_mat)))
KaKs_with_paralog_htmp <- featureHeatmap(mat = KaKs_with_paralog_mat,
  colFun = KaKs_with_paralog_colFun,
  datName = "Ka/Ks with paralog",
  rowOrder = c(1:nrow(KaKs_with_paralog_mat)))
mD_hotspot_htmp <- featureHeatmap(mat = mD_hotspot_mat,
  colFun = mD_hotspot_colFun,
  datName = expression(italic(D) ~ "hotspot"),
  rowOrder = c(1:nrow(mD_hotspot_mat)))
mD_coldspot_htmp <- featureHeatmap(mat = mD_coldspot_mat,
  colFun = mD_coldspot_colFun,
  datName = expression(italic(D) ~ "coldspot"),
  rowOrder = c(1:nrow(mD_coldspot_mat)))

#Min_ACF_htmp <- featureHeatmap(mat = Min_ACF_mat,
#  colFun = Min_ACF_colFun,
#  datName = "ACF min.",
#  rowOrder = c(1:nrow(Min_ACF_mat)))

#AraNet_edges_htmp <- featureHeatmap(mat = AraNet_edges_mat,
#  colFun = AraNet_edges_colFun,
#  datName = "AraNet edges",
#  rowOrder = c(1:nrow(AraNet_edges_mat)))
#PPIs_htmp <- featureHeatmap(mat = PPIs_mat,
#  colFun = PPIs_colFun,
#  datName = "PPIs",
#  rowOrder = c(1:nrow(PPIs_mat)))

#Expression_correlation_Ks_lt2_htmp <- featureHeatmap(mat = Expression_correlation_Ks_lt2_mat,
#  colFun = Expression_correlation_Ks_lt2_colFun,
#  datName = "Expr. corr. Ks<2",
#  rowOrder = c(1:nrow(Expression_correlation_Ks_lt2_mat)))
#Expression_correlation_htmp <- featureHeatmap(mat = Expression_correlation_mat,
#  colFun = Expression_correlation_colFun,
#  datName = "Expr. corr.",
#  rowOrder = c(1:nrow(Expression_correlation_mat)))


htmps <- Kappa_htmp + Stocha_htmp +
         Kappa_C_density_htmp + Stocha_C_density_htmp +
         Mean_mC_htmp + gbM_htmp +
         mD_hotspot_htmp + mD_coldspot_htmp +
         Expression_breadth_htmp + Expression_variation_htmp + Median_expression_htmp + Coexpression_module_size_htmp +
         feature_width_htmp + exons_width_prop_htmp + introns_width_prop_htmp + IRratio_median_htmp +
         Lethal_htmp + Sequence_conservation_htmp + Nucleotide_diversity_htmp

legendGap <- unit(15, "mm")

pdf(paste0(plotDir_kappa,
           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase, "_", context,
           "_NAmax", NAmax, "_all_filt_heatmap_kappa_", paste0(chrName, collapse = "_"), ".pdf"),
    width = 1.5*length(htmps), height = 10)
draw(htmps,
     gap = unit(1, "mm"),
     column_title = paste0("Genes ordered by decreasing among-read agreement (m", context, ") in gene ", featRegion), 
     column_title_gp = gpar(font = 2, fontsize = 16),
     heatmap_legend_side = "bottom",
     legend_gap = legendGap)
dev.off()

