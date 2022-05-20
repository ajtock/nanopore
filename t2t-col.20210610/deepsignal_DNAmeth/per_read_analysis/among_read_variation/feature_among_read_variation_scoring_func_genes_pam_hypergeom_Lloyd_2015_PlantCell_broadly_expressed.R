#!/usr/bin/env Rscript

# Analysis:
# 1. Lloyd et al. (2015) Plant Cell features over- and under-representation analysis of genes grouped by both among-read agreement and mean methylation proportion
# 2. Lloyd et al. (2015) Plant Cell features over- and under-representation analysis of genes grouped by both mean within-read stochasticity and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_pam_hypergeom_Lloyd_2015_PlantCell_broadly_expressed.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions'
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
library(parallel)
library(dplyr)
library(GenomicRanges)
#BiocManager::install("clusterProfiler")
#library(clusterProfiler)
#BiocManager::install("org.At.tair.db", character.only = TRUE)
#library("org.At.tair.db", character.only = TRUE)
#BiocManager::install("pathview") # installation of package ‘Rgraphviz’ had non-zero exit status; installation of package ‘pathview’ had non-zero exit status
#library(pathview)
#BiocManager::install("enrichplot")
#library(enrichplot)
#BiocManager::install("DOSE")
#library(DOSE)
library(ggplot2)

library(methods)
library(plotrix)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

#library(parallel)
#library(GenomicRanges)
#library(irr)
#library(dplyr)
#library(tidyr)
#library(cluster)
#library(fpc)
##library(data.table)
##library(segmentSeq)
#library(ComplexHeatmap)
##library(RColorBrewer)
#library(scales)
##library(circlize)
# 
#library(ggplot2)
#library(cowplot)
##library(ggcorrplot)
#library(viridis)
#library(ggthemes)
#library(tidyquant)
##library(grid)

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_Lloyd_2015_PlantCell_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_Lloyd_2015_PlantCell_", context, "_stocha_mC/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/hypergeom_Lloyd_2015_PlantCell_", context, "_kappa_stocha/")
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

## Load coordinates for mitochondrial insertion on Chr2, in BED format
#mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
#                       header = F)
#colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
#mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
#mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
#mito_ins_chr <- unique(mito_ins$chr)
#mito_ins_start <- min(mito_ins$start)+1
#mito_ins_end <- max(mito_ins$end)
##mito_ins_GR <- GRanges(seqnames = "Chr2",
##                       ranges = IRanges(start = min(mito_ins$start)+1,
##                                        end = max(mito_ins$end)),
##                       strand = "*")

# Read in Lloyd et al. (2015) Plant Cell gene features
ds3 <- read.csv(paste0("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/",
                       "TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv"),
                header = T, na.strings = c("NA", "?"))
nrow(ds3)

range(ds3$Expression.breadth, na.rm = T)
broadly_expressed <- ds3[which(ds3$Expression.breadth >= max(ds3$Expression.breadth, na.rm = T)*.75),]
nrow(broadly_expressed)


# Define hypergeometric test function
# P-value is the probability of drawing >= length(group_feat_query) [x] features
# in a sample size of length(group_feat) [k] from a total feature set consisting of
# length(genome_feat_query) [m] + ( length(genome_feat) - length(genome_feat_query) ) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(kappa_group_feat_query) to length(kappa_group_feat)

# Set class for hypergeometric test results object
setClass("hypergeomTest",
         representation(alternative = "character",
                        alpha0.05 = "numeric",
                        pval = "numeric",
                        observed = "numeric",
                        expected = "numeric",
                        log2obsexp = "numeric",
                        log2alpha = "numeric",
                        group_feat = "numeric",
                        proportion_of_group = "numeric",
                        random_proportions_of_group = "numeric",
                        hypergeomDist = "numeric"))

hgTest <- function(group, group_feat_list, genome_feat, genome_feat_query, samples_num) {

  # Get group_feat from list of groups
  group_feat <- group_feat_list[[group]]
  
  # Get intersection of gene IDs in group and gene IDs of NLRs
  group_feat_query <- intersect(group_feat, genome_feat_query)

  # Calculate the P-values for over-representation and under-representation
  # of query genes among group genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(group_feat_query):length(group_feat),
                             m = length(genome_feat_query),
                             n = length(genome_feat) - length(genome_feat_query),
                             k = length(group_feat)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(group_feat_query)-1)
  print(1 - sum(dhyper(x = 0:(length(group_feat_query)-1),
                       m = length(genome_feat_query),
                       n = length(genome_feat) - length(genome_feat_query),
                       k = length(group_feat))))

  # Under-representation
  Pval_underrep <- phyper(q = length(group_feat_query),
                          m = length(genome_feat_query),
                          n = length(genome_feat) - length(genome_feat_query),
                          k = length(group_feat))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samples_num,
                   m = length(genome_feat_query),
                   n = length(genome_feat) - length(genome_feat_query),
                   k = length(group_feat))

  # Calculate P-values and significance levels
  if(length(group_feat_query) > mean(hgDist)) {
    Pval <- Pval_overrep
    MoreOrLessThanRandom <- "MoreThanRandom"
    alpha0.05 <- quantile(hgDist, probs = 0.95)[[1]]
  } else {
    Pval <- Pval_underrep
    MoreOrLessThanRandom <- "LessThanRandom"
    alpha0.05 <- quantile(hgDist, probs = 0.05)[[1]]
  }

  hgTestResults <- new("hypergeomTest",
                       alternative = MoreOrLessThanRandom,
                       alpha0.05 = alpha0.05,
                       pval = Pval,
                       observed = length(group_feat_query),
                       expected = mean(hgDist),
                       log2obsexp = log2( (length(group_feat_query)+1) / (mean(hgDist)+1) ),
                       log2alpha  = log2( (alpha0.05+1) / (mean(hgDist)+1) ),
                       group_feat = length(group_feat),
                       proportion_of_group = length(group_feat_query) / length(group_feat),
                       random_proportions_of_group = hgDist / length(group_feat),
                       hypergeomDist = hgDist)

  data.frame(group = group,
             group_feat = hgTestResults@group_feat,
             alternative = hgTestResults@alternative,
             alpha0.05 = hgTestResults@alpha0.05,
             observed = hgTestResults@observed,
             expected = hgTestResults@expected,
             log2obsexp = hgTestResults@log2obsexp,
             log2alpha = hgTestResults@log2alpha,
             proportion_of_group = hgTestResults@proportion_of_group,
             pval = hgTestResults@pval)

}


# fk_kappa_all
# Get features that are within query feature set 
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)
featDF$parent <- sub(pattern = "\\.\\d+", replacement = "", x = featDF$name)
featDF$parent <- sub(pattern = "_\\d+", replacement = "", x = featDF$parent)

featID_broadly_expressed <- unique(featDF[which(featDF$parent %in% broadly_expressed[,1]),]$parent)
 

# Load feature groups (defined based on Fleiss' kappa vs mean mC trend plots) to enable enrichment analysis

filt_kappa_mC_groups <- lapply(seq_along(1:4), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_mC_all_cluster", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
  tmp$parent <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
  tmp$parent <- sub(pattern = "_\\d+", replacement = "", x = tmp$parent)
  tmp[ with(tmp, order(fk_kappa_all, decreasing = T)), ]
})

filt_kappa_mC_groups_featID <- lapply(seq_along(filt_kappa_mC_groups), function(x) {
  tmp <- filt_kappa_mC_groups[[x]]$fk_kappa_all
  names(tmp) <- filt_kappa_mC_groups[[x]]$parent
  names(tmp)
#  na.omit(tmp)
})

filt_kappa_mC_groups_featID_universe <- unlist(filt_kappa_mC_groups_featID)
stopifnot(length(filt_kappa_mC_groups_featID_universe) == nrow(featDF))


# Run hypergeometric test on each group of genes to evaluate
# representation of broadly_expressed genes
hgTest_kappa_broadly_expressed_list <- lapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_kappa_mC_groups_featID,
         genome_feat = filt_kappa_mC_groups_featID_universe,
         genome_feat_query = featID_broadly_expressed,
         samples_num = 1e5)
})

hgTest_kappa_broadly_expressed_DF <- dplyr::bind_rows(hgTest_kappa_broadly_expressed_list)
hgTest_kappa_broadly_expressed_DF$group <- as.character(hgTest_kappa_broadly_expressed_DF$group)

hgTest_kappa_broadly_expressed_DF$BHadj_pval <- p.adjust(hgTest_kappa_broadly_expressed_DF$pval, method = "BH")

gg_hgTest_kappa_broadly_expressed <- ggplot(data = hgTest_kappa_broadly_expressed_DF,
                                 mapping = aes(x = group,
                                               y = log2obsexp,
                                               colour = pval,
                                               size = observed)) +
  geom_point() +
  scale_colour_gradient(low = "red", high = "blue") +
  guides(colour = guide_colourbar(reverse = T)) +
  geom_hline(mapping = aes(yintercept = 0),
             linetype = "solid", size = 1, colour = "black") +
#  labs(size = "Count", colour = bquote("BH-adjusted" ~ italic(P))) +
  labs(size = "Count", colour = bquote(italic(P)*"-value")) +
  xlab(bquote(atop("Fleiss' kappa and mean m" * .(context), .(featName) ~ .(featRegion) ~ "group"))) +
  ylab(bquote(atop("Log"[2] * "(observed/expected) broadly expressed genes"))) +
  ylim(-max(abs(hgTest_kappa_broadly_expressed_DF$log2obsexp)), max(abs(hgTest_kappa_broadly_expressed_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_kappa_mC,
              sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
              featName, "_", featRegion,
              "_", paste0(chrName, collapse = "_"), "_", context,
              "_broadly_expressed_genes_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_kappa_broadly_expressed,
       height = 5, width = 6)
write.table(hgTest_kappa_broadly_expressed_DF,
            paste0(outDir,
                   sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_", paste0(chrName, collapse = "_"), "_", context,
                   "_broadly_expressed_genes_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# mean_stocha_all

# Load feature groups (defined based on stochasticity vs mean mC trend plots) to enable enrichment analysis

filt_stocha_mC_groups <- lapply(seq_along(1:4), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_mean_stocha_all_mean_mC_all_cluster", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
  tmp$parent <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
  tmp$parent <- sub(pattern = "_\\d+", replacement = "", x = tmp$parent)
  tmp[ with(tmp, order(mean_stocha_all, decreasing = T)), ]
})

filt_stocha_mC_groups_featID <- lapply(seq_along(filt_stocha_mC_groups), function(x) {
  tmp <- filt_stocha_mC_groups[[x]]$mean_stocha_all
  names(tmp) <- filt_stocha_mC_groups[[x]]$parent
  names(tmp)
#  na.omit(tmp)
})

filt_stocha_mC_groups_featID_universe <- unlist(filt_stocha_mC_groups_featID)
stopifnot(length(filt_stocha_mC_groups_featID_universe) == nrow(featDF))


# Run hypergeometric test on each group of genes to evaluate
# representation of broadly_expressed genes
hgTest_stocha_broadly_expressed_list <- lapply(seq_along(filt_stocha_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_stocha_mC_groups_featID,
         genome_feat = filt_stocha_mC_groups_featID_universe,
         genome_feat_query = featID_broadly_expressed,
         samples_num = 1e5)
})

hgTest_stocha_broadly_expressed_DF <- dplyr::bind_rows(hgTest_stocha_broadly_expressed_list)
hgTest_stocha_broadly_expressed_DF$group <- as.character(hgTest_stocha_broadly_expressed_DF$group)

hgTest_stocha_broadly_expressed_DF$BHadj_pval <- p.adjust(hgTest_stocha_broadly_expressed_DF$pval, method = "BH")

gg_hgTest_stocha_broadly_expressed <- ggplot(data = hgTest_stocha_broadly_expressed_DF,
                                mapping = aes(x = group,
                                              y = log2obsexp,
                                              colour = pval,
                                              size = observed)) +
  geom_point() +
  scale_colour_gradient(low = "red", high = "blue") +
  guides(colour = guide_colourbar(reverse = T)) +
  geom_hline(mapping = aes(yintercept = 0),
             linetype = "solid", size = 1, colour = "black") +
#  labs(size = "Count", colour = bquote("BH-adjusted" ~ italic(P))) +
  labs(size = "Count", colour = bquote(italic(P)*"-value")) +
  xlab(bquote(atop("Stochasticity and mean m" * .(context), .(featName) ~ .(featRegion) ~ "group"))) +
  ylab(bquote(atop("Log"[2] * "(observed/expected) broadly expressed genes"))) +
  ylim(-max(abs(hgTest_stocha_broadly_expressed_DF$log2obsexp)), max(abs(hgTest_stocha_broadly_expressed_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_stocha_mC,
              sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
              featName, "_", featRegion,
              "_", paste0(chrName, collapse = "_"), "_", context,
              "_broadly_expressed_genes_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_stocha_broadly_expressed,
       height = 5, width = 6)
write.table(hgTest_stocha_broadly_expressed_DF,
            paste0(outDir,
                   sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_", paste0(chrName, collapse = "_"), "_", context,
                   "_broadly_expressed_genes_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
