#!/usr/bin/env Rscript

# Analysis:
# 1. TE superfamily over- and under-representation analysis of TEs grouped by both among-read agreement and mean methylation proportion
# 2. TE superfamily over- and under-representation analysis of TEs grouped by both mean within-read stochasticity and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_TEs_quantiles_hypergeom_superfam.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies'
# conda deactivate
 
#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CHG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "TE"
#featRegion <- "bodies"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]

groupNames <- sort(c(paste0("quantile", 1:4, "_lower"), paste0("quantile", 1:4, "_upper")))
groupNames_alpha <- sort(c(paste0("quantile", 1:4, "_lower"), "quantile1_middle", paste0("quantile", 1:4, "_upper")))

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
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_TE_superfam_", context, "_kappa_mC/")
plotDir_alpha_mC <- paste0(outDir, "plots/hypergeom_TE_superfam_", context, "_alpha_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_TE_superfam_", context, "_stocha_mC/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/hypergeom_TE_superfam_", context, "_kappa_stocha/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_kappa_mC, " ] || mkdir -p ", plotDir_kappa_mC))
system(paste0("[ -d ", plotDir_alpha_mC, " ] || mkdir -p ", plotDir_alpha_mC))
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
  
  # Get intersection of gene IDs in group and gene IDs of query genes
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

superfamNames <- sort(unique(featDF$score))
superfamNames <- c(superfamNames[3], superfamNames[1], superfamNames[14],
                   superfamNames[10], superfamNames[7], superfamNames[13],
                   superfamNames[6], superfamNames[2], superfamNames[4],
                   superfamNames[5], superfamNames[8], superfamNames[9],
                   superfamNames[12], superfamNames[11])
superfamNamesPlot <- gsub("Pogo_Tc1_Mariner", "Pogo/Tc1/Mar", superfamNames)
superfamNamesPlot <- gsub("_", " ", superfamNamesPlot)
superfamNamesPlot <- gsub("classified", ".", superfamNamesPlot)

# Load feature groups (defined based on Fleiss' kappa vs mean mC trend plots) to enable enrichment analysis

filt_kappa_mC_groups <- lapply(groupNames, function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_mC_all_", x, "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
  tmp[ with(tmp, order(fk_kappa_all, decreasing = T)), ]
})

filt_kappa_mC_groups_featID <- lapply(seq_along(filt_kappa_mC_groups), function(x) {
  tmp <- filt_kappa_mC_groups[[x]]$fk_kappa_all
  names(tmp) <- filt_kappa_mC_groups[[x]]$name
  names(tmp)
#  na.omit(tmp)
})

filt_kappa_mC_groups_featID_universe <- unlist(filt_kappa_mC_groups_featID)
# For kappa only
featDF <- featDF[which(featDF$name %in% filt_kappa_mC_groups_featID_universe),]
#
stopifnot(length(filt_kappa_mC_groups_featID_universe) == nrow(featDF))

featID_superfamList <- lapply(superfamNames, function(y) {
  featDF[which(featDF$score == y),]$name
})


# Run hypergeometric test on each group of genes to evaluate
# representation of query genes
hgTest_kappa_superfamList_list <- lapply(1:length(superfamNames), function(y) {
  lapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
    hgTest(group = x,
           group_feat_list = filt_kappa_mC_groups_featID,
           genome_feat = filt_kappa_mC_groups_featID_universe,
           genome_feat_query = featID_superfamList[[y]],
           samples_num = 1e5)
  })
})

hgTest_kappa_superfam_DF_list <- lapply(1:length(superfamNames), function(y) {
  hgTest_kappa_superfam_DF <- dplyr::bind_rows(hgTest_kappa_superfamList_list[[y]])
  hgTest_kappa_superfam_DF$group <- as.character(hgTest_kappa_superfam_DF$group)
  hgTest_kappa_superfam_DF$BHadj_pval <- p.adjust(hgTest_kappa_superfam_DF$pval, method = "BH")
  hgTest_kappa_superfam_DF
})

for(y in 1:length(superfamNames)) {
  gg_hgTest_kappa_superfam <- ggplot(data = hgTest_kappa_superfam_DF_list[[y]],
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
    ylab(bquote(atop("Log"[2] * "(observed/expected)" ~ .(superfamNamesPlot[y]) ~ "TEs"))) +
    ylim(-max(abs(hgTest_kappa_superfam_DF_list[[y]]$log2obsexp)), max(abs(hgTest_kappa_superfam_DF_list[[y]]$log2obsexp))) +
    theme_bw()
  ggsave(paste0(plotDir_kappa_mC,
                sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                featName, "_", featRegion,
                "_", paste0(chrName, collapse = "_"), "_", context,
                "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.pdf"),
         plot = gg_hgTest_kappa_superfam,
         height = 5, width = 6)
  write.table(hgTest_kappa_superfam_DF_list[[y]],
              paste0(outDir,
                     sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                     featName, "_", featRegion,
                     "_", paste0(chrName, collapse = "_"), "_", context,
                     "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
}


# ka_alpha_all

# Get features that are within query feature set 
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)

superfamNames <- sort(unique(featDF$score))
superfamNames <- c(superfamNames[3], superfamNames[1], superfamNames[14],
                   superfamNames[10], superfamNames[7], superfamNames[13],
                   superfamNames[6], superfamNames[2], superfamNames[4],
                   superfamNames[5], superfamNames[8], superfamNames[9],
                   superfamNames[12], superfamNames[11])
superfamNamesPlot <- gsub("Pogo_Tc1_Mariner", "Pogo/Tc1/Mar", superfamNames)
superfamNamesPlot <- gsub("_", " ", superfamNamesPlot)
superfamNamesPlot <- gsub("classified", ".", superfamNamesPlot)

# Load feature groups (defined based on Krippendorff's alpha vs mean mC trend plots) to enable enrichment analysis

filt_alpha_mC_groups <- lapply(groupNames_alpha, function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_ka_alpha_all_mean_mC_all_", x, "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
  tmp[ with(tmp, order(ka_alpha_all, decreasing = T)), ]
})

filt_alpha_mC_groups_featID <- lapply(seq_along(filt_alpha_mC_groups), function(x) {
  tmp <- filt_alpha_mC_groups[[x]]$ka_alpha_all
  names(tmp) <- filt_alpha_mC_groups[[x]]$name
  names(tmp)
#  na.omit(tmp)
})

filt_alpha_mC_groups_featID_universe <- unlist(filt_alpha_mC_groups_featID)
# For alpha only
featDF <- featDF[which(featDF$name %in% filt_alpha_mC_groups_featID_universe),]
#
stopifnot(length(filt_alpha_mC_groups_featID_universe) == nrow(featDF))

featID_superfamList <- lapply(superfamNames, function(y) {
  featDF[which(featDF$score == y),]$name
})


# Run hypergeometric test on each group of genes to evaluate
# representation of query genes
hgTest_alpha_superfamList_list <- lapply(1:length(superfamNames), function(y) {
  lapply(seq_along(filt_alpha_mC_groups_featID), function(x) {
    hgTest(group = x,
           group_feat_list = filt_alpha_mC_groups_featID,
           genome_feat = filt_alpha_mC_groups_featID_universe,
           genome_feat_query = featID_superfamList[[y]],
           samples_num = 1e5)
  })
})

hgTest_alpha_superfam_DF_list <- lapply(1:length(superfamNames), function(y) {
  hgTest_alpha_superfam_DF <- dplyr::bind_rows(hgTest_alpha_superfamList_list[[y]])
  hgTest_alpha_superfam_DF$group <- as.character(hgTest_alpha_superfam_DF$group)
  hgTest_alpha_superfam_DF$BHadj_pval <- p.adjust(hgTest_alpha_superfam_DF$pval, method = "BH")
  hgTest_alpha_superfam_DF
})

for(y in 1:length(superfamNames)) {
  gg_hgTest_alpha_superfam <- ggplot(data = hgTest_alpha_superfam_DF_list[[y]],
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
    xlab(bquote(atop("Krippendorff's alpha and mean m" * .(context), .(featName) ~ .(featRegion) ~ "group"))) +
    ylab(bquote(atop("Log"[2] * "(observed/expected)" ~ .(superfamNamesPlot[y]) ~ "TEs"))) +
    ylim(-max(abs(hgTest_alpha_superfam_DF_list[[y]]$log2obsexp)), max(abs(hgTest_alpha_superfam_DF_list[[y]]$log2obsexp))) +
    theme_bw()
  ggsave(paste0(plotDir_alpha_mC,
                sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                featName, "_", featRegion,
                "_", paste0(chrName, collapse = "_"), "_", context,
                "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.pdf"),
         plot = gg_hgTest_alpha_superfam,
         height = 5, width = 6)
  write.table(hgTest_alpha_superfam_DF_list[[y]],
              paste0(outDir,
                     sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                     featName, "_", featRegion,
                     "_", paste0(chrName, collapse = "_"), "_", context,
                     "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
}


# mean_stocha_all

# Get features that are within query feature set 
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)

superfamNames <- sort(unique(featDF$score))
superfamNames <- c(superfamNames[3], superfamNames[1], superfamNames[14],
                   superfamNames[10], superfamNames[7], superfamNames[13],
                   superfamNames[6], superfamNames[2], superfamNames[4],
                   superfamNames[5], superfamNames[8], superfamNames[9],
                   superfamNames[12], superfamNames[11])
superfamNamesPlot <- gsub("PogoTc1Mariner", "Pogo/Tc1/Mar", superfamNames)
superfamNamesPlot <- gsub("_", " ", superfamNamesPlot)
superfamNamesPlot <- gsub("classified", ".", superfamNamesPlot)

# Load feature groups (defined based on stochasticity vs mean mC trend plots) to enable enrichment analysis

filt_stocha_mC_groups <- lapply(groupNames, function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_mean_stocha_all_mean_mC_all_", x, "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
  tmp[ with(tmp, order(mean_stocha_all, decreasing = T)), ]
})

filt_stocha_mC_groups_featID <- lapply(seq_along(filt_stocha_mC_groups), function(x) {
  tmp <- filt_stocha_mC_groups[[x]]$mean_stocha_all
  names(tmp) <- filt_stocha_mC_groups[[x]]$name
  names(tmp)
#  na.omit(tmp)
})

filt_stocha_mC_groups_featID_universe <- unlist(filt_stocha_mC_groups_featID)
## Not required for stocha
#featDF <- featDF[which(featDF$name %in% filt_stocha_mC_groups_featID_universe),]
##
stopifnot(length(filt_stocha_mC_groups_featID_universe) == nrow(featDF))

featID_superfamList <- lapply(superfamNames, function(y) {
  featDF[which(featDF$score == y),]$name
})


# Run hypergeometric test on each group of genes to evaluate
# representation of query genes
hgTest_stocha_superfamList_list <- lapply(1:length(superfamNames), function(y) {
  lapply(seq_along(filt_stocha_mC_groups_featID), function(x) {
    hgTest(group = x,
           group_feat_list = filt_stocha_mC_groups_featID,
           genome_feat = filt_stocha_mC_groups_featID_universe,
           genome_feat_query = featID_superfamList[[y]],
           samples_num = 1e5)
  })
})

hgTest_stocha_superfam_DF_list <- lapply(1:length(superfamNames), function(y) {
  hgTest_stocha_superfam_DF <- dplyr::bind_rows(hgTest_stocha_superfamList_list[[y]])
  hgTest_stocha_superfam_DF$group <- as.character(hgTest_stocha_superfam_DF$group)
  hgTest_stocha_superfam_DF$BHadj_pval <- p.adjust(hgTest_stocha_superfam_DF$pval, method = "BH")
  hgTest_stocha_superfam_DF
})

for(y in 1:length(superfamNames)) {
  gg_hgTest_stocha_superfam <- ggplot(data = hgTest_stocha_superfam_DF_list[[y]],
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
    ylab(bquote(atop("Log"[2] * "(observed/expected)" ~ .(superfamNamesPlot[y]) ~ "TEs"))) +
    ylim(-max(abs(hgTest_stocha_superfam_DF_list[[y]]$log2obsexp)), max(abs(hgTest_stocha_superfam_DF_list[[y]]$log2obsexp))) +
    theme_bw()
  ggsave(paste0(plotDir_stocha_mC,
                sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
                featName, "_", featRegion,
                "_", paste0(chrName, collapse = "_"), "_", context,
                "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.pdf"),
         plot = gg_hgTest_stocha_superfam,
         height = 5, width = 6)
  write.table(hgTest_stocha_superfam_DF_list[[y]],
              paste0(outDir,
                     sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
                     featName, "_", featRegion,
                     "_", paste0(chrName, collapse = "_"), "_", context,
                     "_", superfamNames[y], "_TEs_hypergeomTest_quantiles.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
}
