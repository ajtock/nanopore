#!/usr/bin/env Rscript

# Analysis:
# 1. Epimutation hotspot over- and under-representation analysis of genes grouped by both among-read agreement and mean methylation proportion
# 1. Epimutation hotspot over- and under-representation analysis of genes grouped by both mean within-read stochasticity and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_pam_hypergeom_epimutation_hotspots.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions' 10000 1000
# conda deactivate
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "regions"
#genomeBinSize <- 10000
#genomeStepSize <- 1000

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]
genomeBinSize <- as.numeric(args[8])
genomeStepSize <- as.numeric(args[9])

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

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
  genomeBinNamePlot <- paste0(genomeBinSize, "-bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e3, "-kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e6, "-Mb")
}

if(floor(log10(genomeStepSize)) + 1 < 4) {
  genomeStepName <- paste0(genomeStepSize, "bp")
  genomeStepNamePlot <- paste0(genomeStepSize, "-bp")
} else if(floor(log10(genomeStepSize)) + 1 >= 4 &
          floor(log10(genomeStepSize)) + 1 <= 6) {
  genomeStepName <- paste0(genomeStepSize/1e3, "kb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e3, "-kb")
} else if(floor(log10(genomeStepSize)) + 1 >= 7) {
  genomeStepName <- paste0(genomeStepSize/1e6, "Mb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e6, "-Mb")
}

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_stocha_mC/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_kappa_stocha/")
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

# Read in methylation divergence (mD) files
mDdir <- "/home/ajt200/analysis/BSseq_leaf_Hazarika_Johannes_2022_NatPlants/snakemake_BSseq_t2t-col.20210610/MA1_2/coverage/report/alphabeta/"
MA1_2_mD_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0(mDdir,
                    "mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
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

# Histogram of mD values
gg_dens_mD <- ggplot(tab_mD,
                     mapping = aes(x = MA1_2_mean.D)) +
  geom_histogram(binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), colour = "black", fill = "white") +
  geom_vline(mapping = aes(xintercept = min(mD_hs$MA1_2_mean.D)),
             linetype = "dashed", size = 0.5, colour = "red") +
  geom_vline(mapping = aes(xintercept = min(mD_cs$MA1_2_mean.D)),
             linetype = "dashed", size = 0.5, colour = "dodgerblue") +
  xlab(bquote("m" * .(context) ~ "divergence (" * italic(D) * ")")) +
  theme_bw()
ggsave(paste0(plotDir,
              "mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context, "_density.pdf"),
       plot = gg_dens_mD,
       height = 5, width = 4)


# Define hypergeometric test function
# P-value is the probability of drawing >= length(group_feat_mD_hs) [x] features
# in a sample size of length(group_feat) [k] from a total feature set consisting of
# length(genome_feat_mD_hs) [m] + ( length(genome_feat) - length(genome_feat_mD_hs) ) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(kappa_group_feat_mD_hs) to length(kappa_group_feat)

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

hgTest <- function(group, group_feat_list, genome_feat, genome_feat_mD_loci, samples_num) {

  # Get group_feat from list of groups
  group_feat <- group_feat_list[[group]]
  
  # Get intersection of gene IDs in group and gene IDs of NLRs
  group_feat_mD_loci <- intersect(group_feat, genome_feat_mD_loci)

  # Calculate the P-values for over-representation and under-representation
  # of NLRs among group genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(group_feat_mD_loci):length(group_feat),
                             m = length(genome_feat_mD_loci),
                             n = length(genome_feat) - length(genome_feat_mD_loci),
                             k = length(group_feat)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(group_feat_mD_loci)-1)
  print(1 - sum(dhyper(x = 0:(length(group_feat_mD_loci)-1),
                       m = length(genome_feat_mD_loci),
                       n = length(genome_feat) - length(genome_feat_mD_loci),
                       k = length(group_feat))))

  # Under-representation
  Pval_underrep <- phyper(q = length(group_feat_mD_loci),
                          m = length(genome_feat_mD_loci),
                          n = length(genome_feat) - length(genome_feat_mD_loci),
                          k = length(group_feat))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samples_num,
                   m = length(genome_feat_mD_loci),
                   n = length(genome_feat) - length(genome_feat_mD_loci),
                   k = length(group_feat))

  # Calculate P-values and significance levels
  if(length(group_feat_mD_loci) > mean(hgDist)) {
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
                       observed = length(group_feat_mD_loci),
                       expected = mean(hgDist),
                       log2obsexp = log2( (length(group_feat_mD_loci)+1) / (mean(hgDist)+1) ),
                       log2alpha  = log2( (alpha0.05+1) / (mean(hgDist)+1) ),
                       group_feat = length(group_feat),
                       proportion_of_group = length(group_feat_mD_loci) / length(group_feat),
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
# Get features that overlap epimutation hotspots and coldspots
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)

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
#featID_mD_hs <- sub(pattern = "\\.\\d+", replacement = "", x = featID_mD_hs)
#featID_mD_hs <- sub(pattern = "_\\d+", replacement = "", x = featID_mD_hs)
#featID_mD_hs <- unique(featID_mD_hs)

fOverlaps_feat_mD_cs <- findOverlaps(query = featGR,
                                     subject = mD_cs_GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = T)
featGR_mD_cs <- featGR[unique(queryHits(fOverlaps_feat_mD_cs))]
featID_mD_cs <- unique(featGR_mD_cs$featID)
#featID_mD_cs <- sub(pattern = "\\.\\d+", replacement = "", x = featID_mD_cs)
#featID_mD_cs <- sub(pattern = "_\\d+", replacement = "", x = featID_mD_cs)
#featID_mD_cs <- unique(featID_mD_cs)
 

# Load feature groups (defined based on Fleiss' kappa vs mean mC trend plots) to enable enrichment analysis

filt_kappa_mC_groups <- lapply(seq_along(1:4), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_mC_all_cluster", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
#  tmp$name <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
#  tmp$name <- sub(pattern = "_\\d+", replacement = "", x = tmp$name)
  tmp[ with(tmp, order(fk_kappa_all, decreasing = T)), ]
})

filt_kappa_mC_groups_featID <- lapply(seq_along(filt_kappa_mC_groups), function(x) {
  tmp <- filt_kappa_mC_groups[[x]]$fk_kappa_all
  names(tmp) <- filt_kappa_mC_groups[[x]]$name
  names(tmp)
#  na.omit(tmp)
})

filt_kappa_mC_groups_featID_universe <- unlist(filt_kappa_mC_groups_featID)
stopifnot(length(filt_kappa_mC_groups_featID_universe) == length(featGR))


# Run hypergeometric test on each group of genes to evaluate
# representation of epimutation hotspots
hgTest_kappa_mD_hs_list <- lapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_kappa_mC_groups_featID,
         genome_feat = filt_kappa_mC_groups_featID_universe,
         genome_feat_mD_loci = featID_mD_hs,
         samples_num = 1e5)
})

hgTest_kappa_mD_hs_DF <- dplyr::bind_rows(hgTest_kappa_mD_hs_list)
hgTest_kappa_mD_hs_DF$group <- as.character(hgTest_kappa_mD_hs_DF$group)

hgTest_kappa_mD_hs_DF$BHadj_pval <- p.adjust(hgTest_kappa_mD_hs_DF$pval, method = "BH")

gg_hgTest_kappa_mD_hs <- ggplot(data = hgTest_kappa_mD_hs_DF,
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
  ylab(bquote(atop("Log"[2] * "(observed/expected)", "m" * .(context) ~ "divergence (" * italic(D) * ") hotspot overlap"))) +
  ylim(-max(abs(hgTest_kappa_mD_hs_DF$log2obsexp)), max(abs(hgTest_kappa_mD_hs_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_kappa_mC,
              sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
              featName, "_", featRegion,
              "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
              "_hotspot_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_kappa_mD_hs,
       height = 5, width = 6)
write.table(hgTest_kappa_mD_hs_DF,
            paste0(outDir,
                   sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                   "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                   "_hotspot_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Run hypergeometric test on each group of genes to evaluate
# representation of epimutation coldspots
hgTest_kappa_mD_cs_list <- lapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_kappa_mC_groups_featID,
         genome_feat = filt_kappa_mC_groups_featID_universe,
         genome_feat_mD_loci = featID_mD_cs,
         samples_num = 1e5)
})

hgTest_kappa_mD_cs_DF <- dplyr::bind_rows(hgTest_kappa_mD_cs_list)
hgTest_kappa_mD_cs_DF$group <- as.character(hgTest_kappa_mD_cs_DF$group)

hgTest_kappa_mD_cs_DF$BHadj_pval <- p.adjust(hgTest_kappa_mD_cs_DF$pval, method = "BH")

gg_hgTest_kappa_mD_cs <- ggplot(data = hgTest_kappa_mD_cs_DF,
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
  ylab(bquote(atop("Log"[2] * "(observed/expected)", "m" * .(context) ~ "divergence (" * italic(D) * ") coldspot overlap"))) +
  ylim(-max(abs(hgTest_kappa_mD_cs_DF$log2obsexp)), max(abs(hgTest_kappa_mD_cs_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_kappa_mC,
              sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
              featName, "_", featRegion,
              "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
              "_coldspot_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_kappa_mD_cs,
       height = 5, width = 6)
write.table(hgTest_kappa_mD_cs_DF,
            paste0(outDir,
                   sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                   "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                   "_coldspot_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# mean_stocha_all
# Get features that overlap epimutation hotspots and coldspots

# Load feature groups (defined based on stochasticity vs mean mC trend plots) to enable enrichment analysis

filt_stocha_mC_groups <- lapply(seq_along(1:4), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_mean_stocha_all_mean_mC_all_cluster", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T, sep = "\t")
#  tmp$name <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
#  tmp$name <- sub(pattern = "_\\d+", replacement = "", x = tmp$name)
  tmp[ with(tmp, order(mean_stocha_all, decreasing = T)), ]
})

filt_stocha_mC_groups_featID <- lapply(seq_along(filt_stocha_mC_groups), function(x) {
  tmp <- filt_stocha_mC_groups[[x]]$mean_stocha_all
  names(tmp) <- filt_stocha_mC_groups[[x]]$name
  names(tmp)
#  na.omit(tmp)
})

filt_stocha_mC_groups_featID_universe <- unlist(filt_stocha_mC_groups_featID)
stopifnot(length(filt_stocha_mC_groups_featID_universe) == length(featGR))


# Run hypergeometric test on each group of genes to evaluate
# representation of epimutation hotspots
hgTest_stocha_mD_hs_list <- lapply(seq_along(filt_stocha_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_stocha_mC_groups_featID,
         genome_feat = filt_stocha_mC_groups_featID_universe,
         genome_feat_mD_loci = featID_mD_hs,
         samples_num = 1e5)
})

hgTest_stocha_mD_hs_DF <- dplyr::bind_rows(hgTest_stocha_mD_hs_list)
hgTest_stocha_mD_hs_DF$group <- as.character(hgTest_stocha_mD_hs_DF$group)

hgTest_stocha_mD_hs_DF$BHadj_pval <- p.adjust(hgTest_stocha_mD_hs_DF$pval, method = "BH")

gg_hgTest_stocha_mD_hs <- ggplot(data = hgTest_stocha_mD_hs_DF,
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
  ylab(bquote(atop("Log"[2] * "(observed/expected)", "m" * .(context) ~ "divergence (" * italic(D) * ") hotspot overlap"))) +
  ylim(-max(abs(hgTest_stocha_mD_hs_DF$log2obsexp)), max(abs(hgTest_stocha_mD_hs_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_stocha_mC,
              sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
              featName, "_", featRegion,
              "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
              "_hotspot_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_stocha_mD_hs,
       height = 5, width = 6)
write.table(hgTest_stocha_mD_hs_DF,
            paste0(outDir,
                   sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                   "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                   "_hotspot_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Run hypergeometric test on each group of genes to evaluate
# representation of epimutation coldspots
hgTest_stocha_mD_cs_list <- lapply(seq_along(filt_stocha_mC_groups_featID), function(x) {
  hgTest(group = x,
         group_feat_list = filt_stocha_mC_groups_featID,
         genome_feat = filt_stocha_mC_groups_featID_universe,
         genome_feat_mD_loci = featID_mD_cs,
         samples_num = 1e5)
})

hgTest_stocha_mD_cs_DF <- dplyr::bind_rows(hgTest_stocha_mD_cs_list)
hgTest_stocha_mD_cs_DF$group <- as.character(hgTest_stocha_mD_cs_DF$group)

hgTest_stocha_mD_cs_DF$BHadj_pval <- p.adjust(hgTest_stocha_mD_cs_DF$pval, method = "BH")

gg_hgTest_stocha_mD_cs <- ggplot(data = hgTest_stocha_mD_cs_DF,
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
  ylab(bquote(atop("Log"[2] * "(observed/expected)", "m" * .(context) ~ "divergence (" * italic(D) * ") coldspot overlap"))) +
  ylim(-max(abs(hgTest_stocha_mD_cs_DF$log2obsexp)), max(abs(hgTest_stocha_mD_cs_DF$log2obsexp))) +
  theme_bw()
ggsave(paste0(plotDir_stocha_mC,
              sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
              featName, "_", featRegion,
              "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
              "_coldspot_hypergeomTest_clusters.pdf"),
       plot = gg_hgTest_stocha_mD_cs,
       height = 5, width = 6)
write.table(hgTest_stocha_mD_cs_DF,
            paste0(outDir,
                   sampleName, "_filt_df_mean_stocha_all_mean_mC_all_",
                   featName, "_", featRegion,
                   "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                   "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                   "_coldspot_hypergeomTest_clusters.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
