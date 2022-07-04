#!/usr/bin/env Rscript

# Analysis:
# Plot combined results of over- and under-representation analysis of genes grouped by both among-read agreement and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_quantiles_hypergeom_combined_bargraph_alpha.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions' 10000 1000
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
options(scipen = 999)
library(parallel)
library(dplyr)
library(ggplot2)
library(methods)
library(plotrix)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(RColorBrewer)

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
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_combined_", context, "_kappa_mC/")
plotDir_alpha_mC <- paste0(outDir, "plots/hypergeom_combined_", context, "_alpha_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_combined_", context, "_stocha_mC/")
#plotDir_alpha_stocha <- paste0(outDir, "plots/hypergeom_combined_", context, "_alpha_stocha/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_kappa_mC, " ] || mkdir -p ", plotDir_kappa_mC))
system(paste0("[ -d ", plotDir_alpha_mC, " ] || mkdir -p ", plotDir_alpha_mC))
system(paste0("[ -d ", plotDir_stocha_mC, " ] || mkdir -p ", plotDir_stocha_mC))
#system(paste0("[ -d ", plotDir_alpha_stocha, " ] || mkdir -p ", plotDir_alpha_stocha))


# Load hypergeometric test/distribution results

mD_hs <- read.table(paste0(outDir,
                           sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                           featName, "_", featRegion,
                           "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                           "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                           "_hotspot_hypergeomTest_quantiles.tsv"),
                    header = T)
mD_hs <- data.frame(Feature = rep("Epimutation hotspots", 4),
                    mD_hs)

mD_cs <- read.table(paste0(outDir,
                           sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                           featName, "_", featRegion,
                           "_mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                           "_MA1_2_MappedOn_", refbase, "_", paste0(chrName, collapse = "_"), "_", context,
                           "_coldspot_hypergeomTest_quantiles.tsv"),
                    header = T)
mD_cs <- data.frame(Feature = rep("Epimutation coldspots", 4),
                    mD_cs)

lethal <- read.table(paste0(outDir,
                            sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                            featName, "_", featRegion,
                            "_", paste0(chrName, collapse = "_"), "_", context,
                            "_lethal_genes_hypergeomTest_quantiles.tsv"),
                     header = T)
lethal <- data.frame(Feature = rep("Lethal", 4),
                     lethal)

nonlethal <- read.table(paste0(outDir,
                               sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                               featName, "_", featRegion,
                               "_", paste0(chrName, collapse = "_"), "_", context,
                               "_nonlethal_genes_hypergeomTest_quantiles.tsv"),
                        header = T)
nonlethal <- data.frame(Feature = rep("Nonlethal", 4),
                        nonlethal)

gbM <- read.table(paste0(outDir,
                         sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                         featName, "_", featRegion,
                         "_", paste0(chrName, collapse = "_"), "_", context,
                         "_gbM_genes_hypergeomTest_quantiles.tsv"),
                  header = T)
gbM <- data.frame(Feature = rep("gbM", 4),
                  gbM)

broadly_expressed <- read.table(paste0(outDir,
                                       sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                                       featName, "_", featRegion,
                                       "_", paste0(chrName, collapse = "_"), "_", context,
                                       "_broadly_expressed_genes_hypergeomTest_quantiles.tsv"),
                                header = T)
broadly_expressed <- data.frame(Feature = rep("Broadly expressed", 4),
                                broadly_expressed)

CEG <- read.table(paste0(outDir,
                         sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                         featName, "_", featRegion,
                         "_", paste0(chrName, collapse = "_"), "_", context,
                         "_CEG_genes_hypergeomTest_quantiles.tsv"),
                  header = T)
CEG <- data.frame(Feature = rep("Core eukaryotic", 4),
                  CEG)

betagammaWGD <- read.table(paste0(outDir,
                                  sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                                  featName, "_", featRegion,
                                  "_", paste0(chrName, collapse = "_"), "_", context,
                                  "_betagammaWGD_genes_hypergeomTest_quantiles.tsv"),
                           header = T)
betagammaWGD <- data.frame(Feature = rep("Beta gamma WGD", 4),
                           betagammaWGD)

alphaWGD <- read.table(paste0(outDir,
                              sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                              featName, "_", featRegion,
                              "_", paste0(chrName, collapse = "_"), "_", context,
                              "_alphaWGD_genes_hypergeomTest_quantiles.tsv"),
                       header = T)
alphaWGD <- data.frame(Feature = rep("Alpha WGD", 4),
                       alphaWGD)

tandemdup <- read.table(paste0(outDir,
                               sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
                               featName, "_", featRegion,
                               "_", paste0(chrName, collapse = "_"), "_", context,
                               "_tandemdup_genes_hypergeomTest_quantiles.tsv"),
                        header = T)
tandemdup <- data.frame(Feature = rep("Tandem duplicated", 4),
                        tandemdup)

combined <- rbind(
                  gbM, broadly_expressed,
                  mD_hs, mD_cs,
                  lethal, nonlethal,
#                  betagammaWGD, alphaWGD,
                  CEG, tandemdup
                 )

combined$group <- gsub("1", "Q1 lower", combined$group)
combined$group <- gsub("2", "Q1 upper", combined$group)
combined$group <- gsub("3", "Q2 lower", combined$group)
combined$group <- gsub("4", "Q2 upper", combined$group)
combined$group <- gsub("5", "Q3 lower", combined$group)
combined$group <- gsub("6", "Q3 upper", combined$group)
combined$group <- gsub("7", "Q4 lower", combined$group)
combined$group <- gsub("8", "Q4 upper", combined$group)

colnames(combined)[which(colnames(combined) == "group")] <- "Group"
combined$Group <- factor(combined$Group,
                         levels = sort(unique(combined$Group)))
combined$Feature <- factor(combined$Feature,
                           levels = c(
                                      "gbM",
                                      "Broadly expressed",
                                      "Epimutation hotspots", "Epimutation coldspots",
                                      "Lethal", "Nonlethal",
#                                      "Beta gamma WGD", "Alpha WGD",
                                      "Core eukaryotic", "Tandem duplicated"
                                     ))

bp <- ggplot(data = combined,
             mapping = aes(x = Group,
                           y = log2obsexp,
                           fill = Feature)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_brewer(name = "",
                    palette = "Paired") +
  geom_point(mapping = aes(x = Group,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey30", size = 7) +
  geom_segment(mapping = aes(x = 0.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 1.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[1],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 1.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 2.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[2],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 2.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 3.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[3],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 3.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 4.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[4],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 4.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 5.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[5],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 5.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 6.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[6],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 6.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 7.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[7],
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 7.55, y = min(c(log2obsexp, log2alpha))-0.50,
                             xend = 8.45, yend = min(c(log2obsexp, log2alpha))-0.50),
               colour = brewer.pal(name = "Dark2", n = 8)[8],
               inherit.aes = F, size = 2) +

  xlab(bquote(atop("Among-read agreement (alpha) and mean m" * .(context), .(featName) ~ .(featRegion) ~ "group"))) +
  ylab(bquote("Log"[2] * "(observed/expected) genes in group")) +
#  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "vertical",
                             label.position = "right",
                             label.theme = element_text(size = 16, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = length(unique(combined$Feature)),
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 20, colour = "black", hjust = 0.5)) +
  ggtitle(bquote(
                 .(prettyNum(1e5,
                             big.mark = ",",
                             trim = T)) ~ "samples from hypergeometric distribution"))
ggsave(paste0(plotDir_alpha_mC,
              sampleName, "_filt_df_ka_alpha_all_mean_mC_all_",
              featName, "_", featRegion,
              "_", paste0(chrName, collapse = "_"), "_", context,
              "_combined_bargraph_hypergeomTest_quantiles.pdf"),
       plot = bp,
       height = 8, width = 16)
