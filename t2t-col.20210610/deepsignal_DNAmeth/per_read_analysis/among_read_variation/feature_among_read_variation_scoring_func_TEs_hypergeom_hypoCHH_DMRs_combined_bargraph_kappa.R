#!/usr/bin/env Rscript

# Analysis:
# Plot combined results of DMR-overlapping TE over- and under-representation analysis of TEs grouped by both among-read agreement and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_TEs_hypergeom_hypoCHH_DMRs_combined_bargraph_kappa.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies' 'drd1_BSseq_Rep1_hypoCHH,cmt2_BSseq_Rep1_hypoCHH,cmt3_BSseq_Rep1_hypoCHH,kss_BSseq_Rep1_hypoCHH,suvh4_BSseq_Rep1_hypoCHH,suvh5_BSseq_Rep1_hypoCHH,suvh6_BSseq_Rep1_hypoCHH,met1_BSseq_Rep1_hypoCHH,ddm1_BSseq_Rep1_hypoCHH'
# ./feature_among_read_variation_scoring_func_TEs_hypergeom_hypoCHH_DMRs_combined_bargraph_kappa.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies' 'drd1_BSseq_Rep1_hypoCHH,drm1_drm2_BSseq_Rep1_hypoCHH,cmt2_BSseq_Rep1_hypoCHH,cmt3_BSseq_Rep1_hypoCHH,cmt2_cmt3_BSseq_Rep1_hypoCHH,kss_BSseq_Rep1_hypoCHH,suvh4_BSseq_Rep1_hypoCHH,suvh5_BSseq_Rep1_hypoCHH,suvh6_BSseq_Rep1_hypoCHH,met1_BSseq_Rep1_hypoCHH,met1_cmt3_BSseq_Rep1_hypoCHH,met1het_BSseq_Rep1_hypoCHH,ddm1_BSseq_Rep1_hypoCHH'
# conda deactivate
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CHG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "TE"
#featRegion <- "bodies"
#DMRsNames <- unlist(strsplit("kss_BSseq_Rep1_hypoCHH,cmt3_BSseq_Rep1_hypoCHH,kss_BSseq_Rep1_hypoCHH,cmt3_BSseq_Rep1_hypoCHH",
#                              split = ",")) 
#DMRsNamesPlot <- sub("_BSseq_Rep1_", " ", DMRsNames)

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]
DMRsNames <- unlist(strsplit(args[8], split = ","))
DMRsNamesPlot <- sub("_BSseq_Rep1_", " ", DMRsNames)

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
library(pals)

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_combined_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_combined_", context, "_stocha_mC/")
#plotDir_kappa_stocha <- paste0(outDir, "plots/hypergeom_combined_", context, "_kappa_stocha/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_kappa_mC, " ] || mkdir -p ", plotDir_kappa_mC))
system(paste0("[ -d ", plotDir_stocha_mC, " ] || mkdir -p ", plotDir_stocha_mC))
#system(paste0("[ -d ", plotDir_kappa_stocha, " ] || mkdir -p ", plotDir_kappa_stocha))

score_colFun <- cols25(n = 25)[-c(7:16, 25)][1:length(DMRsNamesPlot)]
stopifnot(length(score_colFun) == length(DMRsNamesPlot))
#names(score_colFun) <- DMRsNamesPlot

DMRsList <- lapply(1:length(DMRsNames), function(y) {
  DMRs <-  read.table(paste0(outDir,
                             sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
                             featName, "_", featRegion,
                             "_", paste0(chrName, collapse = "_"), "_", context,
                             "_", DMRsNames[y], "_TEs_hypergeomTest.tsv"),
                      header = T)
  DMRs <- data.frame(Feature = rep(DMRsNamesPlot[y], 8),
                     DMRs)
  DMRs
})

combined <- dplyr::bind_rows(DMRsList)

combined$group <- paste0("Group ", combined$group)
colnames(combined)[which(colnames(combined) == "group")] <- "Group"
combined$Group <- factor(combined$Group,
                         levels = sort(unique(combined$Group)))
combined$Feature <- factor(combined$Feature,
                           levels = c(
                                      DMRsNamesPlot
                                     ))

bp <- ggplot(data = combined,
             mapping = aes(x = Group,
                           y = log2obsexp,
                           fill = Feature)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = score_colFun) +
  geom_point(mapping = aes(x = Group,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 6) +
  geom_segment(mapping = aes(x = 0.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 1.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "dodgerblue1",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 1.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 2.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "dodgerblue3",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 2.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 3.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "magenta1",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 3.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 4.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "magenta3",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 4.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 5.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "darkorange1",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 5.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 6.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "darkorange3",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 6.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 7.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "firebrick1",
               inherit.aes = F, size = 2) +
  geom_segment(mapping = aes(x = 7.55, y = min(c(combined$log2obsexp, combined$log2alpha))-0.50,
                             xend = 8.45, yend = min(c(combined$log2obsexp, combined$log2alpha))-0.50),
               colour = "firebrick3",
               inherit.aes = F, size = 2) +

  xlab(bquote(atop("Among-read agreement and mean m" * .(context), .(featName) ~ .(featRegion) ~ "group"))) +
  ylab(bquote("Log"[2] * "(observed/expected) TEs in group")) +
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
ggsave(paste0(plotDir_kappa_mC,
              sampleName, "_filt_df_fk_kappa_all_mean_mC_all_",
              featName, "_", featRegion,
              "_", paste0(chrName, collapse = "_"), "_", context,
              "_TE_hypoCHH_DMRs_combined_bargraph_hypergeomTest.pdf"),
       plot = bp,
       height = 8, width = 20)
