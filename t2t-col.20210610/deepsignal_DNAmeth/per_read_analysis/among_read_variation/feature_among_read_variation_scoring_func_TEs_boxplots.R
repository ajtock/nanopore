#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Plot among-read variation/agreement (e.g., Fleiss' kappa) and stochasticity for each TE superfamily (e.g., as boxplots or violin plots)

# Usage:
# /applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_TEs_boxplots.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies'
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "TE"
#featRegion <- "bodies"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]
featRegion <- args[8]

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
#library(GenomicRanges)
#library(irr)
library(dplyr)
#library(tidyr)
#library(cluster)
#library(fpc)
##library(data.table)
##library(segmentSeq)
#library(ComplexHeatmap)
##library(RColorBrewer)
library(scales)
##library(circlize)
 
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggthemes)
##library(ggcorrplot)
library(viridis)
#library(ggthemes)
#library(tidyquant)
##library(grid)

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Load among-read and within-read mC data for featName featRegion
con_fk_df_all <- read.table(paste0(outDir,
                                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                                   "_", context,
                                   "_NAmax", NAmax,
                                   "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                                   paste0(chrName, collapse = "_"), ".tsv"),
                            header = T)

con_fk_df_all_filt <- read.table(paste0(outDir,
                                        featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                                        "_", context,
                                        "_NAmax", NAmax,
                                        "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                                        paste0(chrName, collapse = "_"), ".tsv"),
                                 header = T)

dataFrame = con_fk_df_all_filt
mapping = 
violinPlotFun <- function(dataFrame) {
  vplot <- ggplot(data = dataFrame,
                  mapping = aes(x = score,
                                y = metric 
## RF plotting function
RFplotFun <- function(dataFrame, intervalName, Pval, genoColours) {
  RFplot <- ggplot(data = dataFrame[dataFrame$interval == intervalName,],
                   mapping = aes(x = genotype,
                                 y = RF,
                                 colour = genotype)) +
            scale_colour_manual(values = genoColours) +
            #geom_boxplot(mapping = aes(colour = genotype),
            #             varwidth = TRUE) +
            geom_violin(scale = "count",
                        trim = FALSE,
                        draw_quantiles = c(0.25, 0.50, 0.75)) +
            geom_beeswarm(cex = 6,
                          size = 4) +
            labs(x = "",
                 y = "Recombination frequency (cM)") +
            theme_bw() +
            theme(axis.line.y = element_line(size = 1.0, colour = "black"),
                  axis.ticks.y = element_line(size = 1.0, colour = "black"),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(size = 18, colour = "black"),
                  axis.text.x = element_text(size = 18, colour = genoColours),
                  axis.title = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"),
                  plot.title = element_text(hjust = 0.5, size = 30)) +
            ggtitle(bquote(italic(.(intervalName)) ~ ~ ~ ~
                           "MWW" ~ italic("P") ~ .(Pval)))
  ggsave(RFplot,
         file = paste0("./", intervalName, "_",
                       geno1Name, "_vs_", geno2Name,
                       "_RF_Utest_and_plot.pdf"),
         width = 18, height = 18, units = "cm")
}


  # Plot orderingFactor means and LSDs for IDs vs annoGOIDs
  popgen_stats_meanLSDs <- function(dataFrame,
                                    parameterLab,
                                    featureGroup,
                                    featureNamePlot) {
    ggplot(data = dataFrame,
           mapping = aes(x = get(featureGroup),
                         y = mean,
                         colour = get(featureGroup))) +
    labs(colour = "") +
    geom_point(shape = "-", size = 18, position = position_dodge(width = 0.2), colour = "black") +
    geom_errorbar(mapping = aes(ymin = mean-(lsd/2),
                                ymax = mean+(lsd/2)),
                  width = 0.5, size = 1, position = position_dodge(width = 0.2), colour = "black") +
    geom_beeswarm(data = IDsDF_annoGOIDsDF,
                  mapping = aes(x = get(featureGroup),
                                y = get(orderingFactor),
                                colour = get(featureGroup)),
                  priority = "ascending",
                  cex = 1,
                  size = 1) +
    scale_colour_manual(values = quantileColours) +
    scale_y_continuous(
#                       limits = c(summary_stats_min, summary_stats_max),
                       labels = function(x) sprintf(yDec, x)) +
#    scale_x_discrete(breaks = as.vector(dataFrame$quantile),
#                     labels = as.vector(dataFrame$quantile)) +
    labs(x = "",
         y = parameterLab) +
    theme_bw() +
    theme(axis.line.y = element_line(size = 2.0, colour = "black"),
          axis.ticks.y = element_line(size = 2.0, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.text.x = element_text(size = 22, colour = quantileColours, hjust = 1.0, vjust = 1.0, angle = 45),
          axis.title = element_text(size = 26, colour = "black"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0.8,1.2,0.1,0.3),"cm"),
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote(atop(.(featureNamePlot),
                        atop(italic("t") * "-test" ~ italic("P") ~ .(ttestPvalChar),
                             "Yuen's test (10% trimmed mean)"  ~ italic("P") ~ .(yuenttestPvalChar)))))
#                             "MWW test"  ~ italic("P") ~ .(UtestPvalChar)))))
  }


  ggObjGA_feature_mean <- popgen_stats_meanLSDs(dataFrame = estimates,
                                                parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                                featureGroup = "quantile",
                                                featureNamePlot = gsub("_", " ", featureNamePlot)
                                               )



# Plot relationships and define groups
violinPlot <- function(dataFrame, mapping, xvar, yvar, xlab, ylab, yaxtrans, ybreaks, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = mapping,
         colour = score) +
#  scale_colour_manual(values = superfamColours) +
#  geom_boxplot(mapping = aes(colour = genotype),
#               varwidth = TRUE) +
  geom_violin(scale = "count",
              trim = FALSE,
              draw_quantiles = c(0.25, 0.50, 0.75)) +
  geom_beeswarm(cex = 6,
                size = 4) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  labs(x = xlab,
       y = ylab) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 1.5, colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        panel.grid = element_blank(),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18))
}

if(featRegion %in% c("bodies", "regions")) {
  if(context == "CpG") {
    fk_kappa_all_high <- 0.55
    fk_kappa_all_mid  <- 0.35
    fk_kappa_all_low  <- 0.04
    mean_stocha_all_high <- 0.28
    mean_stocha_all_mid  <- 0.17
    mean_stocha_all_low  <- 0.08
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.75
    mean_mC_all_mid   <- 0.25
    mean_mC_all_low   <- 0.10
  } else if(context == "CHG") {
    fk_kappa_all_high <- 0.05623413
    fk_kappa_all_mid  <- 0.01778279 
    fk_kappa_all_low  <- 0.01000000 
    mean_stocha_all_high <- 0.28183830
    mean_stocha_all_mid  <- 0.05011872
    mean_stocha_all_low  <- 0.03162278
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.05623413
    mean_mC_all_mid   <- 0.02511886
    mean_mC_all_low   <- 0.01778279
  } else if(context == "CHH") {
    fk_kappa_all_high <- 0.05623413
    fk_kappa_all_mid  <- 0.02511886
    fk_kappa_all_low  <- 0.01778279 
    mean_stocha_all_high <- 0.03162278
    mean_stocha_all_mid  <- 0.01737801
    mean_stocha_all_low  <- 0.01000000
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.02500000
    mean_mC_all_mid   <- 0.01737801
    mean_mC_all_low   <- 0.01000000
  }
} else if(featRegion %in% c("promoters", "terminators")) {
  if(context == "CpG") {
    fk_kappa_all_high <- 0.55
    fk_kappa_all_mid  <- 0.35
    fk_kappa_all_low  <- 0.04
    mean_stocha_all_high <- 0.22
    mean_stocha_all_mid  <- 0.17
    mean_stocha_all_low  <- 0.08
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.75
    mean_mC_all_mid   <- 0.25
    mean_mC_all_low   <- 0.10
  } else if(context == "CHG") {
    fk_kappa_all_high <- 0.05623413
    fk_kappa_all_mid  <- 0.01778279 
    fk_kappa_all_low  <- 0.01000000 
    mean_stocha_all_high <- 0.22387210
    mean_stocha_all_mid  <- 0.05011872
    mean_stocha_all_low  <- 0.03162278
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.05623413
    mean_mC_all_mid   <- 0.02511886
    mean_mC_all_low   <- 0.01778279
  } else if(context == "CHH") {
    fk_kappa_all_high <- 0.05623413
    fk_kappa_all_mid  <- 0.02511886
    fk_kappa_all_low  <- 0.01778279 
    mean_stocha_all_high <- 0.03162278
    mean_stocha_all_mid  <- 0.01737801
    mean_stocha_all_low  <- 0.01000000
    mean_min_acf_all_high <- -0.05
    mean_min_acf_all_mid  <- -0.10
    mean_min_acf_all_low  <- -0.15
    mean_mC_all_high  <- 0.02500000
    mean_mC_all_mid   <- 0.01737801
    mean_mC_all_low   <- 0.01000000
  }
}


ggTrend_mean_mC_all_fk_kappa_all <- trendPlot(dataFrame = con_fk_df_all,
                                              mapping = aes(x = mean_mC_all, y = fk_kappa_all),
                                              xvar = mean_mC_all,
                                              yvar = fk_kappa_all,
                                              xlab = bquote(.(featName)*" mean m"*.(context)),
                                              ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                              xaxtrans = log10_trans(),
                                              yaxtrans = log10_trans(),
                                              xbreaks = trans_breaks("log10", function(x) 10^x),
                                              ybreaks = trans_breaks("log10", function(x) 10^x),
                                              xlabels = trans_format("log10", math_format(10^.x)),
                                              ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_fk_kappa_all <- ggTrend_mean_mC_all_fk_kappa_all +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                   mapping = aes(x = mean_mC_all, y = fk_kappa_all),
                                                   xvar = mean_mC_all,
                                                   yvar = fk_kappa_all,
                                                   xlab = bquote(.(featName)*" mean m"*.(context)),
                                                   ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                                   xaxtrans = log10_trans(),
                                                   yaxtrans = log10_trans(),
                                                   xbreaks = trans_breaks("log10", function(x) 10^x),
                                                   ybreaks = trans_breaks("log10", function(x) 10^x),
                                                   xlabels = trans_format("log10", math_format(10^.x)),
                                                   ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_fk_kappa_all_filt <- ggTrend_mean_mC_all_fk_kappa_all_filt +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_mean_mC_all_mean_stocha_all <- trendPlot(dataFrame = con_fk_df_all,
                                                 mapping = aes(x = mean_mC_all, y = mean_stocha_all),
                                                 xvar = mean_mC_all,
                                                 yvar = mean_stocha_all,
                                                 xlab = bquote(.(featName)*" mean m"*.(context)),
                                                 ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                 xaxtrans = log10_trans(),
                                                 yaxtrans = log10_trans(),
                                                 xbreaks = trans_breaks("log10", function(x) 10^x),
                                                 ybreaks = trans_breaks("log10", function(x) 10^x),
                                                 xlabels = trans_format("log10", math_format(10^.x)),
                                                 ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_mean_stocha_all <- ggTrend_mean_mC_all_mean_stocha_all +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_mean_stocha_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                      mapping = aes(x = mean_mC_all, y = mean_stocha_all),
                                                      xvar = mean_mC_all,
                                                      yvar = mean_stocha_all,
                                                      xlab = bquote(.(featName)*" mean m"*.(context)),
                                                      ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                      xaxtrans = log10_trans(),
                                                      yaxtrans = log10_trans(),
                                                      xbreaks = trans_breaks("log10", function(x) 10^x),
                                                      ybreaks = trans_breaks("log10", function(x) 10^x),
                                                      xlabels = trans_format("log10", math_format(10^.x)),
                                                      ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_mean_stocha_all_filt <- ggTrend_mean_mC_all_mean_stocha_all_filt +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_fk_reads_all_fk_kappa_all <- trendPlot(dataFrame = con_fk_df_all,
                                               mapping = aes(x = fk_reads_all, y = fk_kappa_all),
                                               xvar = fk_reads_all,
                                               yvar = fk_kappa_all,
                                               xlab = bquote(.(featName)*" read coverage"),
                                               ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                               xaxtrans = log10_trans(),
                                               yaxtrans = log10_trans(),
                                               xbreaks = trans_breaks("log10", function(x) 10^x),
                                               ybreaks = trans_breaks("log10", function(x) 10^x),
                                               xlabels = trans_format("log10", math_format(10^.x)),
                                               ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_reads_all_fk_kappa_all <- ggTrend_fk_reads_all_fk_kappa_all +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_reads_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                    mapping = aes(x = fk_reads_all, y = fk_kappa_all),
                                                    xvar = fk_reads_all,
                                                    yvar = fk_kappa_all,
                                                    xlab = bquote(.(featName)*" read coverage"),
                                                    ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                                    xaxtrans = log10_trans(),
                                                    yaxtrans = log10_trans(),
                                                    xbreaks = trans_breaks("log10", function(x) 10^x),
                                                    ybreaks = trans_breaks("log10", function(x) 10^x),
                                                    xlabels = trans_format("log10", math_format(10^.x)),
                                                    ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_reads_all_fk_kappa_all_filt <- ggTrend_fk_reads_all_fk_kappa_all_filt +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_fk_reads_all_mean_stocha_all <- trendPlot(dataFrame = con_fk_df_all,
                                                  mapping = aes(x = fk_reads_all, y = mean_stocha_all),
                                                  xvar = fk_reads_all,
                                                  yvar = mean_stocha_all,
                                                  xlab = bquote(.(featName)*" read coverage"),
                                                  ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                  xaxtrans = log10_trans(),
                                                  yaxtrans = log10_trans(),
                                                  xbreaks = trans_breaks("log10", function(x) 10^x),
                                                  ybreaks = trans_breaks("log10", function(x) 10^x),
                                                  xlabels = trans_format("log10", math_format(10^.x)),
                                                  ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_reads_all_mean_stocha_all <- ggTrend_fk_reads_all_mean_stocha_all +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_reads_all_mean_stocha_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                       mapping = aes(x = fk_reads_all, y = mean_stocha_all),
                                                       xvar = fk_reads_all,
                                                       yvar = mean_stocha_all,
                                                       xlab = bquote(.(featName)*" read coverage"),
                                                       ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_reads_all_mean_stocha_all_filt <- ggTrend_fk_reads_all_mean_stocha_all_filt +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_fk_Cs_all_fk_kappa_all <- trendPlot(dataFrame = con_fk_df_all,
                                            mapping = aes(x = fk_Cs_all, y = fk_kappa_all),
                                            xvar = fk_Cs_all,
                                            yvar = fk_kappa_all,
                                            xlab = bquote(.(featName)*" # cytosines (m"*.(context)*")"),
                                            ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                            xaxtrans = log10_trans(),
                                            yaxtrans = log10_trans(),
                                            xbreaks = trans_breaks("log10", function(x) 10^x),
                                            ybreaks = trans_breaks("log10", function(x) 10^x),
                                            xlabels = trans_format("log10", math_format(10^.x)),
                                            ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_Cs_all_fk_kappa_all <- ggTrend_fk_Cs_all_fk_kappa_all +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_Cs_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                 mapping = aes(x = fk_Cs_all, y = fk_kappa_all),
                                                 xvar = fk_Cs_all,
                                                 yvar = fk_kappa_all,
                                                 xlab = bquote(.(featName)*" # cytosines (m"*.(context)*")"),
                                                 ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                                 xaxtrans = log10_trans(),
                                                 yaxtrans = log10_trans(),
                                                 xbreaks = trans_breaks("log10", function(x) 10^x),
                                                 ybreaks = trans_breaks("log10", function(x) 10^x),
                                                 xlabels = trans_format("log10", math_format(10^.x)),
                                                 ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_Cs_all_fk_kappa_all_filt <- ggTrend_fk_Cs_all_fk_kappa_all_filt +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_fk_Cs_all_mean_stocha_all <- trendPlot(dataFrame = con_fk_df_all,
                                               mapping = aes(x = fk_Cs_all, y = mean_stocha_all),
                                               xvar = fk_Cs_all,
                                               yvar = mean_stocha_all,
                                               xlab = bquote(.(featName)*" # cytosines (m"*.(context)*")"),
                                               ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                               xaxtrans = log10_trans(),
                                               yaxtrans = log10_trans(),
                                               xbreaks = trans_breaks("log10", function(x) 10^x),
                                               ybreaks = trans_breaks("log10", function(x) 10^x),
                                               xlabels = trans_format("log10", math_format(10^.x)),
                                               ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_Cs_all_mean_stocha_all <- ggTrend_fk_Cs_all_mean_stocha_all +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_Cs_all_mean_stocha_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                    mapping = aes(x = fk_Cs_all, y = mean_stocha_all),
                                                    xvar = fk_Cs_all,
                                                    yvar = mean_stocha_all,
                                                    xlab = bquote(.(featName)*" # cytosines (m"*.(context)*")"),
                                                    ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                    xaxtrans = log10_trans(),
                                                    yaxtrans = log10_trans(),
                                                    xbreaks = trans_breaks("log10", function(x) 10^x),
                                                    ybreaks = trans_breaks("log10", function(x) 10^x),
                                                    xlabels = trans_format("log10", math_format(10^.x)),
                                                    ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_Cs_all_mean_stocha_all_filt <- ggTrend_fk_Cs_all_mean_stocha_all_filt +
  geom_hline(yintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


ggTrend_mean_stocha_all_fk_kappa_all <- trendPlot(dataFrame = con_fk_df_all,
                                                  mapping = aes(x = mean_stocha_all, y = fk_kappa_all),
                                                  xvar = mean_stocha_all,
                                                  yvar = fk_kappa_all,
                                                  xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                  ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                                  xaxtrans = log10_trans(),
                                                  yaxtrans = log10_trans(),
                                                  xbreaks = trans_breaks("log10", function(x) 10^x),
                                                  ybreaks = trans_breaks("log10", function(x) 10^x),
                                                  xlabels = trans_format("log10", math_format(10^.x)),
                                                  ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_stocha_all_fk_kappa_all <- ggTrend_mean_stocha_all_fk_kappa_all +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all),
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_stocha_all_fk_kappa_all_filt <- ggTrend_mean_stocha_all_fk_kappa_all_filt +
  geom_hline(yintercept = fk_kappa_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_hline(yintercept = fk_kappa_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_hline(yintercept = fk_kappa_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  geom_vline(xintercept = mean_stocha_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
  geom_vline(xintercept = mean_stocha_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
  geom_vline(xintercept = mean_stocha_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
  facet_grid(cols = vars(chr), scales = "free_x")


gg_cow_list1 <- list(
                     ggTrend_mean_mC_all_fk_kappa_all,
                     ggTrend_mean_mC_all_fk_kappa_all_filt,
                     ggTrend_mean_mC_all_mean_stocha_all,
                     ggTrend_mean_mC_all_mean_stocha_all_filt,
                     ggTrend_fk_reads_all_fk_kappa_all,
                     ggTrend_fk_reads_all_fk_kappa_all_filt,
                     ggTrend_fk_reads_all_mean_stocha_all,
                     ggTrend_fk_reads_all_mean_stocha_all_filt,
                     ggTrend_fk_Cs_all_fk_kappa_all,
                     ggTrend_fk_Cs_all_fk_kappa_all_filt,
                     ggTrend_fk_Cs_all_mean_stocha_all,
                     ggTrend_fk_Cs_all_mean_stocha_all_filt,
                     ggTrend_mean_stocha_all_fk_kappa_all,
                     ggTrend_mean_stocha_all_fk_kappa_all_filt
                    )

gg_cow1 <- plot_grid(plotlist = gg_cow_list1,
                     labels = "AUTO", label_size = 30,
                     align = "hv",
                     axis = "l",
                     nrow = length(gg_cow_list1), ncol = 1)

ggsave(paste0(plotDir,
              featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase, "_", context,
              "_NAmax", NAmax, "_all_trendPlot_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow1,
       height = 5*length(gg_cow_list1), width = 5*length(chrName), limitsize = F)


# Extract feature groups (based on trend plots) to enable enrichment analysis

# Filter by fk_kappa_all and mean_mC_all
if(context == "CpG") {
  con_fk_df_all_filt_kappa_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all > fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  > mean_mC_all_high)
} else if(context == "CHG") {
  con_fk_df_all_filt_kappa_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all > fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  > mean_mC_all_high)
} else if(context == "CHH") {
  con_fk_df_all_filt_kappa_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_low)
  
  con_fk_df_all_filt_kappa_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_mid)
  
  con_fk_df_all_filt_kappa_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all >  fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all  <= mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all <= fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  >  mean_mC_all_high)
  
  con_fk_df_all_filt_kappa_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all > fk_kappa_all_high) %>%
    dplyr::filter(mean_mC_all  > mean_mC_all_high)
}

write.table(con_fk_df_all_filt_kappa_mC_group1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group5,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group6,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group7,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_group8,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all
if(context == "CpG") {
  con_fk_df_all_filt_stocha_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)
  
  con_fk_df_all_filt_stocha_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)
  
  con_fk_df_all_filt_stocha_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)
  
  con_fk_df_all_filt_stocha_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)
  
  con_fk_df_all_filt_stocha_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     > mean_mC_all_high)
} else if(context == "CHG") {
  con_fk_df_all_filt_stocha_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)
  
  con_fk_df_all_filt_stocha_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)
  
  con_fk_df_all_filt_stocha_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)
  
  con_fk_df_all_filt_stocha_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)
  
  con_fk_df_all_filt_stocha_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_high)
  
  con_fk_df_all_filt_stocha_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     > mean_mC_all_high)
} else if(context == "CHH") {
  con_fk_df_all_filt_stocha_mC_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)

  con_fk_df_all_filt_stocha_mC_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_low)

  con_fk_df_all_filt_stocha_mC_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

  con_fk_df_all_filt_stocha_mC_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

  con_fk_df_all_filt_stocha_mC_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)

  con_fk_df_all_filt_stocha_mC_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
    dplyr::filter(mean_mC_all     <= mean_mC_all_high)

  con_fk_df_all_filt_stocha_mC_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     >  mean_mC_all_high)

  con_fk_df_all_filt_stocha_mC_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_high) %>%
    dplyr::filter(mean_mC_all     > mean_mC_all_high)
}


write.table(con_fk_df_all_filt_stocha_mC_group1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group5,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group6,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group7,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_group8,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by fk_kappa_all and mean_stocha_all
if(context == "CpG") {
  con_fk_df_all_filt_kappa_stocha_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    > fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_high)
} else if(context == "CHG") {
  con_fk_df_all_filt_kappa_stocha_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    > fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_high)
} else if(context == "CHH") {
  con_fk_df_all_filt_kappa_stocha_group1 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group2 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_low)
  
  con_fk_df_all_filt_kappa_stocha_group3 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group4 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_mid) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_mid)
  
  con_fk_df_all_filt_kappa_stocha_group5 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group6 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    >  fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
    dplyr::filter(mean_stocha_all <= mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group7 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    <= fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all >  mean_stocha_all_high)
  
  con_fk_df_all_filt_kappa_stocha_group8 <- con_fk_df_all_filt %>%
    dplyr::filter(fk_kappa_all    > fk_kappa_all_high) %>%
    dplyr::filter(mean_stocha_all > mean_stocha_all_high)
}

write.table(con_fk_df_all_filt_kappa_stocha_group1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group5,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group6,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group7,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_stocha_group8,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_stocha_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
