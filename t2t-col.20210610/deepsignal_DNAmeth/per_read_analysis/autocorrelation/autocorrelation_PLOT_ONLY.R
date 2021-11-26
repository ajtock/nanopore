#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Genome-wide autocorrelation of context-specific cytosine methylation status
# at increasing physical distances (e.g., 1 to 10,000 nucleotides)

# Usage on hydrogen node7:
# csmit -m 500G -c 47 "/applications/R/R-4.0.0/bin/Rscript autocorrelation_PLOT_ONLY.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 1e4 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 500 CEN180"

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#nperm <- 1e4
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#maxDist <- 500
#featName <- "CEN180"
#min_pval <- 1 - ( (nperm - 1) / nperm )

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
nperm <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
maxDist <- as.numeric(args[7])
featName <- args[8]
min_pval <- 1 - ( (nperm - 1) / nperm )

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(tidyquant)
library(fastmatch)

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

tabGR_feat_all_acf_df <- read.table(paste0(sampleName, "_MappedOn_", refbase, "_", context,
                                           "_all_autocorrelation_", featName, "_", paste0(chrName, collapse = "_"),
                                           ".tsv"), header = T)
tabGR_feat_all_acf_df <- data.frame(tabGR_feat_all_acf_df,
                                    adj_pval = -log10(p.adjust(10^-tabGR_feat_all_acf_df$pval, method = "BH")))

# Chromosome-scale plotting function
chrPlot <- function(dataFrame, xvar, yvar, xlab, ylab, colour) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) +
#  geom_line(colour = colour, size = 1) +
  geom_ma(ma_fun = SMA, n = 6, colour = colour, linetype = 1, size = 2) +
  scale_x_continuous(
                     labels = function(x) x) +
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
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}


gg_tabGR_feat_all_acf <- chrPlot(dataFrame = tabGR_feat_all_acf_df,
                                 xvar = distance,
                                 yvar = acf,
                                 xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                 ylab = bquote("Observed correlation (m"*.(context)*")"),
                                 colour = "dodgerblue")
gg_tabGR_feat_all_acf <- gg_tabGR_feat_all_acf +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_feat_all_adj_pval <- chrPlot(dataFrame = tabGR_feat_all_acf_df,
                                      xvar = distance,
                                      yvar = adj_pval,
                                      xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                      ylab = bquote("-"*Log[10]*"(BH-adj. "*italic(P)*"-value) (m"*.(context)*")"),
                                      colour = "red") +
  geom_hline(yintercept = -log10(0.05), colour = "black", size = 1, linetype = "dashed")
gg_tabGR_feat_all_adj_pval <- gg_tabGR_feat_all_adj_pval +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_feat_all_exp <- chrPlot(dataFrame = tabGR_feat_all_acf_df,
                                 xvar = distance,
                                 yvar = exp,
                                 xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                 ylab = bquote("Mean permuted correlation (m"*.(context)*")"),
                                 colour = "lightseagreen")
gg_tabGR_feat_all_exp <- gg_tabGR_feat_all_exp +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_all_list <- list(
                        gg_tabGR_feat_all_acf,
                        gg_tabGR_feat_all_adj_pval,
                        gg_tabGR_feat_all_exp
                       )
gg_cow_all <- plot_grid(plotlist = gg_cow_all_list,
                        labels = c("AUTO"), label_size = 30,
                        align = "hv",
                        axis = "l",
                        nrow = length(gg_cow_all_list), ncol = 1)

ggsave(paste0(plotDir,
              sampleName, "_MappedOn_", refbase, "_", context,
              "_all_autocorrelation_", featName, "_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow_all,
       height = 5*length(gg_cow_all_list), width = 10*length(chrName), limitsize = F)




#tabGR_feat_fwd_chr_dist_bool_list <- mclapply(seq_along(tabGR_feat_fwd_chr_dist_list), function(y) {
#  sapply(1:maxDist, function(z) {
#    z %fin% tabGR_feat_fwd_chr_dist_list[[y]]
#  })
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
#
#tabGR_feat_fwd_chr_dist_bool_list <- lapply(seq_alonglapply(1:maxDist, function(z) {
#  unlist(mclapply(seq_along(tabGR_feat_fwd_chr_dist_list), function(y) {
#    z %fin% tabGR_feat_fwd_chr_dist_list[[y]]
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
#})
#
#tabGR_feat_fwd_chr_dist_bool_list <- mclapply(1:maxDist, function(z) {
#  sapply(seq_along(tabGR_feat_fwd_chr_dist_list), function(y) {
#    z %fin% tabGR_feat_fwd_chr_dist_list[[y]]
#  })
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)
