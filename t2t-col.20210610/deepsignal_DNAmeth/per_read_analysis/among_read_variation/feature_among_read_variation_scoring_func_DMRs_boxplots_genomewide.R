#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Plot among-read variation/agreement (e.g., Fleiss' kappa) and stochasticity for each TE superfamily (e.g., as boxplots or violin plots)

# Usage:
# /applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_DMRs_boxplots_genomewide.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'cmt2_BSseq_Rep1_hypoCHG,cmt3_BSseq_Rep1_hypoCHG,kss_BSseq_Rep1_hypoCHG' 'bodies'
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CHG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- unlist(strsplit("cmt2_BSseq_Rep1_hypoCHG,cmt3_BSseq_Rep1_hypoCHG,kss_BSseq_Rep1_hypoCHG", split = ","))
#featRegion <- "bodies"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- unlist(strsplit(args[7], split = ","))
featRegion <- args[8]

featNamePlot <- sub("_BSseq_Rep1_", " ", featName)
featNamePlot <- sub("_", " ", featNamePlot)

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggthemes)
library(viridis)
library(pals)

inDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
outDir <- paste0("hypoCHG_DMRs_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

score_colFun <- cols25(n = 25)[-c(7:16, 25)][1:length(featNamePlot)]
stopifnot(length(score_colFun) == length(featNamePlot))
#names(score_colFun) <- featNamePlot

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Load among-read and within-read mC data for featName featRegion
con_fk_df_all_list <- lapply(1:length(featName), function(x) {
  tmp <- read.table(paste0(inDir[x],
                           featName[x], "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T)
  tmp$Feature <- featNamePlot[x]  
  tmp
})
con_fk_df_all <- dplyr::bind_rows(con_fk_df_all_list)
con_fk_df_all$Feature <- factor(con_fk_df_all$Feature,
                                levels = c(featNamePlot))

con_fk_df_all_filt_list <- lapply(1:length(featName), function(x) {
  tmp <- read.table(paste0(inDir[x],
                           featName[x], "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T)
  tmp$Feature <- featNamePlot[x]  
  tmp
})
con_fk_df_all_filt <- dplyr::bind_rows(con_fk_df_all_filt_list)
con_fk_df_all_filt$Feature <- factor(con_fk_df_all_filt$Feature,
                                     levels = c(featNamePlot))


# Plot relationships and define groups
boxPlot <- function(dataFrame, mapping, xvar, yvar, xlab, ylab, yaxtrans, ybreaks, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = mapping,
         colour = Feature) +
  scale_colour_manual(values = score_colFun) +
  geom_boxplot(notch = T,
               varwidth = T) +
#  geom_violin(scale = "area",
#              trim = T,
#              draw_quantiles = c(0.25, 0.50, 0.75)) +
#  geom_beeswarm(cex = 6,
#                size = 4) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  labs(x = xlab,
       y = ylab) +
#  geom_hline(yintercept = mean(dplyr::mutate(dataFrame, yvar = !!yvar)$yvar, na.rm = T), linetype = "dashed", size = 1, colour = "darkorange1") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black", angle = 45, vjust = 1.0, hjust = 1.0),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        panel.grid = element_blank(),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18))
}

gg_fk_kappa_all <- boxPlot(dataFrame = con_fk_df_all,
                           mapping = aes(x = Feature, y = fk_kappa_all, colour = Feature),
                           xvar = score,
                           yvar = fk_kappa_all,
                           xlab = "DMRs",
                           ylab  = bquote("DMR agreement (m"*.(context)*")"),
                           yaxtrans = "identity",
                           ybreaks = waiver(),
                           ylabels = waiver())

gg_fk_kappa_all_filt <- boxPlot(dataFrame = con_fk_df_all,
                                mapping = aes(x = Feature, y = fk_kappa_all, colour = Feature),
                                xvar = score,
                                yvar = fk_kappa_all,
                                xlab = "DMRs",
                                ylab  = bquote("DMR agreement (m"*.(context)*")"),
                                yaxtrans = "identity",
                                ybreaks = waiver(),
                                ylabels = waiver())

                                 ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                             ylab = bquote(.(featName)*" mean m"*.(context)),

gg_cow_list1 <- list(
                     gg_fk_kappa_all,
                     gg_fk_kappa_all_filt
#                     gg_mean_stocha_all,
                     gg_mean_stocha_all_filt,
#                     gg_mean_mC_all,
                     gg_mean_mC_all_filt
                    )

gg_cow1 <- plot_grid(plotlist = gg_cow_list1,
#                     labels = "AUTO", label_size = 30,
                     align = "hv",
                     axis = "l",
                     nrow = length(gg_cow_list1), ncol = 1)

ggsave(paste0(plotDir,
              "hypoCHG_DMRs_", featRegion, "_", sampleName, "_MappedOn_", refbase, "_", context,
              "_NAmax", NAmax, "_all_boxPlot_", paste0(chrName, collapse = "_"),
              "_genomewide.pdf"),
       plot = gg_cow1,
       height = 6*length(gg_cow_list1), width = 10, limitsize = F)
