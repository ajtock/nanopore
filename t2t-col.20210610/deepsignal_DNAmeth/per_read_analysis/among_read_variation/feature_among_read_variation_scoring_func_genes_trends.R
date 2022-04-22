#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Score among-read variation/agreement (e.g., Fleiss' kappa) for each feature
# 2. Examine relationships between feature among-read agreement and other metrics

# Usage:
# /applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_genes_trends.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'bodies'
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
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
library(stringr)
library(data.table)
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
con_fk_df_all$parent <- sub(pattern = "\\.\\d+", replacement = "", x = con_fk_df_all$name) 
con_fk_df_all$parent <- sub(pattern = "_\\d+", replacement = "", x = con_fk_df_all$parent) 

con_fk_df_all_filt <- read.table(paste0(outDir,
                                        featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                                        "_", context,
                                        "_NAmax", NAmax,
                                        "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                                        paste0(chrName, collapse = "_"), ".tsv"),
                                 header = T)
con_fk_df_all_filt$parent <- sub(pattern = "\\.\\d+", replacement = "", x = con_fk_df_all_filt$name) 
con_fk_df_all_filt$parent <- sub(pattern = "_\\d+", replacement = "", x = con_fk_df_all_filt$parent) 

# Append intron retention ratio (calculated with IRFinder)
Col_Rep3_IRFinder <- fread(paste0("/home/ajt200/analysis/RNAseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_RNAseq_IRFinder_TAIR10_chr_all/REF/TAIR10_chr_all/",
                                  "Col_0_RNAseq_Rep3_ERR966159/IRFinder-IR-nondir.txt"),
                           sep = "\t", data.table = F)
Col_Rep3_IRFinder <- Col_Rep3_IRFinder[grep("clean", Col_Rep3_IRFinder$Name),]
nrow(Col_Rep3_IRFinder[which(Col_Rep3_IRFinder$Warnings == "-"),])
#[1] 45907
#[1] 22136
Col_Rep3_IRFinder$Name <- str_extract(Col_Rep3_IRFinder$Name, "AT\\wG\\d+")
Col_Rep3_IRFinder <- Col_Rep3_IRFinder[-which(is.na(Col_Rep3_IRFinder$Name)),]

library(doFuture)
registerDoFuture()
plan(multicore)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

parentIDs <- unique(Col_Rep3_IRFinder$Name)

Col_Rep3_IRratio <- foreach(i = iter(parentIDs),
                            .combine = "rbind",
                            .multicombine = T,
                            .maxcombine = length(parentIDs)+1e1,
                            .inorder = F,
                            .errorhandling = "pass") %dopar% {
  tmpDF <- Col_Rep3_IRFinder[which(Col_Rep3_IRFinder$Name == i),]
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

con_fk_df_all_tab <- base::merge(x = con_fk_df_all, y = Col_Rep3_IRratio,
                                 by.x = "parent", by.y = "parent")
con_fk_df_all_filt_tab <- base::merge(x = con_fk_df_all_filt, y = Col_Rep3_IRratio,
                                      by.x = "parent", by.y = "parent")

print(cor.test(con_fk_df_all_tab$fk_kappa_all, con_fk_df_all_tab$IRratio_mean, method = "spearman"))
print(cor.test(con_fk_df_all_tab$fk_kappa_all, con_fk_df_all_tab$IRratio_median, method = "spearman"))
print(cor.test(con_fk_df_all_tab$mean_stocha_all, con_fk_df_all_tab$IRratio_mean, method = "spearman"))
print(cor.test(con_fk_df_all_tab$mean_stocha_all, con_fk_df_all_tab$IRratio_median, method = "spearman"))
print(cor.test(con_fk_df_all_filt_tab$fk_kappa_all, con_fk_df_all_filt_tab$IRratio_mean, method = "spearman"))
print(cor.test(con_fk_df_all_filt_tab$fk_kappa_all, con_fk_df_all_filt_tab$IRratio_median, method = "spearman"))
print(cor.test(con_fk_df_all_filt_tab$mean_stocha_all, con_fk_df_all_filt_tab$IRratio_mean, method = "spearman"))
print(cor.test(con_fk_df_all_filt_tab$mean_stocha_all, con_fk_df_all_filt_tab$IRratio_median, method = "spearman"))


# Plot relationships and define groups
trendPlot <- function(dataFrame, mapping, xvar, yvar, xlab, ylab, xaxtrans, yaxtrans, xbreaks, ybreaks, xlabels, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = mapping) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = xaxtrans,
                     breaks = xbreaks,
                     labels = xlabels) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y ~ s(x, bs = "cs")) +
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
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame, !!enquo(xvar))[,1], select(dataFrame, !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame, !!enquo(xvar))[,1], select(dataFrame, !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame)[1]/100) )),
                         digits = 5))))
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
