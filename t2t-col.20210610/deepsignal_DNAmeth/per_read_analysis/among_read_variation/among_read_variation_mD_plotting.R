#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Among-read variation/agreement (Fleiss' kappa) [this script]
# 2. Variance of read stochasticity (variance of per-read, per-window proportion of inter-C intervals that represent a methylation state change [variance across all reads that overlap window examined])
# 3. Per-read autocorrelation, but could be hard to derive a general measure across reads overlapping a given region (maybe look at variance of pairwise read correlation coefficients between autocorrelation values, either standardised to e.g. 100 autocorrelation values per read, or use only reads with info across the same Cs) 

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# /applications/R/R-4.0.0/bin/Rscript among_read_variation_mD_plotting.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5'
# /applications/R/R-4.0.0/bin/Rscript among_read_variation_mD_plotting.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5'

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
genomeBinSize <- as.integer(args[3])
genomeStepSize <- as.integer(args[4])
context <- args[5]
NAmax <- as.numeric(args[6])
CPUpc <- as.numeric(args[7])
chrName <- unlist(strsplit(args[8], split = ","))

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
#library(GenomicRanges)
#library(irr)
library(dplyr)
library(tidyr)
#library(fpc)
##library(data.table)
##library(segmentSeq)
#library(ComplexHeatmap)
#library(RColorBrewer)
#library(circlize)
 
library(ggplot2)
library(cowplot)
library(scales)
#library(ggcorrplot)
library(viridis)
library(ggthemes)
library(tidyquant)
#library(grid)

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

outDir <- paste0("genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "/")
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

# Read in among-read (Fliess' kappa) and within-read (stochasticity) methylation variation files
tab_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0(outDir,
                    sampleName, "_MappedOn_", refbase, "_", context,
                    "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                    "_NAmax", NAmax, "_per_read_var_df_", chrName[x], ".tsv"),
             header = T) 
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  tab <- dplyr::bind_rows(tab_list)
} else {
  tab <- tab_list[[1]]
}
rm(tab_list); gc()

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

tab_mD <- merge(x = tab, y = MA1_2_mD,
                by.x = c("chr", "start", "end"),
                by.y = c("chr", "start", "end"))
tab_mD <- tab_mD[
  with(tab_mD, order(chr, start, end)),
]
stopifnot(identical(tab_mD$fk_kappa_all, tab$fk_kappa_all))
stopifnot(identical(tab_mD$MA1_2_mean.D, MA1_2_mD$MA1_2_mean.D))

# Remove genomic bins tab_mD overlapping the mitochondrial insertion on Chr2,
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


# Load features chromosome-scale profiles
genes <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/genes/", refbase, "_gene_frequency_per_",
                           genomeBinName, "_unsmoothed.tsv"),
                    header = T)
genes$midpoint <- (genes$window-1) + (genomeBinSize/2)
genes <- genes[genes$chr %in% chrName,]

gypsy <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR_frequency_per_",
                           genomeBinName, "_unsmoothed.tsv"),
                    header = T)
gypsy$midpoint <- (gypsy$window-1) + (genomeBinSize/2)
gypsy <- gypsy[gypsy$chr %in% chrName,]

CEN180 <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/CEN180/", refbase, "_CEN180_frequency_per_",
                           genomeBinName, "_unsmoothed.tsv"),
                    header = T)
CEN180$midpoint <- (CEN180$window-1) + (genomeBinSize/2)
CEN180 <- CEN180[CEN180$chr %in% chrName,]


chrPlot <- function(dataFrame, xvar, yvar, xlab, ylab, colour) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) +
#  geom_line(size = 1, colour = colour) +
  geom_ma(ma_fun = SMA, n = 10, colour = colour, linetype = 1, size = 1) +
  scale_x_continuous(
                     labels = function(x) x/1e6) +
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
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}

chrPlot2 <- function(dataFrame, xvar, yvar1, yvar2, yvar1lab, yvar2lab, ylims, xlab, ylab, colour1, colour2) {
  xvar <- enquo(xvar)
  yvar1 <- enquo(yvar1)
  yvar2 <- enquo(yvar2)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar)) +
  geom_ma(ma_fun = SMA, mapping = aes(y = !!yvar1, colour = colour1), n = 10, linetype = 1, size = 1) +
  geom_ma(ma_fun = SMA, mapping = aes(y = !!yvar2, colour = colour2), n = 10, linetype = 1, size = 1) +
  scale_colour_identity(name = NULL,
                        breaks = c(colour1, colour2),
                        labels = c(yvar1lab, yvar2lab),
                        guide = "legend") +
  coord_cartesian(ylim = ylims) +
  scale_x_continuous(labels = function(x) x/1e6) +
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
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        legend.position = c(0.01, 0.975),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.box.background = element_rect(colour = "black", size = 1),
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}

chrPlot3 <- function(dataFrame, xvar, yvar, yvarmin, yvarmax, xlab, ylab, colour) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  yvarmin <- enquo(yvarmin)
  yvarmax <- enquo(yvarmax)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) +
#  geom_line(size = 1, colour = colour) +
  geom_bbands(data = dataFrame,
              mapping = aes(low = !!yvarmin,
                            high = !!yvarmax,
                            close = !!yvar),
              ma_fun = SMA, n = 10,
              color_ma = colour,
              color_bands = NA,
              fill = colour,
              alpha = 0.4,
              linetype = 1, size = 1) +
#  geom_ma(ma_fun = SMA, n = 10, colour = colour, linetype = 1, size = 1) +
#  geom_ribbon(data = dataFrame,
#              mapping = aes(ymin = !!xvarmin,
#                            ymax = !!xvarmax,
#                            fill = colour),
#              alpha = 0.4) + 
  scale_x_continuous(
                     labels = function(x) x/1e6) +
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
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}


mito_ins_bounds <- data.frame(chr = chrs,
                              xmin = c(NA, mito_ins_start, rep(NA, 3)),
                              xmax = c(NA, mito_ins_end, rep(NA, 3)),
                              ymin = c(NA, -Inf, rep(NA, 3)),
                              ymax = c(NA, Inf, rep(NA, 3)))

tab2 <- tab %>%
  mutate(fk_total_reads_all = tab$fk_num_reads_all/tab$fk_prop_reads_all,
         fk_total_Cs_all = tab$fk_num_Cs_all/tab$fk_prop_Cs_all)

gg_fk_num_reads_all <- chrPlot2(dataFrame = tab2,
                                xvar = midpoint,
                                yvar1 = fk_num_reads_all,
                                yvar2 = fk_total_reads_all,
                                yvar1lab = bquote("Kept (NAs at <" * .(NAmax*100) * "% of sites)"),
                                yvar2lab = "Total",
                                ylims = c(0, 100),
                                xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                                ylab = bquote("No. of reads (m"*.(context)*")"),
                                colour1 = "blue",
                                colour2 = "grey50") 
gg_fk_num_reads_all <- gg_fk_num_reads_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_fk_prop_reads_all <- chrPlot(dataFrame = tab,
                                xvar = midpoint,
                                yvar = fk_prop_reads_all,
                                xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                                ylab = bquote("Proportion of reads kept (m"*.(context)*")"),
                                colour = "blue") 
gg_fk_prop_reads_all <- gg_fk_prop_reads_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")


gg_fk_num_Cs_all <- chrPlot2(dataFrame = tab2,
                             xvar = midpoint,
                             yvar1 = fk_num_Cs_all,
                             yvar2 = fk_total_Cs_all,
                             yvar1lab = "Kept for Fleiss' kappa",
                             yvar2lab = "Total",
                             ylims = NULL,
                             xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                             ylab = bquote("No. of sites (m"*.(context)*")"),
                             colour1 = "red",
                             colour2 = "grey50") 
gg_fk_num_Cs_all <- gg_fk_num_Cs_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_fk_prop_Cs_all <- chrPlot(dataFrame = tab,
                             xvar = midpoint,
                             yvar = fk_prop_Cs_all,
                             xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                             ylab = bquote("Proportion of sites kept (m"*.(context)*")"),
                             colour = "red") 
gg_fk_prop_Cs_all <- gg_fk_prop_Cs_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_MA1_2_D <- chrPlot3(dataFrame = tab_mD,
                       xvar = midpoint,
                       yvar = MA1_2_mean.D,
                       yvarmin = MA1_2_min.D,
                       yvarmax = MA1_2_max.D,
                       xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                       ylab = bquote("MA1_2" ~ italic("D") ~ "(m"*.(context)*")"),
                       colour = "darkgreen") 
gg_MA1_2_D <- gg_MA1_2_D +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_MA1_2_mean.D <- chrPlot(dataFrame = tab_mD,
                           xvar = midpoint,
                           yvar = MA1_2_mean.D,
                           xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                           ylab = bquote("MA1_2" ~ italic("D") ~ italic(bar(x)) ~ "(m"*.(context)*")"),
                           colour = "darkgreen")
gg_MA1_2_mean.D <- gg_MA1_2_mean.D +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_fk_kappa_all <- chrPlot(dataFrame = tab,
                           xvar = midpoint,
                           yvar = fk_kappa_all,
                           xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                           ylab = bquote("Fleiss' kappa (m"*.(context)*")"),
                           colour = "red") 
gg_fk_kappa_all <- gg_fk_kappa_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")


gg_mean_stocha_all <- chrPlot(dataFrame = tab,
                              xvar = midpoint,
                              yvar = mean_stocha_all,
                              xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                              ylab = bquote("Mean stochasticity (m"*.(context)*")"),
                              colour = "dodgerblue") 
gg_mean_stocha_all <- gg_mean_stocha_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_mean_mean_acf_all <- chrPlot(dataFrame = tab,
                                xvar = midpoint,
                                yvar = mean_mean_acf_all,
                                xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                                ylab = bquote("Mean mean ACF (m"*.(context)*")"),
                                colour = "darkgreen") 
gg_mean_mean_acf_all <- gg_mean_mean_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_mean_min_acf_all <- chrPlot(dataFrame = tab,
                               xvar = midpoint,
                               yvar = mean_min_acf_all,
                               xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                               ylab = bquote("Mean min. ACF (m"*.(context)*")"),
                               colour = "forestgreen") 
gg_mean_min_acf_all <- gg_mean_min_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_genes <- chrPlot(dataFrame = genes,
                    xvar = midpoint,
                    yvar = features,
                    xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                    ylab = bquote("Genes"),
                    colour = "lightseagreen") 
gg_genes <- gg_genes +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_gypsy <- chrPlot(dataFrame = gypsy,
                    xvar = midpoint,
                    yvar = features,
                    xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                    ylab = bquote(italic(GYPSY)),
                    colour = "purple4") 
gg_gypsy <- gg_gypsy +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")

gg_CEN180 <- chrPlot(dataFrame = CEN180,
                     xvar = midpoint,
                     yvar = features,
                     xlab = paste0("Coordinates (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
                     ylab = bquote(italic("CEN180")),
                     colour = "darkorange") 
gg_CEN180 <- gg_CEN180 +
  facet_grid(cols = vars(chr), scales = "free_x") +
  geom_rect(data = mito_ins_bounds,
            mapping = aes(x = NULL, y = NULL,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
            alpha = 0.4, fill = "grey50")


gg_cow_list <- list(
                    gg_MA1_2_D, gg_MA1_2_mean.D, gg_fk_kappa_all, gg_mean_stocha_all, gg_mean_mean_acf_all, gg_mean_min_acf_all,
                    gg_genes, gg_gypsy, gg_CEN180,
                    gg_fk_num_reads_all, gg_fk_num_Cs_all
                   )
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = toupper(letters[1:length(gg_cow_list)]), label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDir,
              sampleName, "_MappedOn_", refbase, "_", context,
              "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_NAmax", NAmax, "_all_chrPlot_", paste0(chrName, collapse = "_"),
              "_incl_MA1_2_mD.pdf"),
       plot = gg_cow,
       height = 5*length(gg_cow_list), width = 15*length(chrName), limitsize = F)


trendPlot <- function(dataFrame, xvar, yvar, xlab, ylab, xtrans, ytrans, xbreaks, ybreaks, xlabels, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) + 
  geom_hex(bins = 60) + 
  scale_x_continuous(trans = xtrans,
                     breaks = xbreaks,
                     labels = xlabels) +
  scale_y_continuous(trans = ytrans,
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
                         digits = 5)) ~
                 "(" * .(genomeBinNamePlot) * " windows)"))
}

#ggTrend_fk_kappa_all_mean_stocha_all <- trendPlot(dataFrame = tab,
#                                                  xvar = fk_kappa_all,
#                                                  yvar = mean_stocha_all,
#                                                  xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
#                                                  ylab = bquote("Mean stochasticity (m"*.(context)*")"),
#                                                  xtrans = log10_trans(),
#                                                  ytrans = log10_trans(),
#                                                  xbreaks = trans_breaks("log10", function(x) 10^x),
#                                                  ybreaks = trans_breaks("log10", function(x) 10^x),
#                                                  xlabels = trans_format("log10", math_format(10^.x)),
#                                                  ylabels = trans_format("log10", math_format(10^.x)))

ggTrend_fk_kappa_all_mean_stocha_all <- trendPlot(dataFrame = tab,
                                                  xvar = fk_kappa_all,
                                                  yvar = mean_stocha_all,
                                                  xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
                                                  ylab = bquote("Mean stochasticity (m"*.(context)*")"),
                                                  xtrans = "identity",
                                                  ytrans = "identity",
                                                  xbreaks = waiver(),
                                                  ybreaks = waiver(),
                                                  xlabels = waiver(),
                                                  ylabels = waiver())
ggTrend_fk_kappa_all_mean_stocha_all <- ggTrend_fk_kappa_all_mean_stocha_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_mean_mean_acf_all <- trendPlot(dataFrame = tab,
                                                    xvar = fk_kappa_all,
                                                    yvar = mean_mean_acf_all,
                                                    xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
                                                    ylab = bquote("Mean mean ACF (m"*.(context)*")"),
                                                    xtrans = "identity",
                                                    ytrans = "identity",
                                                    xbreaks = waiver(),
                                                    ybreaks = waiver(),
                                                    xlabels = waiver(),
                                                    ylabels = waiver())
ggTrend_fk_kappa_all_mean_mean_acf_all <- ggTrend_fk_kappa_all_mean_mean_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_mean_min_acf_all <- trendPlot(dataFrame = tab,
                                                    xvar = fk_kappa_all,
                                                    yvar = mean_min_acf_all,
                                                    xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
                                                    ylab = bquote("Mean min. ACF (m"*.(context)*")"),
                                                    xtrans = "identity",
                                                    ytrans = "identity",
                                                    xbreaks = waiver(),
                                                    ybreaks = waiver(),
                                                    xlabels = waiver(),
                                                    ylabels = waiver())
ggTrend_fk_kappa_all_mean_min_acf_all <- ggTrend_fk_kappa_all_mean_min_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_mean_mean_acf_all <- trendPlot(dataFrame = tab,
                                                    xvar = mean_stocha_all,
                                                    yvar = mean_mean_acf_all,
                                                    xlab = bquote("Mean stochasticity (m"*.(context)*")"),
                                                    ylab = bquote("Mean mean ACF (m"*.(context)*")"),
                                                    xtrans = "identity",
                                                    ytrans = "identity",
                                                    xbreaks = waiver(),
                                                    ybreaks = waiver(),
                                                    xlabels = waiver(),
                                                    ylabels = waiver())
ggTrend_mean_stocha_all_mean_mean_acf_all <- ggTrend_mean_stocha_all_mean_mean_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_mean_min_acf_all <- trendPlot(dataFrame = tab,
                                                    xvar = mean_stocha_all,
                                                    yvar = mean_min_acf_all,
                                                    xlab = bquote("Mean stochasticity (m"*.(context)*")"),
                                                    ylab = bquote("Mean min. ACF (m"*.(context)*")"),
                                                    xtrans = "identity",
                                                    ytrans = "identity",
                                                    xbreaks = waiver(),
                                                    ybreaks = waiver(),
                                                    xlabels = waiver(),
                                                    ylabels = waiver())
ggTrend_mean_stocha_all_mean_min_acf_all <- ggTrend_mean_stocha_all_mean_min_acf_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_genes <- trendPlot(dataFrame = cbind(tab, genes[,3:4]),
                                        xvar = fk_kappa_all,
                                        yvar = features,
                                        xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
                                        ylab = bquote("Genes"),
                                        xtrans = "identity",
                                        ytrans = "identity",
                                        xbreaks = waiver(),
                                        ybreaks = waiver(),
                                        xlabels = waiver(),
                                        ylabels = waiver())
ggTrend_fk_kappa_all_genes <- ggTrend_fk_kappa_all_genes +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_gypsy <- trendPlot(dataFrame = cbind(tab, gypsy[,3:4]),
                                        xvar = fk_kappa_all,
                                        yvar = features,
                                        xlab = bquote("Fleiss' kappa (m"*.(context)*")"),
                                        ylab = bquote(italic("GYPSY")),
                                        xtrans = "identity",
                                        ytrans = "identity",
                                        xbreaks = waiver(),
                                        ybreaks = waiver(),
                                        xlabels = waiver(),
                                        ylabels = waiver())
ggTrend_fk_kappa_all_gypsy <- ggTrend_fk_kappa_all_gypsy +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_genes <- trendPlot(dataFrame = cbind(tab, genes[,3:4]),
                                           xvar = mean_stocha_all,
                                           yvar = features,
                                           xlab = bquote("Mean stochasticity (m"*.(context)*")"),
                                           ylab = bquote("Genes"),
                                           xtrans = "identity",
                                           ytrans = "identity",
                                           xbreaks = waiver(),
                                           ybreaks = waiver(),
                                           xlabels = waiver(),
                                           ylabels = waiver())
ggTrend_mean_stocha_all_genes <- ggTrend_mean_stocha_all_genes +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_gypsy <- trendPlot(dataFrame = cbind(tab, gypsy[,3:4]),
                                           xvar = mean_stocha_all,
                                           yvar = features,
                                           xlab = bquote("Mean stochasticity (m"*.(context)*")"),
                                           ylab = bquote(italic("GYPSY")),
                                           xtrans = "identity",
                                           ytrans = "identity",
                                           xbreaks = waiver(),
                                           ybreaks = waiver(),
                                           xlabels = waiver(),
                                           ylabels = waiver())
ggTrend_mean_stocha_all_gypsy <- ggTrend_mean_stocha_all_gypsy +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mean_acf_all_genes <- trendPlot(dataFrame = cbind(tab, genes[,3:4]),
                                           xvar = mean_mean_acf_all,
                                           yvar = features,
                                           xlab = bquote("Mean mean ACF (m"*.(context)*")"),
                                           ylab = bquote("Genes"),
                                           xtrans = "identity",
                                           ytrans = "identity",
                                           xbreaks = waiver(),
                                           ybreaks = waiver(),
                                           xlabels = waiver(),
                                           ylabels = waiver())
ggTrend_mean_mean_acf_all_genes <- ggTrend_mean_mean_acf_all_genes +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mean_acf_all_gypsy <- trendPlot(dataFrame = cbind(tab, gypsy[,3:4]),
                                           xvar = mean_mean_acf_all,
                                           yvar = features,
                                           xlab = bquote("Mean mean ACF (m"*.(context)*")"),
                                           ylab = bquote(italic("GYPSY")),
                                           xtrans = "identity",
                                           ytrans = "identity",
                                           xbreaks = waiver(),
                                           ybreaks = waiver(),
                                           xlabels = waiver(),
                                           ylabels = waiver())
ggTrend_mean_mean_acf_all_gypsy <- ggTrend_mean_mean_acf_all_gypsy +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_min_acf_all_genes <- trendPlot(dataFrame = cbind(tab, genes[,3:4]),
                                            xvar = mean_min_acf_all,
                                            yvar = features,
                                            xlab = bquote("Mean min. ACF (m"*.(context)*")"),
                                            ylab = bquote("Genes"),
                                            xtrans = "identity",
                                            ytrans = "identity",
                                            xbreaks = waiver(),
                                            ybreaks = waiver(),
                                            xlabels = waiver(),
                                            ylabels = waiver())
ggTrend_mean_min_acf_all_genes <- ggTrend_mean_min_acf_all_genes +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_min_acf_all_gypsy <- trendPlot(dataFrame = cbind(tab, gypsy[,3:4]),
                                            xvar = mean_min_acf_all,
                                            yvar = features,
                                            xlab = bquote("Mean min. ACF (m"*.(context)*")"),
                                            ylab = bquote(italic("GYPSY")),
                                            xtrans = "identity",
                                            ytrans = "identity",
                                            xbreaks = waiver(),
                                            ybreaks = waiver(),
                                            xlabels = waiver(),
                                            ylabels = waiver())
ggTrend_mean_min_acf_all_gypsy <- ggTrend_mean_min_acf_all_gypsy +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list2 <- list(
                    ggTrend_fk_kappa_all_mean_stocha_all, ggTrend_fk_kappa_all_mean_mean_acf_all, ggTrend_fk_kappa_all_mean_min_acf_all,
                    ggTrend_mean_stocha_all_mean_mean_acf_all, ggTrend_mean_stocha_all_mean_min_acf_all,
                    ggTrend_fk_kappa_all_genes, ggTrend_fk_kappa_all_gypsy,
                    ggTrend_mean_stocha_all_genes, ggTrend_mean_stocha_all_gypsy,
                    ggTrend_mean_mean_acf_all_genes, ggTrend_mean_mean_acf_all_gypsy,
                    ggTrend_mean_min_acf_all_genes, ggTrend_mean_min_acf_all_gypsy
                   )
 
gg_cow2 <- plot_grid(plotlist = gg_cow_list2,
                     labels = "AUTO", label_size = 30,
                     align = "hv",
                     axis = "l",
                     nrow = length(gg_cow_list2), ncol = 1)

ggsave(paste0(plotDir,
              sampleName, "_MappedOn_", refbase, "_", context,
              "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
              "_NAmax", NAmax, "_all_trendPlot_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow2,
       height = 5*length(gg_cow_list2), width = 5*length(chrName), limitsize = F)

#ggsave(paste0(plotDir,
#              sampleName, "_MappedOn_", refbase, "_", context,
#              "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#              "_NAmax", NAmax, "_test_", paste0(chrName, collapse = "_"),
#              ".pdf"),
#       plot = ggTrend_fk_kappa_all_genes,
#       height = 5, width = 5*length(chrName), limitsize = F)

