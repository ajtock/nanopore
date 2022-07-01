#!/usr/bin/env Rscript

# Analysis:
# Group genes into quantiles
# based on among-read agreement or site-to-site variability in DNA methylation patterns,
# and mean DNA methylation

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_quantiles.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions' 0.975 0.25
# conda deactivate

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "regions"
#topQthresh <- 0.975
#botQthresh <- 0.25

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]
featRegion <- args[8]
topQthresh <- as.numeric(args[9])
botQthresh <- as.numeric(args[10])

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(stringr)
library(data.table)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(fpc)
library(cluster)
library(ggfortify)
library(Rtsne)
library(umap)
library(RColorBrewer)
library(colorspace)
library(viridis)

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
featDF_filt <- read.table(paste0(outDir,
                                 featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                                 "_", context,
                                 "_NAmax", NAmax,
                                 "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                                 paste0(chrName, collapse = "_"), ".tsv"),
                          header = T)
featDF_filt$kappa_C_density <- featDF_filt$fk_Cs_all / ( (featDF_filt$end - featDF_filt$start + 1) / 1e3)
featDF_filt$stocha_C_density <- featDF_filt$stocha_Cs_all / ( (featDF_filt$end - featDF_filt$start + 1) / 1e3)
featDF_filt$parent <- sub(pattern = "\\.\\d+", replacement = "", x = featDF_filt$name) 
featDF_filt$parent <- sub(pattern = "_\\d+", replacement = "", x = featDF_filt$parent) 

colnames(featDF_filt)[ which(colnames(featDF_filt) %in%
  c("mean_mC_all", "fk_kappa_all", "ka_alpha_all", "mean_stocha_all")
) ] <- c("mC", "fkAgreement", "kaAgreement", "Stochasticity")

featDF_filt <- data.frame(featDF_filt,
                          mC_quantile = "",
                          fkAgreement_quantile = "",
                          kaAgreement_quantile = "",
                          Stochasticity_quantile = "")

featDF_filt[ which(!is.na(featDF_filt[,which(colnames(featDF_filt) == "mC")]) &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) <=
                   1 &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) >
                   topQthresh), ]$mC_quantile <- "Quantile 4"
featDF_filt[ which(!is.na(featDF_filt[,which(colnames(featDF_filt) == "mC")]) &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) <=
                   topQthresh &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) >
                   0.50), ]$mC_quantile <- "Quantile 3"
featDF_filt[ which(!is.na(featDF_filt[,which(colnames(featDF_filt) == "mC")]) &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) <=
                   0.50 &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) >
                   botQthresh), ]$mC_quantile <- "Quantile 2"
featDF_filt[ which(!is.na(featDF_filt[,which(colnames(featDF_filt) == "mC")]) &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) <=
                   botQthresh &
                   rank(featDF_filt[,which(colnames(featDF_filt) == "mC")]) /
                   length(featDF_filt[,which(colnames(featDF_filt) == "mC")]) >=
                   0.0), ]$mC_quantile <- "Quantile 1"

mC_quantiles <- sort(unique(featDF_filt$mC_quantile))

featDF_filt_mC_quantiles <- data.frame() 
for(i in mC_quantiles) {
  featDF_filt_mC_quantile_i <- featDF_filt[ which(featDF_filt$mC_quantile == i), ]

  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) <=
                                   1 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) >
                                   0.50), ]$fkAgreement_quantile <- paste0(i, " upper")
  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) <=
                                   0.50 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "fkAgreement")]) >=
                                   0), ]$fkAgreement_quantile <- paste0(i, " lower")

  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) <=
                                   1 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) >
                                   0.50), ]$kaAgreement_quantile <- paste0(i, " upper")
  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) <=
                                   0.50 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "kaAgreement")]) >=
                                   0), ]$kaAgreement_quantile <- paste0(i, " lower")

  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) <=
                                   1 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) >
                                   0.50), ]$Stochasticity_quantile <- paste0(i, " upper")
  featDF_filt_mC_quantile_i[ which(!is.na(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) <=
                                   0.50 &
                                   rank(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) /
                                   length(featDF_filt_mC_quantile_i[,which(colnames(featDF_filt_mC_quantile_i) == "Stochasticity")]) >=
                                   0), ]$Stochasticity_quantile <- paste0(i, " lower")

  featDF_filt_mC_quantiles <- rbind(featDF_filt_mC_quantiles, featDF_filt_mC_quantile_i)
}


# Dimension reduction:

# PCA of fkAgreement and mean methylation for the given context

fkAgreement_mat_filt_pca_n_dim <- 2

fkAgreement_mat_filt_pca <- prcomp(featDF_filt_mC_quantiles[, which(colnames(featDF_filt_mC_quantiles) %in% c("mC", "fkAgreement")), drop = F], center = T, scale = T)
fkAgreement_mat_filt_pca_summ <- summary(fkAgreement_mat_filt_pca)
fkAgreement_mat_filt_pca_PC1_varexp <- round(fkAgreement_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
fkAgreement_mat_filt_pca_PC2_varexp <- round(fkAgreement_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
fkAgreement_mat_filt_pca_dim <- fkAgreement_mat_filt_pca$x[, seq_len(fkAgreement_mat_filt_pca_n_dim)]
colnames(fkAgreement_mat_filt_pca_dim) <- c("PC1", "PC2")
head(fkAgreement_mat_filt_pca_dim)

stopifnot(nrow(fkAgreement_mat_filt_pca_dim) == length(featDF_filt_mC_quantiles$fkAgreement_quantile))
fkAgreement_mat_filt_pca_dim <- cbind(as.data.frame(fkAgreement_mat_filt_pca_dim),
                                      featDF_filt_mC_quantiles,
                                      type = "PCA")
head(fkAgreement_mat_filt_pca_dim)

#fkAgreement_mat_filt_pam_cl_colours <- rainbow(fkAgreement_mat_filt_pamk_n_cl) 
#fkAgreement_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = fkAgreement_mat_filt_pamk_n_cl)

fkAgreement_mat_filt_pca_loadings <- data.frame(variables = rownames(fkAgreement_mat_filt_pca$rotation),
                                                fkAgreement_mat_filt_pca$rotation)

gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                               mapping = aes(x = PC1, y = PC2, colour = fkAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                               mapping = aes(x = PC1, y = PC2, colour = kaAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                                 mapping = aes(x = PC1, y = PC2, colour = Stochasticity_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_fkAgreement_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = fkAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "fkAgreement",
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = kaAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
  labs(colour = "kaAgreement",
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                                        mapping = aes(x = PC1, y = PC2, colour = Stochasticity)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "plasma") +
  labs(colour = "Stochasticity",
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_mC_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
                                             mapping = aes(x = PC1, y = PC2, colour = mC)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = bquote("m" * .(context) ~ "mean"),
       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_fkAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_fkAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_fkAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_fkAgreement_mat_filt_pca_dim_chr <- gg_mC_fkAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_fkAgreement_mat_filt_pca_dim_chr <- list(
                                                     gg_fkAgreement_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_kaAgreement_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_mC_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim_chr
                                                    )
gg_cow_fkAgreement_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_fkAgreement_mat_filt_pca_dim_chr,
                                                   align = "hv",
                                                   axis = "l",
                                                   nrow = length(gg_cow_list_fkAgreement_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_fkAgreement_mat_filt_pca_dim_chr_quantiles.pdf"),
       plot = gg_cow_fkAgreement_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_fkAgreement_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_fkAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_fkAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_fkAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_mC_fkAgreement_mat_filt_pca_dim_loadings <- gg_mC_fkAgreement_mat_filt_pca_dim +
  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings <- list(
                                                          gg_fkAgreement_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_kaAgreement_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Stochasticity_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_mC_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_fkAgreement_quantile_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_kaAgreement_quantile_fkAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Stochasticity_quantile_fkAgreement_mat_filt_pca_dim_loadings
                                                         )
gg_cow_fkAgreement_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings,
                                                          align = "hv",
                                                          axis = "l",
                                                          nrow = 1, ncol = length(gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_fkAgreement_mat_filt_pca_dim_loadings_quantiles.pdf"),
       plot = gg_cow_fkAgreement_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings), limitsize = F)



# PCA of kaAgreement and mean methylation for the given context

kaAgreement_mat_filt_pca_n_dim <- 2

kaAgreement_mat_filt_pca <- prcomp(featDF_filt_mC_quantiles[, which(colnames(featDF_filt_mC_quantiles) %in% c("mC", "kaAgreement")), drop = F], center = T, scale = T)
kaAgreement_mat_filt_pca_summ <- summary(kaAgreement_mat_filt_pca)
kaAgreement_mat_filt_pca_PC1_varexp <- round(kaAgreement_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
kaAgreement_mat_filt_pca_PC2_varexp <- round(kaAgreement_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
kaAgreement_mat_filt_pca_dim <- kaAgreement_mat_filt_pca$x[, seq_len(kaAgreement_mat_filt_pca_n_dim)]
colnames(kaAgreement_mat_filt_pca_dim) <- c("PC1", "PC2")
head(kaAgreement_mat_filt_pca_dim)

stopifnot(nrow(kaAgreement_mat_filt_pca_dim) == length(featDF_filt_mC_quantiles$kaAgreement_quantile))
kaAgreement_mat_filt_pca_dim <- cbind(as.data.frame(kaAgreement_mat_filt_pca_dim),
                                      featDF_filt_mC_quantiles,
                                      type = "PCA")
head(kaAgreement_mat_filt_pca_dim)

#kaAgreement_mat_filt_pam_cl_colours <- rainbow(kaAgreement_mat_filt_pamk_n_cl) 
#kaAgreement_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = kaAgreement_mat_filt_pamk_n_cl)

kaAgreement_mat_filt_pca_loadings <- data.frame(variables = rownames(kaAgreement_mat_filt_pca$rotation),
                                                kaAgreement_mat_filt_pca$rotation)

gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                               mapping = aes(x = PC1, y = PC2, colour = fkAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                               mapping = aes(x = PC1, y = PC2, colour = kaAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                                 mapping = aes(x = PC1, y = PC2, colour = Stochasticity_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_fkAgreement_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = fkAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "fkAgreement",
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = kaAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "kaAgreement",
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                        mapping = aes(x = PC1, y = PC2, colour = Stochasticity)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "plasma") +
  labs(colour = "Stochasticity",
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_mC_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                             mapping = aes(x = PC1, y = PC2, colour = mC)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = bquote("m" * .(context) ~ "mean"),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_kaAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_kaAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_kaAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_kaAgreement_mat_filt_pca_dim_chr <- gg_mC_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_kaAgreement_mat_filt_pca_dim_chr <- list(
                                                     gg_kaAgreement_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_fkAgreement_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_mC_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim_chr
                                                    )
gg_cow_kaAgreement_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_kaAgreement_mat_filt_pca_dim_chr,
                                                     align = "hv",
                                                     axis = "l",
                                                     nrow = length(gg_cow_list_kaAgreement_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_kaAgreement_mat_filt_pca_dim_chr_quantiles.pdf"),
       plot = gg_cow_kaAgreement_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_kaAgreement_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_kaAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_kaAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_kaAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_mC_kaAgreement_mat_filt_pca_dim_loadings <- gg_mC_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings <- list(
                                                          gg_kaAgreement_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_fkAgreement_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Stochasticity_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_mC_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_kaAgreement_quantile_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_fkAgreement_quantile_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Stochasticity_quantile_kaAgreement_mat_filt_pca_dim_loadings
                                                         )
gg_cow_kaAgreement_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings,
                                                          align = "hv",
                                                          axis = "l",
                                                          nrow = 1, ncol = length(gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_kaAgreement_mat_filt_pca_dim_loadings_quantiles.pdf"),
       plot = gg_cow_kaAgreement_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings), limitsize = F)



# PCA of Stochasticity and mean methylation for the given context

Stochasticity_mat_filt_pca_n_dim <- 2

Stochasticity_mat_filt_pca <- prcomp(featDF_filt_mC_quantiles[, which(colnames(featDF_filt_mC_quantiles) %in% c("mC", "Stochasticity")), drop = F], center = T, scale = T)
Stochasticity_mat_filt_pca_summ <- summary(Stochasticity_mat_filt_pca)
Stochasticity_mat_filt_pca_PC1_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
Stochasticity_mat_filt_pca_PC2_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
Stochasticity_mat_filt_pca_dim <- Stochasticity_mat_filt_pca$x[, seq_len(Stochasticity_mat_filt_pca_n_dim)]
colnames(Stochasticity_mat_filt_pca_dim) <- c("PC1", "PC2")
head(Stochasticity_mat_filt_pca_dim)

stopifnot(nrow(Stochasticity_mat_filt_pca_dim) == length(featDF_filt_mC_quantiles$Stochasticity_quantile))
Stochasticity_mat_filt_pca_dim <- cbind(as.data.frame(Stochasticity_mat_filt_pca_dim),
                                        featDF_filt_mC_quantiles,
                                        type = "PCA")
head(Stochasticity_mat_filt_pca_dim)

#Stochasticity_mat_filt_pam_cl_colours <- rainbow(Stochasticity_mat_filt_pamk_n_cl) 
#Stochasticity_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = Stochasticity_mat_filt_pamk_n_cl)

Stochasticity_mat_filt_pca_loadings <- data.frame(variables = rownames(Stochasticity_mat_filt_pca$rotation),
                                                  Stochasticity_mat_filt_pca$rotation)

gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                mapping = aes(x = PC1, y = PC2, colour = fkAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                mapping = aes(x = PC1, y = PC2, colour = kaAgreement_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                   mapping = aes(x = PC1, y = PC2, colour = Stochasticity_quantile)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_fkAgreement_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                        mapping = aes(x = PC1, y = PC2, colour = fkAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "fkAgreement",
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                        mapping = aes(x = PC1, y = PC2, colour = kaAgreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "kaAgreement",
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                          mapping = aes(x = PC1, y = PC2, colour = Stochasticity)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "plasma") +
  labs(colour = "Stochasticity",
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_mC_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                               mapping = aes(x = PC1, y = PC2, colour = mC)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = bquote("m" * .(context) ~ "mean"),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim_chr <- gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim_chr <- gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_Stochasticity_mat_filt_pca_dim_chr <- gg_fkAgreement_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_Stochasticity_mat_filt_pca_dim_chr <- gg_kaAgreement_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_Stochasticity_mat_filt_pca_dim_chr <- gg_mC_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_Stochasticity_mat_filt_pca_dim_chr <- list(
                                                       gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_fkAgreement_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_kaAgreement_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_mC_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim_chr
                                                      )
gg_cow_Stochasticity_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_chr,
                                                       align = "hv",
                                                       axis = "l",
                                                       nrow = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_chr_quantiles.pdf"),
       plot = gg_cow_Stochasticity_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim_loadings <- gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim_loadings <- gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim_loadings <- gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_Stochasticity_mat_filt_pca_dim_loadings <- gg_Stochasticity_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_Stochasticity_mat_filt_pca_dim_loadings <- gg_fkAgreement_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_Stochasticity_mat_filt_pca_dim_loadings <- gg_kaAgreement_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_mC_Stochasticity_mat_filt_pca_dim_loadings <- gg_mC_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings <- list(
                                                            gg_Stochasticity_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_fkAgreement_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_kaAgreement_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_mC_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_Stochasticity_quantile_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_fkAgreement_quantile_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_kaAgreement_quantile_Stochasticity_mat_filt_pca_dim_loadings
                                                           )
gg_cow_Stochasticity_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings,
                                                            align = "hv",
                                                            axis = "l",
                                                            nrow = 1, ncol = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_loadings_quantiles.pdf"),
       plot = gg_cow_Stochasticity_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings), limitsize = F)



# Append intron retention ratio (calculated with IRFinder)
Col_Rep1_IRFinder <- fread(paste0("/home/ajt200/analysis/RNAseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_RNAseq_IRFinder_TAIR10_chr_all/REF/TAIR10_chr_all/",
                                  "Col_0_RNAseq_pooled_ERR96615/IRFinder-IR-dir.txt"),
                           sep = "\t", data.table = F)
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[grep("clean", Col_Rep1_IRFinder$Name),]
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[-which(Col_Rep1_IRFinder$Warnings %in% c("LowCover")),]
#nrow(Col_Rep1_IRFinder[which(Col_Rep1_IRFinder$Warnings == "-"),])
#[1] 45907
#[1] 22136
Col_Rep1_IRFinder$Name <- str_extract(Col_Rep1_IRFinder$Name, "AT\\wG\\d+")
Col_Rep1_IRFinder <- Col_Rep1_IRFinder[-which(is.na(Col_Rep1_IRFinder$Name)),]

library(doParallel)
library(doFuture)
registerDoFuture()
plan(multicore)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

parentIDs <- unique(Col_Rep1_IRFinder$Name)

Col_Rep1_IRratio <- foreach(i = iter(parentIDs),
                            .combine = "rbind",
                            .multicombine = T,
                            .maxcombine = length(parentIDs)+1e1,
                            .inorder = F,
                            .errorhandling = "pass") %dopar% {
  tmpDF <- Col_Rep1_IRFinder[which(Col_Rep1_IRFinder$Name == i),]
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

featDF_filt_tab <- base::merge(x = featDF_filt_mC_quantiles, y = Col_Rep1_IRratio,
                               by.x = "parent", by.y = "parent")

colnames(featDF_filt_tab)[which(colnames(featDF_filt_tab) %in%
  c("mC", "fkAgreement", "kaAgreement", "Stochasticity"))] <- c("mean_mC_all", "fk_kappa_all", "ka_alpha_all", "mean_stocha_all")

colnames(featDF_filt_mC_quantiles)[which(colnames(featDF_filt_mC_quantiles) %in%
  c("mC", "fkAgreement", "kaAgreement", "Stochasticity"))] <- c("mean_mC_all", "fk_kappa_all", "ka_alpha_all", "mean_stocha_all")


print(cor.test(featDF_filt_tab$fk_kappa_all, featDF_filt_tab$IRratio_mean, method = "spearman"))
#-0.1034838
print(cor.test(featDF_filt_tab$fk_kappa_all, featDF_filt_tab$IRratio_median, method = "spearman"))
#-0.2840318

print(cor.test(featDF_filt_tab$ka_alpha_all, featDF_filt_tab$IRratio_mean, method = "spearman"))
#-0.1052147
print(cor.test(featDF_filt_tab$ka_alpha_all, featDF_filt_tab$IRratio_median, method = "spearman"))
#-0.2763196

print(cor.test(featDF_filt_tab$mean_stocha_all, featDF_filt_tab$IRratio_mean, method = "spearman"))
#-0.1120429
print(cor.test(featDF_filt_tab$mean_stocha_all, featDF_filt_tab$IRratio_median, method = "spearman"))
#-0.3074874


# Plot relationships and define groups
trendPlot <- function(dataFrame, mapping, paletteName, xvar, yvar, quantilelab, xlab, ylab, xaxtrans, yaxtrans, xbreaks, ybreaks, xlabels, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = mapping) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_x_continuous(trans = xaxtrans,
                     breaks = xbreaks,
                     labels = xlabels) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  scale_colour_brewer(palette = paletteName) +
  geom_smooth(colour = "black", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y ~ s(x, bs = "cs")) +
  labs(colour = quantilelab,
       x = xlab,
       y = ylab) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        panel.grid = element_blank(),
        legend.key.height = unit(4, "mm"),
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


ggTrend_mean_mC_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                   mapping = aes(x = mean_mC_all, y = fk_kappa_all, colour = fkAgreement_quantile),
                                                   paletteName = "Dark2",
                                                   xvar = mean_mC_all,
                                                   yvar = fk_kappa_all,
                                                   quantilelab = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
                                                   xlab = bquote(.(featName)*" mean m"*.(context)),
                                                   ylab = bquote(.(featName)*" fkAgreement (m"*.(context)*")"),
                                                   xaxtrans = log10_trans(),
                                                   yaxtrans = log10_trans(),
                                                   xbreaks = trans_breaks("log10", function(x) 10^x),
                                                   ybreaks = trans_breaks("log10", function(x) 10^x),
                                                   xlabels = trans_format("log10", math_format(10^.x)),
                                                   ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_fk_kappa_all_filt <- ggTrend_mean_mC_all_fk_kappa_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                   mapping = aes(x = mean_mC_all, y = ka_alpha_all, colour = kaAgreement_quantile),
                                                   paletteName = "Dark2",
                                                   xvar = mean_mC_all,
                                                   yvar = ka_alpha_all,
                                                   quantilelab = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
                                                   xlab = bquote(.(featName)*" mean m"*.(context)),
                                                   ylab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                   xaxtrans = log10_trans(),
                                                   yaxtrans = log10_trans(),
                                                   xbreaks = trans_breaks("log10", function(x) 10^x),
                                                   ybreaks = trans_breaks("log10", function(x) 10^x),
                                                   xlabels = trans_format("log10", math_format(10^.x)),
                                                   ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_ka_alpha_all_filt <- ggTrend_mean_mC_all_ka_alpha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                      mapping = aes(x = mean_mC_all, y = mean_stocha_all, colour = Stochasticity_quantile),
                                                      paletteName = "Set1",
                                                      xvar = mean_mC_all,
                                                      yvar = mean_stocha_all,
                                                      quantilelab = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
                                                      xlab = bquote(.(featName)*" mean m"*.(context)),
                                                      ylab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                      xaxtrans = log10_trans(),
                                                      yaxtrans = log10_trans(),
                                                      xbreaks = trans_breaks("log10", function(x) 10^x),
                                                      ybreaks = trans_breaks("log10", function(x) 10^x),
                                                      xlabels = trans_format("log10", math_format(10^.x)),
                                                      ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_mean_stocha_all_filt <- ggTrend_mean_mC_all_mean_stocha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_ka_alpha_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                    mapping = aes(x = ka_alpha_all, y = fk_kappa_all, colour = fkAgreement_quantile),
                                                    paletteName = "Dark2",
                                                    xvar = ka_alpha_all,
                                                    yvar = fk_kappa_all,
                                                    quantilelab = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
                                                    xlab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                    ylab = bquote(.(featName)*" fkAgreement (m"*.(context)*")"),
                                                    xaxtrans = log10_trans(),
                                                    yaxtrans = log10_trans(),
                                                    xbreaks = trans_breaks("log10", function(x) 10^x),
                                                    ybreaks = trans_breaks("log10", function(x) 10^x),
                                                    xlabels = trans_format("log10", math_format(10^.x)),
                                                    ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_ka_alpha_all_fk_kappa_all_filt <- ggTrend_ka_alpha_all_fk_kappa_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                    mapping = aes(x = ka_alpha_all, y = fk_kappa_all, colour = kaAgreement_quantile),
                                                    paletteName = "Dark2",
                                                    xvar = ka_alpha_all,
                                                    yvar = fk_kappa_all,
                                                    quantilelab = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
                                                    xlab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                    ylab = bquote(.(featName)*" fkAgreement (m"*.(context)*")"),
                                                    xaxtrans = log10_trans(),
                                                    yaxtrans = log10_trans(),
                                                    xbreaks = trans_breaks("log10", function(x) 10^x),
                                                    ybreaks = trans_breaks("log10", function(x) 10^x),
                                                    xlabels = trans_format("log10", math_format(10^.x)),
                                                    ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_kappa_all_ka_alpha_all_filt <- ggTrend_fk_kappa_all_ka_alpha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = fkAgreement_quantile),
                                                       paletteName = "Dark2",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       quantilelab = bquote(atop("fkAgreement &", "m"*.(context)~"quantile")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" fkAgreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_stocha_all_fk_kappa_all_filt <- ggTrend_mean_stocha_all_fk_kappa_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = Stochasticity_quantile),
                                                       paletteName = "Set1",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       quantilelab = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" fkAgreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_kappa_all_mean_stocha_all_filt <- ggTrend_fk_kappa_all_mean_stocha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_stocha_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                       mapping = aes(x = mean_stocha_all, y = ka_alpha_all, colour = fkAgreement_quantile),
                                                       paletteName = "Dark2",
                                                       xvar = mean_stocha_all,
                                                       yvar = ka_alpha_all,
                                                       quantilelab = bquote(atop("kaAgreement &", "m"*.(context)~"quantile")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_stocha_all_ka_alpha_all_filt <- ggTrend_mean_stocha_all_ka_alpha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_ka_alpha_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt_mC_quantiles,
                                                       mapping = aes(x = mean_stocha_all, y = ka_alpha_all, colour = Stochasticity_quantile),
                                                       paletteName = "Set1",
                                                       xvar = mean_stocha_all,
                                                       yvar = ka_alpha_all,
                                                       quantilelab = bquote(atop("Stochasticity &", "m"*.(context)~"quantile")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_ka_alpha_all_mean_stocha_all_filt <- ggTrend_ka_alpha_all_mean_stocha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list1 <- list(
                     ggTrend_mean_mC_all_fk_kappa_all_filt,
                     ggTrend_mean_mC_all_ka_alpha_all_filt,
                     ggTrend_mean_mC_all_mean_stocha_all_filt,
                     ggTrend_ka_alpha_all_fk_kappa_all_filt,
                     ggTrend_fk_kappa_all_ka_alpha_all_filt,
                     ggTrend_mean_stocha_all_fk_kappa_all_filt,
                     ggTrend_fk_kappa_all_mean_stocha_all_filt,
                     ggTrend_mean_stocha_all_ka_alpha_all_filt,
                     ggTrend_ka_alpha_all_mean_stocha_all_filt
                    )

gg_cow1 <- plot_grid(plotlist = gg_cow_list1,
#                     labels = "AUTO", label_size = 30,
                     align = "hv",
                     axis = "l",
                     nrow = length(gg_cow_list1), ncol = 1)

ggsave(paste0(plotDir,
              featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase, "_", context,
              "_NAmax", NAmax, "_all_trendPlot_", paste0(chrName, collapse = "_"),
              "_quantiles.pdf"),
       plot = gg_cow1,
       height = 5*length(gg_cow_list1), width = 5*length(chrName), limitsize = F)



# Extract feature quantiles to enable enrichment analysis


# Filter by fk_kappa_all and mean_mC_all quantile
featDF_filt_kappa_mC_quantile1_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 1 upper")
featDF_filt_kappa_mC_quantile2_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 2 upper")
featDF_filt_kappa_mC_quantile3_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 3 upper")
featDF_filt_kappa_mC_quantile4_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 4 upper")

write.table(featDF_filt_kappa_mC_quantile1_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile1_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile2_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile2_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile3_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile3_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile4_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile4_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

featDF_filt_kappa_mC_quantile1_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 1 lower")
featDF_filt_kappa_mC_quantile2_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 2 lower")
featDF_filt_kappa_mC_quantile3_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 3 lower")
featDF_filt_kappa_mC_quantile4_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(fkAgreement_quantile == "Quantile 4 lower")

write.table(featDF_filt_kappa_mC_quantile1_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile1_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile2_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile2_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile3_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile3_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_quantile4_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_quantile4_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by ka_alpha_all and mean_mC_all quantile
featDF_filt_alpha_mC_quantile1_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 1 upper")
featDF_filt_alpha_mC_quantile2_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 2 upper")
featDF_filt_alpha_mC_quantile3_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 3 upper")
featDF_filt_alpha_mC_quantile4_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 4 upper")

write.table(featDF_filt_alpha_mC_quantile1_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile1_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile2_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile2_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile3_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile3_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile4_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile4_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

featDF_filt_alpha_mC_quantile1_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 1 lower")
featDF_filt_alpha_mC_quantile2_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 2 lower")
featDF_filt_alpha_mC_quantile3_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 3 lower")
featDF_filt_alpha_mC_quantile4_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(kaAgreement_quantile == "Quantile 4 lower")

write.table(featDF_filt_alpha_mC_quantile1_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile1_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile2_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile2_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile3_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile3_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_quantile4_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_quantile4_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all quantile
featDF_filt_stocha_mC_quantile1_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 1 upper")
featDF_filt_stocha_mC_quantile2_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 2 upper")
featDF_filt_stocha_mC_quantile3_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 3 upper")
featDF_filt_stocha_mC_quantile4_upper <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 4 upper")

write.table(featDF_filt_stocha_mC_quantile1_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile1_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile2_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile2_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile3_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile3_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile4_upper,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile4_upper_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

featDF_filt_stocha_mC_quantile1_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 1 lower")
featDF_filt_stocha_mC_quantile2_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 2 lower")
featDF_filt_stocha_mC_quantile3_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 3 lower")
featDF_filt_stocha_mC_quantile4_lower <- featDF_filt_mC_quantiles %>%
  dplyr::filter(Stochasticity_quantile == "Quantile 4 lower")

write.table(featDF_filt_stocha_mC_quantile1_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile1_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile2_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile2_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile3_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile3_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_quantile4_lower,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_quantile4_lower_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

