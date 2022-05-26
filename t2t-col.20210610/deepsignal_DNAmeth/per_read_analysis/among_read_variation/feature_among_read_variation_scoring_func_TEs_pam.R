#!/usr/bin/env Rscript

# Analysis:
# Group TEs into clusters using partitioning around medoids (PAM),
# based on among-read agreement or site-to-site variability in DNA methylation patterns,
# and mean DNA methylation

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_TEs_pam.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'TE' 'bodies'
# conda deactivate

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CHG"
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
library(pals)

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
colnames(featDF_filt)[which(colnames(featDF_filt) == "score")] <- "Superfamily"

mat_filt <- featDF_filt[,which(colnames(featDF_filt) %in%
                               c("mean_mC_all", "fk_kappa_all", "ka_alpha_all", "mean_stocha_all")), drop = F]
colnames(mat_filt) <- c("mC", "fkAgreement", "kaAgreement", "Stochasticity")

set.seed(4849345)
fkAgreement_mat_filt_pamk <- fpc::pamk(data = mat_filt[ , c(1,2), drop = F],
                                       krange = 3:10,
                                       criterion = "multiasw",
                                       usepam = F,
                                       scaling = T,
                                       alpha = 0.001,
                                       diss = F,
                                       critout = T,
                                       ns = 10)
fkAgreement_mat_filt_pamk_n_cl <- fkAgreement_mat_filt_pamk$nc
print(fkAgreement_mat_filt_pamk_n_cl)
fkAgreement_mat_filt_pam <- cluster::pam(x = mat_filt[ , c(1,2), drop = F],
                                         k = 4,
                                         diss = F,
                                         metric = "manhattan",
                                         stand = T,
                                         cluster.only = F,
                                         do.swap = T,
                                         pamonce = 0)
 
featDF_filt$fkAgreement_cluster <- paste0("Cluster ", fkAgreement_mat_filt_pam$clustering)
#featDF_filt$fkAgreement_cluster <- fkAgreement_mat_filt_pam$clustering
#featDF_filt$fkAgreement_cluster[which(featDF_filt$fkAgreement_cluster == 2)] <- "Cluster 1"
#featDF_filt$fkAgreement_cluster[which(featDF_filt$fkAgreement_cluster == 1)] <- "Cluster 2"
#featDF_filt$fkAgreement_cluster[which(featDF_filt$fkAgreement_cluster == 3)] <- "Cluster 3"
#featDF_filt$fkAgreement_cluster[which(featDF_filt$fkAgreement_cluster == 4)] <- "Cluster 4"

set.seed(4849345)
kaAgreement_mat_filt_pamk <- fpc::pamk(data = mat_filt[ , c(1,3), drop = F],
                                       krange = 3:10,
                                       criterion = "multiasw",
                                       usepam = F,
                                       scaling = T,
                                       alpha = 0.001,
                                       diss = F,
                                       critout = T,
                                       ns = 10)
kaAgreement_mat_filt_pamk_n_cl <- kaAgreement_mat_filt_pamk$nc
print(kaAgreement_mat_filt_pamk_n_cl)
kaAgreement_mat_filt_pam <- cluster::pam(x = mat_filt[ , c(1,3), drop = F],
                                         k = kaAgreement_mat_filt_pamk_n_cl,
                                         diss = F,
                                         metric = "manhattan",
                                         stand = T,
                                         cluster.only = F,
                                         do.swap = T,
                                         pamonce = 0)
 
featDF_filt$kaAgreement_cluster <- paste0("Cluster ", kaAgreement_mat_filt_pam$clustering)
#featDF_filt$kaAgreement_cluster <- kaAgreement_mat_filt_pam$clustering
#featDF_filt$kaAgreement_cluster[which(featDF_filt$kaAgreement_cluster == 2)] <- "Cluster 1"
#featDF_filt$kaAgreement_cluster[which(featDF_filt$kaAgreement_cluster == 1)] <- "Cluster 2"
#featDF_filt$kaAgreement_cluster[which(featDF_filt$kaAgreement_cluster == 3)] <- "Cluster 3"
#featDF_filt$kaAgreement_cluster[which(featDF_filt$kaAgreement_cluster == 4)] <- "Cluster 4"

set.seed(4849345)
Stochasticity_mat_filt_pamk <- fpc::pamk(data = mat_filt[ , c(1,4), drop = F],
                                         krange = 3:10,
                                         criterion = "multiasw",
                                         usepam = F,
                                         scaling = T,
                                         alpha = 0.001,
                                         diss = F,
                                         critout = T,
                                         ns = 10)
Stochasticity_mat_filt_pamk_n_cl <- Stochasticity_mat_filt_pamk$nc
print(Stochasticity_mat_filt_pamk_n_cl)
Stochasticity_mat_filt_pam <- cluster::pam(x = mat_filt[ , c(1,4), drop = F],
                                           k = 4,
                                           diss = F,
                                           metric = "manhattan",
                                           stand = T,
                                           cluster.only = F,
                                           do.swap = T,
                                           pamonce = 0)

featDF_filt$Stochasticity_cluster <- paste0("Cluster ", Stochasticity_mat_filt_pam$clustering)
#featDF_filt$Stochasticity_cluster <- Stochasticity_mat_filt_pam$clustering
#featDF_filt$Stochasticity_cluster[which(featDF_filt$Stochasticity_cluster == 3)] <- "Cluster 1"
#featDF_filt$Stochasticity_cluster[which(featDF_filt$Stochasticity_cluster == 1)] <- "Cluster 2"
#featDF_filt$Stochasticity_cluster[which(featDF_filt$Stochasticity_cluster == 2)] <- "Cluster 3"
#featDF_filt$Stochasticity_cluster[which(featDF_filt$Stochasticity_cluster == 4)] <- "Cluster 4"


# Dimension reduction:

## PCA of fkAgreement and mean methylation for the given context
#
#fkAgreement_mat_filt_pca_n_dim <- 2
#
##fkAgreement_mat_filt_pca <- prcomp(mat_filt[ , c(1,2), drop = F], center = T, scale = T)
#fkAgreement_mat_filt_pca <- prcomp(fkAgreement_mat_filt_pam$data, center = F, scale = F)
#fkAgreement_mat_filt_pca_summ <- summary(fkAgreement_mat_filt_pca)
#fkAgreement_mat_filt_pca_PC1_varexp <- round(fkAgreement_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
#fkAgreement_mat_filt_pca_PC2_varexp <- round(fkAgreement_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
#fkAgreement_mat_filt_pca_dim <- fkAgreement_mat_filt_pca$x[, seq_len(fkAgreement_mat_filt_pca_n_dim)]
#colnames(fkAgreement_mat_filt_pca_dim) <- c("PC1", "PC2")
#head(fkAgreement_mat_filt_pca_dim)
#
#stopifnot(nrow(fkAgreement_mat_filt_pca_dim) == length(featDF_filt$fkAgreement_cluster))
#fkAgreement_mat_filt_pca_dim <- cbind(as.data.frame(fkAgreement_mat_filt_pca_dim),
#                                      mat_filt,
#                                      chr = featDF_filt$chr,
#                                      fkAgreement_cluster = featDF_filt$fkAgreement_cluster,
#                                      kaAgreement_cluster = featDF_filt$kaAgreement_cluster,
#                                      Stochasticity_cluster = featDF_filt$Stochasticity_cluster,
#                                      Superfamily = featDF_filt$Superfamily,
#                                      type = "PCA")
#head(fkAgreement_mat_filt_pca_dim)
#
##fkAgreement_mat_filt_pam_cl_colours <- rainbow(fkAgreement_mat_filt_pamk_n_cl) 
##fkAgreement_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = fkAgreement_mat_filt_pamk_n_cl)
#
#fkAgreement_mat_filt_pca_loadings <- data.frame(variables = rownames(fkAgreement_mat_filt_pca$rotation),
#                                                fkAgreement_mat_filt_pca$rotation)
#
#ap_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim <- autoplot(fkAgreement_mat_filt_pam,
#                                                                frame = T,
#                                                                frame.type = "norm",
#                                                                loadings = T,
#                                                                loadings.colour = "black",
#                                                                loadings.label = T,
#                                                                loadings.label.colour = "black",
#                                                                label = F)
#aps_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim <- autoplot(silhouette(fkAgreement_mat_filt_pam))
#ggsave(paste0(plotDir, "ap_fkAgreement_mat_filt_pca_dim.pdf"),
#       plot = ap_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim,
#       height = 4, width = 5, limitsize = F)
#ggsave(paste0(plotDir, "aps_fkAgreement_mat_filt_pca_dim.pdf"),
#       plot = aps_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim,
#       height = 4, width = 5, limitsize = F)
#
#gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                              mapping = aes(x = PC1, y = PC2, colour = fkAgreement_cluster)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_brewer(palette = "Dark2") +
##  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
#  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                              mapping = aes(x = PC1, y = PC2, colour = kaAgreement_cluster)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_brewer(palette = "Dark2") +
##  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
#  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                                mapping = aes(x = PC1, y = PC2, colour = Stochasticity_cluster)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_brewer(palette = "Set1") +
##  scale_colour_manual(values = fkAgreement_mat_filt_pam_cl_colours) +
#  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_fkAgreement_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                      mapping = aes(x = PC1, y = PC2, colour = fkAgreement)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_viridis_c(option = "turbo") +
##  scale_colour_gradient(low = "red", high = "yellow") +
##  scale_colour_gradient2(low = "blue", high = "red") +
#  labs(colour = "fkAgreement",
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_kaAgreement_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                      mapping = aes(x = PC1, y = PC2, colour = kaAgreement)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_viridis_c(option = "turbo") +
#  labs(colour = "kaAgreement",
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_Stochasticity_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                                        mapping = aes(x = PC1, y = PC2, colour = Stochasticity)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_viridis_c(option = "plasma") +
#  labs(colour = "Stochasticity",
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
#gg_mC_fkAgreement_mat_filt_pca_dim <- ggplot(fkAgreement_mat_filt_pca_dim,
#                                             mapping = aes(x = PC1, y = PC2, colour = mC)) +
#  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_viridis_c(option = "viridis") +
#  labs(colour = bquote("m" * .(context) ~ "mean"),
#       x = bquote("PC1 (" * .(fkAgreement_mat_filt_pca_PC1_varexp) * "%)"),
#       y = bquote("PC2 (" * .(fkAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        legend.key.height = unit(4, "mm"))
#
## Plot chr
#gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_fkAgreement_fkAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_kaAgreement_fkAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_Stochasticity_fkAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_mC_fkAgreement_mat_filt_pca_dim_chr <- gg_mC_fkAgreement_mat_filt_pca_dim +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_cow_list_fkAgreement_mat_filt_pca_dim_chr <- list(
#                                                     gg_fkAgreement_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_kaAgreement_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_Stochasticity_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_mC_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim_chr,
#                                                     gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim_chr
#                                                    )
#gg_cow_fkAgreement_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_fkAgreement_mat_filt_pca_dim_chr,
#                                                   align = "hv",
#                                                   axis = "l",
#                                                   nrow = length(gg_cow_list_fkAgreement_mat_filt_pca_dim_chr), ncol = 1)
#ggsave(paste0(plotDir, "gg_fkAgreement_mat_filt_pca_dim_chr.pdf"),
#       plot = gg_cow_fkAgreement_mat_filt_pca_dim_chr,
#       height = 4*length(gg_cow_list_fkAgreement_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)
#
#
## Overlay loadings
#gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_fkAgreement_fkAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_kaAgreement_fkAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_Stochasticity_fkAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_mC_fkAgreement_mat_filt_pca_dim_loadings <- gg_mC_fkAgreement_mat_filt_pca_dim +
#  geom_segment(data = fkAgreement_mat_filt_pca_loadings,
#               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
#               arrow = arrow(length = unit(1/2, "picas")),
#               colour = "grey0") +
#  annotate("text", x = (fkAgreement_mat_filt_pca_loadings$PC1*3.2), y = (fkAgreement_mat_filt_pca_loadings$PC2*3.2),
#           label = fkAgreement_mat_filt_pca_loadings$variables, colour = "grey0")
#
#gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings <- list(
#                                                        gg_fkAgreement_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_kaAgreement_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_Stochasticity_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_mC_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_fkAgreement_cluster_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_kaAgreement_cluster_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        gg_Stochasticity_cluster_fkAgreement_mat_filt_pca_dim_loadings
#                                                       )
#gg_cow_fkAgreement_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings,
#                                                        align = "hv",
#                                                        axis = "l",
#                                                        nrow = 1, ncol = length(gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings))
#ggsave(paste0(plotDir, "gg_fkAgreement_mat_filt_pca_dim_loadings.pdf"),
#       plot = gg_cow_fkAgreement_mat_filt_pca_dim_loadings,
#       height = 4, width = 4*length(gg_cow_list_fkAgreement_mat_filt_pca_dim_loadings), limitsize = F)



# PCA of kaAgreement and mean methylation for the given context

kaAgreement_mat_filt_pca_n_dim <- 2

#kaAgreement_mat_filt_pca <- prcomp(mat_filt[ , c(1,3), drop = F], center = T, scale = T)
kaAgreement_mat_filt_pca <- prcomp(kaAgreement_mat_filt_pam$data, center = F, scale = F)
kaAgreement_mat_filt_pca_summ <- summary(kaAgreement_mat_filt_pca)
kaAgreement_mat_filt_pca_PC1_varexp <- round(kaAgreement_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
kaAgreement_mat_filt_pca_PC2_varexp <- round(kaAgreement_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
kaAgreement_mat_filt_pca_dim <- kaAgreement_mat_filt_pca$x[, seq_len(kaAgreement_mat_filt_pca_n_dim)]
colnames(kaAgreement_mat_filt_pca_dim) <- c("PC1", "PC2")
head(kaAgreement_mat_filt_pca_dim)

stopifnot(nrow(kaAgreement_mat_filt_pca_dim) == length(featDF_filt$kaAgreement_cluster))
kaAgreement_mat_filt_pca_dim <- cbind(as.data.frame(kaAgreement_mat_filt_pca_dim),
                                      mat_filt,
                                      chr = featDF_filt$chr,
                                      fkAgreement_cluster = featDF_filt$fkAgreement_cluster,
                                      kaAgreement_cluster = featDF_filt$kaAgreement_cluster,
                                      Stochasticity_cluster = featDF_filt$Stochasticity_cluster,
                                      Superfamily = featDF_filt$Superfamily,
                                      type = "PCA")
head(kaAgreement_mat_filt_pca_dim)

#kaAgreement_mat_filt_pam_cl_colours <- rainbow(kaAgreement_mat_filt_pamk_n_cl) 
#kaAgreement_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = kaAgreement_mat_filt_pamk_n_cl)

kaAgreement_mat_filt_pca_loadings <- data.frame(variables = rownames(kaAgreement_mat_filt_pca$rotation),
                                                kaAgreement_mat_filt_pca$rotation)

ap_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim <- autoplot(kaAgreement_mat_filt_pam,
                                                                frame = T,
                                                                frame.type = "norm",
                                                                loadings = T,
                                                                loadings.colour = "black",
                                                                loadings.label = T,
                                                                loadings.label.colour = "black",
                                                                label = F)
aps_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim <- autoplot(silhouette(kaAgreement_mat_filt_pam))
ggsave(paste0(plotDir, "ap_kaAgreement_mat_filt_pca_dim.pdf"),
       plot = ap_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)
ggsave(paste0(plotDir, "aps_kaAgreement_mat_filt_pca_dim.pdf"),
       plot = aps_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)

gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                              mapping = aes(x = PC1, y = PC2, colour = fkAgreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                              mapping = aes(x = PC1, y = PC2, colour = kaAgreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                                mapping = aes(x = PC1, y = PC2, colour = Stochasticity_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
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

gg_Superfamily_kaAgreement_mat_filt_pca_dim <- ggplot(kaAgreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = Superfamily)) +
  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = kaAgreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(kaAgreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(kaAgreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_kaAgreement_mat_filt_pca_dim_chr <- gg_kaAgreement_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_kaAgreement_mat_filt_pca_dim_chr <- gg_fkAgreement_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_kaAgreement_mat_filt_pca_dim_chr <- gg_Stochasticity_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_kaAgreement_mat_filt_pca_dim_chr <- gg_mC_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Superfamily_kaAgreement_mat_filt_pca_dim_chr <- gg_Superfamily_kaAgreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_kaAgreement_mat_filt_pca_dim_chr <- list(
                                                     gg_kaAgreement_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_fkAgreement_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_mC_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim_chr,
                                                     gg_Superfamily_kaAgreement_mat_filt_pca_dim_chr
                                                    )
gg_cow_kaAgreement_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_kaAgreement_mat_filt_pca_dim_chr,
                                                     align = "hv",
                                                     axis = "l",
                                                     nrow = length(gg_cow_list_kaAgreement_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_kaAgreement_mat_filt_pca_dim_chr.pdf"),
       plot = gg_cow_kaAgreement_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_kaAgreement_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim_loadings <- gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim_loadings <- gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim +
  geom_segment(data = kaAgreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (kaAgreement_mat_filt_pca_loadings$PC1*3.2), y = (kaAgreement_mat_filt_pca_loadings$PC2*3.2),
           label = kaAgreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim +
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

gg_Superfamily_kaAgreement_mat_filt_pca_dim_loadings <- gg_Superfamily_kaAgreement_mat_filt_pca_dim +
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
                                                          gg_kaAgreement_cluster_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_fkAgreement_cluster_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Stochasticity_cluster_kaAgreement_mat_filt_pca_dim_loadings,
                                                          gg_Superfamily_kaAgreement_mat_filt_pca_dim_loadings
                                                         )
gg_cow_kaAgreement_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings,
                                                          align = "hv",
                                                          axis = "l",
                                                          nrow = 1, ncol = length(gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_kaAgreement_mat_filt_pca_dim_loadings.pdf"),
       plot = gg_cow_kaAgreement_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_kaAgreement_mat_filt_pca_dim_loadings), limitsize = F)



# PCA of Stochasticity and mean methylation for the given context

Stochasticity_mat_filt_pca_n_dim <- 2

#Stochasticity_mat_filt_pca <- prcomp(mat_filt[ , c(1,4), drop = F], center = T, scale = T)
Stochasticity_mat_filt_pca <- prcomp(Stochasticity_mat_filt_pam$data, center = F, scale = F)
Stochasticity_mat_filt_pca_summ <- summary(Stochasticity_mat_filt_pca)
Stochasticity_mat_filt_pca_PC1_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
Stochasticity_mat_filt_pca_PC2_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
Stochasticity_mat_filt_pca_dim <- Stochasticity_mat_filt_pca$x[, seq_len(Stochasticity_mat_filt_pca_n_dim)]
colnames(Stochasticity_mat_filt_pca_dim) <- c("PC1", "PC2")
head(Stochasticity_mat_filt_pca_dim)

stopifnot(nrow(Stochasticity_mat_filt_pca_dim) == length(featDF_filt$Stochasticity_cluster))
Stochasticity_mat_filt_pca_dim <- cbind(as.data.frame(Stochasticity_mat_filt_pca_dim),
                                                      mat_filt,
                                                      chr = featDF_filt$chr,
                                                      fkAgreement_cluster = featDF_filt$fkAgreement_cluster,
                                                      kaAgreement_cluster = featDF_filt$kaAgreement_cluster,
                                                      Stochasticity_cluster = featDF_filt$Stochasticity_cluster,
                                                      Superfamily = featDF_filt$Superfamily,
                                                      type = "PCA")
head(Stochasticity_mat_filt_pca_dim)

#Stochasticity_mat_filt_pam_cl_colours <- rainbow(Stochasticity_mat_filt_pamk_n_cl) 
#Stochasticity_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = Stochasticity_mat_filt_pamk_n_cl)

Stochasticity_mat_filt_pca_loadings <- data.frame(variables = rownames(Stochasticity_mat_filt_pca$rotation),
                                                  Stochasticity_mat_filt_pca$rotation)

ap_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim <- autoplot(Stochasticity_mat_filt_pam,
                                                                    frame = T,
                                                                    frame.type = "norm",
                                                                    loadings = T,
                                                                    loadings.colour = "black",
                                                                    loadings.label = T,
                                                                    loadings.label.colour = "black",
                                                                    label = F)
aps_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim <- autoplot(silhouette(Stochasticity_mat_filt_pam))
ggsave(paste0(plotDir, "ap_Stochasticity_mat_filt_pca_dim.pdf"),
       plot = ap_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)
ggsave(paste0(plotDir, "aps_Stochasticity_mat_filt_pca_dim.pdf"),
       plot = aps_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)

gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                mapping = aes(x = PC1, y = PC2, colour = fkAgreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                mapping = aes(x = PC1, y = PC2, colour = kaAgreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                                   mapping = aes(x = PC1, y = PC2, colour = Stochasticity_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
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

gg_Superfamily_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                        mapping = aes(x = PC1, y = PC2, colour = Superfamily)) +
  geom_point(size = 0.7, alpha = 0.5) +
#  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim_chr <- gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim_chr <- gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_fkAgreement_Stochasticity_mat_filt_pca_dim_chr <- gg_fkAgreement_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_kaAgreement_Stochasticity_mat_filt_pca_dim_chr <- gg_kaAgreement_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_Stochasticity_mat_filt_pca_dim_chr <- gg_mC_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Superfamily_Stochasticity_mat_filt_pca_dim_chr <- gg_Superfamily_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_Stochasticity_mat_filt_pca_dim_chr <- list(
                                                       gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_fkAgreement_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_kaAgreement_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_mC_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Superfamily_Stochasticity_mat_filt_pca_dim_chr
                                                      )
gg_cow_Stochasticity_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_chr,
                                                       align = "hv",
                                                       axis = "l",
                                                       nrow = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_chr.pdf"),
       plot = gg_cow_Stochasticity_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_loadings <- gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim_loadings <- gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim_loadings <- gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim +
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

gg_Superfamily_Stochasticity_mat_filt_pca_dim_loadings <- gg_Superfamily_Stochasticity_mat_filt_pca_dim +
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
                                                            gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_fkAgreement_cluster_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_kaAgreement_cluster_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_Superfamily_Stochasticity_mat_filt_pca_dim_loadings
                                                           )
gg_cow_Stochasticity_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings,
                                                            align = "hv",
                                                            axis = "l",
                                                            nrow = 1, ncol = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_loadings.pdf"),
       plot = gg_cow_Stochasticity_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings), limitsize = F)


# Plot relationships and show groups
trendPlot <- function(dataFrame, mapping, paletteName, xvar, yvar, clusterlab, xlab, ylab, xaxtrans, yaxtrans, xbreaks, ybreaks, xlabels, ylabels) {
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
  labs(colour = clusterlab,
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

trendPlot2 <- function(dataFrame, mapping, paletteName, xvar, yvar, clusterlab, xlab, ylab, xaxtrans, yaxtrans, xbreaks, ybreaks, xlabels, ylabels) {
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
  scale_colour_manual(values = paletteName) +
  geom_smooth(colour = "black", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y ~ s(x, bs = "cs")) +
  labs(colour = clusterlab,
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

superfamNames <- sort(unique(featDF_filt$Superfamily))
superfamNames <- c(superfamNames[3], superfamNames[1], superfamNames[14],
                   superfamNames[10], superfamNames[7], superfamNames[13],
                   superfamNames[6], superfamNames[2], superfamNames[4],
                   superfamNames[5], superfamNames[8], superfamNames[9],
                   superfamNames[12], superfamNames[11])
superfamNames <- superfamNames[-grep("Unclassified", superfamNames)]
superfamNamesPlot <- gsub("Pogo_Tc1_Mariner", "Pogo/Tc1/Mar", superfamNames)
superfamNamesPlot <- gsub("_", " ", superfamNamesPlot)
superfamNamesPlot <- gsub("classified", ".", superfamNamesPlot)

Superfamily_colFun <- cols25(n = 25)[-c(7:16, 25)][1:length(superfamNames)]
stopifnot(length(Superfamily_colFun) == length(superfamNames))
names(Superfamily_colFun) <- superfamNames

ggTrend_mean_mC_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                   mapping = aes(x = mean_mC_all, y = fk_kappa_all, colour = fkAgreement_cluster),
                                                   paletteName = "Dark2",
                                                   xvar = mean_mC_all,
                                                   yvar = fk_kappa_all,
                                                   clusterlab = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_mean_mC_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                   mapping = aes(x = mean_mC_all, y = ka_alpha_all, colour = kaAgreement_cluster),
                                                   paletteName = "Dark2",
                                                   xvar = mean_mC_all,
                                                   yvar = ka_alpha_all,
                                                   clusterlab = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_mean_mC_all_ka_alpha_all_filt_Superfamily <- trendPlot2(dataFrame = featDF_filt,
                                                                mapping = aes(x = mean_mC_all, y = ka_alpha_all, colour = Superfamily),
                                                                paletteName = Superfamily_colFun,
                                                                xvar = mean_mC_all,
                                                                yvar = ka_alpha_all,
                                                                clusterlab = "Superfamily",
                                                                xlab = bquote(.(featName)*" mean m"*.(context)),
                                                                ylab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                                xaxtrans = log10_trans(),
                                                                yaxtrans = log10_trans(),
                                                                xbreaks = trans_breaks("log10", function(x) 10^x),
                                                                ybreaks = trans_breaks("log10", function(x) 10^x),
                                                                xlabels = trans_format("log10", math_format(10^.x)),
                                                                ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_ka_alpha_all_filt_Superfamily <- ggTrend_mean_mC_all_ka_alpha_all_filt_Superfamily +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_feature_width_ka_alpha_all_filt_Superfamily <- trendPlot2(dataFrame = featDF_filt,
                                                                  mapping = aes(x = feature_width, y = ka_alpha_all, colour = Superfamily),
                                                                  paletteName = Superfamily_colFun,
                                                                  xvar = feature_width,
                                                                  yvar = ka_alpha_all,
                                                                  clusterlab = "Superfamily",
                                                                  xlab = bquote(.(featName)*" length (bp)"),
                                                                  ylab = bquote(.(featName)*" kaAgreement (m"*.(context)*")"),
                                                                  xaxtrans = log10_trans(),
                                                                  yaxtrans = log10_trans(),
                                                                  xbreaks = trans_breaks("log10", function(x) 10^x),
                                                                  ybreaks = trans_breaks("log10", function(x) 10^x),
                                                                  xlabels = trans_format("log10", math_format(10^.x)),
                                                                  ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_feature_width_ka_alpha_all_filt_Superfamily <- ggTrend_feature_width_ka_alpha_all_filt_Superfamily +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                      mapping = aes(x = mean_mC_all, y = mean_stocha_all, colour = Stochasticity_cluster),
                                                      paletteName = "Set1",
                                                      xvar = mean_mC_all,
                                                      yvar = mean_stocha_all,
                                                      clusterlab = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
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

ggTrend_ka_alpha_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                    mapping = aes(x = ka_alpha_all, y = fk_kappa_all, colour = fkAgreement_cluster),
                                                    paletteName = "Dark2",
                                                    xvar = ka_alpha_all,
                                                    yvar = fk_kappa_all,
                                                    clusterlab = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_fk_kappa_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                    mapping = aes(x = ka_alpha_all, y = fk_kappa_all, colour = kaAgreement_cluster),
                                                    paletteName = "Dark2",
                                                    xvar = ka_alpha_all,
                                                    yvar = fk_kappa_all,
                                                    clusterlab = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_mean_stocha_all_fk_kappa_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = fkAgreement_cluster),
                                                       paletteName = "Dark2",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       clusterlab = bquote(atop("fkAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_fk_kappa_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = Stochasticity_cluster),
                                                       paletteName = "Set1",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       clusterlab = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
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

ggTrend_mean_stocha_all_ka_alpha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                       mapping = aes(x = mean_stocha_all, y = ka_alpha_all, colour = fkAgreement_cluster),
                                                       paletteName = "Dark2",
                                                       xvar = mean_stocha_all,
                                                       yvar = ka_alpha_all,
                                                       clusterlab = bquote(atop("kaAgreement &", "m"*.(context)~"cluster")),
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

ggTrend_ka_alpha_all_mean_stocha_all_filt <- trendPlot(dataFrame = featDF_filt,
                                                       mapping = aes(x = mean_stocha_all, y = ka_alpha_all, colour = Stochasticity_cluster),
                                                       paletteName = "Set1",
                                                       xvar = mean_stocha_all,
                                                       yvar = ka_alpha_all,
                                                       clusterlab = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
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
                     ggTrend_mean_mC_all_ka_alpha_all_filt_Superfamily,
                     ggTrend_feature_width_ka_alpha_all_filt_Superfamily,
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
              "_pam_clusters.pdf"),
       plot = gg_cow1,
       height = 5*length(gg_cow_list1), width = 5*length(chrName), limitsize = F)


# Extract feature clusters to enable enrichment analysis

# Filter by fk_kappa_all and mean_mC_all cluster
featDF_filt_kappa_mC_cluster1 <- featDF_filt %>%
  dplyr::filter(fkAgreement_cluster == "Cluster 1")
featDF_filt_kappa_mC_cluster2 <- featDF_filt %>%
  dplyr::filter(fkAgreement_cluster == "Cluster 2")
featDF_filt_kappa_mC_cluster3 <- featDF_filt %>%
  dplyr::filter(fkAgreement_cluster == "Cluster 3")
featDF_filt_kappa_mC_cluster4 <- featDF_filt %>%
  dplyr::filter(fkAgreement_cluster == "Cluster 4")

write.table(featDF_filt_kappa_mC_cluster1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_cluster2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_cluster3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_kappa_mC_cluster4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by ka_alpha_all and mean_mC_all cluster
featDF_filt_alpha_mC_cluster1 <- featDF_filt %>%
  dplyr::filter(kaAgreement_cluster == "Cluster 1")
featDF_filt_alpha_mC_cluster2 <- featDF_filt %>%
  dplyr::filter(kaAgreement_cluster == "Cluster 2")
featDF_filt_alpha_mC_cluster3 <- featDF_filt %>%
  dplyr::filter(kaAgreement_cluster == "Cluster 3")
featDF_filt_alpha_mC_cluster4 <- featDF_filt %>%
  dplyr::filter(kaAgreement_cluster == "Cluster 4")

write.table(featDF_filt_alpha_mC_cluster1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_cluster1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_cluster2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_cluster2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_cluster3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_cluster3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_alpha_mC_cluster4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_ka_alpha_all_mean_mC_all_cluster4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all cluster
featDF_filt_stocha_mC_cluster1 <- featDF_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 1")
featDF_filt_stocha_mC_cluster2 <- featDF_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 2")
featDF_filt_stocha_mC_cluster3 <- featDF_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 3")
featDF_filt_stocha_mC_cluster4 <- featDF_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 4")

write.table(featDF_filt_stocha_mC_cluster1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_cluster2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_cluster3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(featDF_filt_stocha_mC_cluster4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
