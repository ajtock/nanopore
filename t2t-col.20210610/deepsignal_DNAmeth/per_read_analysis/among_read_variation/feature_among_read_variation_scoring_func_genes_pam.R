#!/usr/bin/env Rscript

# Analysis:
# Group genes into clusters using partitioning around medoids (PAM),
# based on among-read agreement or site-to-site variability in DNA methylation patterns,
# and mean DNA methylation

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_pam.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions'
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
con_fk_df_all_filt <- read.table(paste0(outDir,
                                        featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                                        "_", context,
                                        "_NAmax", NAmax,
                                        "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                                        paste0(chrName, collapse = "_"), ".tsv"),
                                 header = T)
con_fk_df_all_filt$kappa_C_density <- con_fk_df_all_filt$fk_Cs_all / ( (con_fk_df_all_filt$end - con_fk_df_all_filt$start + 1) / 1e3)
con_fk_df_all_filt$stocha_C_density <- con_fk_df_all_filt$stocha_Cs_all / ( (con_fk_df_all_filt$end - con_fk_df_all_filt$start + 1) / 1e3)
con_fk_df_all_filt$parent <- sub(pattern = "\\.\\d+", replacement = "", x = con_fk_df_all_filt$name) 
con_fk_df_all_filt$parent <- sub(pattern = "_\\d+", replacement = "", x = con_fk_df_all_filt$parent) 

mat_filt <- con_fk_df_all_filt[,which(colnames(con_fk_df_all_filt) %in%
                               c("mean_mC_all", "fk_kappa_all", "mean_stocha_all")), drop = F]
colnames(mat_filt) <- c(paste0("m", context), "Agreement", "Stochasticity")

set.seed(4849345)
Agreement_mat_filt_pamk <- fpc::pamk(data = mat_filt[ , c(1,2), drop = F],
                                     krange = 3:10,
                                     criterion = "multiasw",
                                     usepam = F,
                                     scaling = T,
                                     alpha = 0.001,
                                     diss = F,
                                     critout = T,
                                     ns = 10)
Agreement_mat_filt_pamk_n_cl <- Agreement_mat_filt_pamk$nc
print(Agreement_mat_filt_pamk_n_cl)
Agreement_mat_filt_pam <- cluster::pam(x = mat_filt[ , c(1,2), drop = F],
                                       k = Agreement_mat_filt_pamk_n_cl,
                                       diss = F,
                                       metric = "manhattan",
                                       stand = T,
                                       cluster.only = F,
                                       do.swap = T,
                                       pamonce = 0)
 
con_fk_df_all_filt$Agreement_cluster <- paste0("Cluster ", Agreement_mat_filt_pam$clustering)
#con_fk_df_all_filt$Agreement_cluster <- Agreement_mat_filt_pam$clustering
#con_fk_df_all_filt$Agreement_cluster[which(con_fk_df_all_filt$Agreement_cluster == 2)] <- "Cluster 1"
#con_fk_df_all_filt$Agreement_cluster[which(con_fk_df_all_filt$Agreement_cluster == 1)] <- "Cluster 2"
#con_fk_df_all_filt$Agreement_cluster[which(con_fk_df_all_filt$Agreement_cluster == 3)] <- "Cluster 3"
#con_fk_df_all_filt$Agreement_cluster[which(con_fk_df_all_filt$Agreement_cluster == 4)] <- "Cluster 4"

set.seed(4849345)
Stochasticity_mat_filt_pamk <- fpc::pamk(data = mat_filt[ , c(1,3), drop = F],
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
Stochasticity_mat_filt_pam <- cluster::pam(x = mat_filt[ , c(1,3), drop = F],
                                           k = Stochasticity_mat_filt_pamk_n_cl,
                                           diss = F,
                                           metric = "manhattan",
                                           stand = T,
                                           cluster.only = F,
                                           do.swap = T,
                                           pamonce = 0)

con_fk_df_all_filt$Stochasticity_cluster <- paste0("Cluster ", Stochasticity_mat_filt_pam$clustering)
#con_fk_df_all_filt$Stochasticity_cluster <- Stochasticity_mat_filt_pam$clustering
#con_fk_df_all_filt$Stochasticity_cluster[which(con_fk_df_all_filt$Stochasticity_cluster == 3)] <- "Cluster 1"
#con_fk_df_all_filt$Stochasticity_cluster[which(con_fk_df_all_filt$Stochasticity_cluster == 1)] <- "Cluster 2"
#con_fk_df_all_filt$Stochasticity_cluster[which(con_fk_df_all_filt$Stochasticity_cluster == 2)] <- "Cluster 3"
#con_fk_df_all_filt$Stochasticity_cluster[which(con_fk_df_all_filt$Stochasticity_cluster == 4)] <- "Cluster 4"

# Dimension reduction:
# PCA of Agreement and mean methylation for the given context

Agreement_mat_filt_pca_n_dim <- 2

#Agreement_mat_filt_pca <- prcomp(mat_filt[ , c(1,2), drop = F], center = T, scale = T)
Agreement_mat_filt_pca <- prcomp(Agreement_mat_filt_pam$data, center = F, scale = F)
Agreement_mat_filt_pca_summ <- summary(Agreement_mat_filt_pca)
Agreement_mat_filt_pca_PC1_varexp <- round(Agreement_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
Agreement_mat_filt_pca_PC2_varexp <- round(Agreement_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
Agreement_mat_filt_pca_dim <- Agreement_mat_filt_pca$x[, seq_len(Agreement_mat_filt_pca_n_dim)]
colnames(Agreement_mat_filt_pca_dim) <- c("PC1", "PC2")
head(Agreement_mat_filt_pca_dim)

stopifnot(nrow(Agreement_mat_filt_pca_dim) == length(con_fk_df_all_filt$Agreement_cluster))
Agreement_mat_filt_pca_dim <- cbind(as.data.frame(Agreement_mat_filt_pca_dim),
                                    mat_filt,
                                    chr = con_fk_df_all_filt$chr,
                                    Agreement_cluster = con_fk_df_all_filt$Agreement_cluster,
                                    Stochasticity_cluster = con_fk_df_all_filt$Stochasticity_cluster,
                                    type = "PCA")
head(Agreement_mat_filt_pca_dim)

#Agreement_mat_filt_pam_cl_colours <- rainbow(Agreement_mat_filt_pamk_n_cl) 
#Agreement_mat_filt_pam_cl_colours <- brewer.pal(name = "Dark2", n = Agreement_mat_filt_pamk_n_cl)

Agreement_mat_filt_pca_loadings <- data.frame(variables = rownames(Agreement_mat_filt_pca$rotation),
                                              Agreement_mat_filt_pca$rotation)

ap_Agreement_cluster_Agreement_mat_filt_pca_dim <- autoplot(Agreement_mat_filt_pam,
                                                            frame = T,
                                                            frame.type = "norm",
                                                            loadings = T,
                                                            loadings.colour = "black",
                                                            loadings.label = T,
                                                            loadings.label.colour = "black",
                                                            label = F)
aps_Agreement_cluster_Agreement_mat_filt_pca_dim <- autoplot(silhouette(Agreement_mat_filt_pam))
ggsave(paste0(plotDir, "ap_Agreement_mat_filt_pca_dim.pdf"),
       plot = ap_Agreement_cluster_Agreement_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)
ggsave(paste0(plotDir, "aps_Agreement_mat_filt_pca_dim.pdf"),
       plot = aps_Agreement_cluster_Agreement_mat_filt_pca_dim,
       height = 4, width = 5, limitsize = F)

gg_Agreement_cluster_Agreement_mat_filt_pca_dim <- ggplot(Agreement_mat_filt_pca_dim,
                                                          mapping = aes(x = PC1, y = PC2, colour = Agreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Agreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Agreement &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(Agreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Agreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim <- ggplot(Agreement_mat_filt_pca_dim,
                                                              mapping = aes(x = PC1, y = PC2, colour = Stochasticity_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
#  scale_colour_manual(values = Agreement_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
       x = bquote("PC1 (" * .(Agreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Agreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Agreement_Agreement_mat_filt_pca_dim <- ggplot(Agreement_mat_filt_pca_dim,
                                                  mapping = aes(x = PC1, y = PC2, colour = Agreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "Agreement",
       x = bquote("PC1 (" * .(Agreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Agreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_Stochasticity_Agreement_mat_filt_pca_dim <- ggplot(Agreement_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = Stochasticity)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "plasma") +
  labs(colour = "Stochasticity",
       x = bquote("PC1 (" * .(Agreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Agreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

gg_mC_Agreement_mat_filt_pca_dim <- ggplot(Agreement_mat_filt_pca_dim,
                                           mapping = aes(x = PC1, y = PC2, colour = mCpG)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = bquote("m" * .(context) ~ "mean"),
       x = bquote("PC1 (" * .(Agreement_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Agreement_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_Agreement_cluster_Agreement_mat_filt_pca_dim_chr <- gg_Agreement_cluster_Agreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim_chr <- gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Agreement_Agreement_mat_filt_pca_dim_chr <- gg_Agreement_Agreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_Agreement_mat_filt_pca_dim_chr <- gg_Stochasticity_Agreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_Agreement_mat_filt_pca_dim_chr <- gg_mC_Agreement_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_Agreement_mat_filt_pca_dim_chr <- list(
                                                   gg_Agreement_Agreement_mat_filt_pca_dim_chr,
                                                   gg_Stochasticity_Agreement_mat_filt_pca_dim_chr,
                                                   gg_mC_Agreement_mat_filt_pca_dim_chr,
                                                   gg_Agreement_cluster_Agreement_mat_filt_pca_dim_chr,
                                                   gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim_chr
                                                  )
gg_cow_Agreement_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_Agreement_mat_filt_pca_dim_chr,
                                                   align = "hv",
                                                   axis = "l",
                                                   nrow = length(gg_cow_list_Agreement_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_Agreement_mat_filt_pca_dim_chr.pdf"),
       plot = gg_cow_Agreement_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_Agreement_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_Agreement_cluster_Agreement_mat_filt_pca_dim_loadings <- gg_Agreement_cluster_Agreement_mat_filt_pca_dim +
  geom_segment(data = Agreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Agreement_mat_filt_pca_loadings$PC1*3.2), y = (Agreement_mat_filt_pca_loadings$PC2*3.2),
           label = Agreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim +
  geom_segment(data = Agreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Agreement_mat_filt_pca_loadings$PC1*3.2), y = (Agreement_mat_filt_pca_loadings$PC2*3.2),
           label = Agreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Agreement_Agreement_mat_filt_pca_dim_loadings <- gg_Agreement_Agreement_mat_filt_pca_dim +
  geom_segment(data = Agreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Agreement_mat_filt_pca_loadings$PC1*3.2), y = (Agreement_mat_filt_pca_loadings$PC2*3.2),
           label = Agreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_Agreement_mat_filt_pca_dim_loadings <- gg_Stochasticity_Agreement_mat_filt_pca_dim +
  geom_segment(data = Agreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Agreement_mat_filt_pca_loadings$PC1*3.2), y = (Agreement_mat_filt_pca_loadings$PC2*3.2),
           label = Agreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_mC_Agreement_mat_filt_pca_dim_loadings <- gg_mC_Agreement_mat_filt_pca_dim +
  geom_segment(data = Agreement_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Agreement_mat_filt_pca_loadings$PC1*3.2), y = (Agreement_mat_filt_pca_loadings$PC2*3.2),
           label = Agreement_mat_filt_pca_loadings$variables, colour = "grey0")

gg_cow_list_Agreement_mat_filt_pca_dim_loadings <- list(
                                                        gg_Agreement_Agreement_mat_filt_pca_dim_loadings,
                                                        gg_Stochasticity_Agreement_mat_filt_pca_dim_loadings,
                                                        gg_mC_Agreement_mat_filt_pca_dim_loadings,
                                                        gg_Agreement_cluster_Agreement_mat_filt_pca_dim_loadings,
                                                        gg_Stochasticity_cluster_Agreement_mat_filt_pca_dim_loadings
                                                       )
gg_cow_Agreement_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_Agreement_mat_filt_pca_dim_loadings,
                                                        align = "hv",
                                                        axis = "l",
                                                        nrow = 1, ncol = length(gg_cow_list_Agreement_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_Agreement_mat_filt_pca_dim_loadings.pdf"),
       plot = gg_cow_Agreement_mat_filt_pca_dim_loadings,
       height = 4, width = 4*length(gg_cow_list_Agreement_mat_filt_pca_dim_loadings), limitsize = F)


# PCA of Stochasticity and mean methylation for the given context

Stochasticity_mat_filt_pca_n_dim <- 2

#Stochasticity_mat_filt_pca <- prcomp(mat_filt[ , c(1,3), drop = F], center = T, scale = T)
Stochasticity_mat_filt_pca <- prcomp(Stochasticity_mat_filt_pam$data, center = F, scale = F)
Stochasticity_mat_filt_pca_summ <- summary(Stochasticity_mat_filt_pca)
Stochasticity_mat_filt_pca_PC1_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,1] * 100, digits = 2)
Stochasticity_mat_filt_pca_PC2_varexp <- round(Stochasticity_mat_filt_pca_summ$importance[2,2] * 100, digits = 2)
Stochasticity_mat_filt_pca_dim <- Stochasticity_mat_filt_pca$x[, seq_len(Stochasticity_mat_filt_pca_n_dim)]
colnames(Stochasticity_mat_filt_pca_dim) <- c("PC1", "PC2")
head(Stochasticity_mat_filt_pca_dim)

stopifnot(nrow(Stochasticity_mat_filt_pca_dim) == length(con_fk_df_all_filt$Stochasticity_cluster))
Stochasticity_mat_filt_pca_dim <- cbind(as.data.frame(Stochasticity_mat_filt_pca_dim),
                                                      mat_filt,
                                                      chr = con_fk_df_all_filt$chr,
                                                      Agreement_cluster = con_fk_df_all_filt$Agreement_cluster,
                                                      Stochasticity_cluster = con_fk_df_all_filt$Stochasticity_cluster,
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

gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                              mapping = aes(x = PC1, y = PC2, colour = Agreement_cluster)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
#  scale_colour_manual(values = Stochasticity_mat_filt_pam_cl_colours) +
  labs(colour = bquote(atop("Agreement &", "m"*.(context)~"cluster")),
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

gg_Agreement_Stochasticity_mat_filt_pca_dim <- ggplot(Stochasticity_mat_filt_pca_dim,
                                                      mapping = aes(x = PC1, y = PC2, colour = Agreement)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "turbo") +
#  scale_colour_gradient(low = "red", high = "yellow") +
#  scale_colour_gradient2(low = "blue", high = "red") +
  labs(colour = "Agreement",
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
                                               mapping = aes(x = PC1, y = PC2, colour = mCpG)) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = bquote("m" * .(context) ~ "mean"),
       x = bquote("PC1 (" * .(Stochasticity_mat_filt_pca_PC1_varexp) * "%)"),
       y = bquote("PC2 (" * .(Stochasticity_mat_filt_pca_PC2_varexp) * "%)")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.key.height = unit(4, "mm"))

# Plot chr
gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim_chr <- gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Agreement_Stochasticity_mat_filt_pca_dim_chr <- gg_Agreement_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr <- gg_Stochasticity_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_mC_Stochasticity_mat_filt_pca_dim_chr <- gg_mC_Stochasticity_mat_filt_pca_dim +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list_Stochasticity_mat_filt_pca_dim_chr <- list(
                                                       gg_Stochasticity_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Agreement_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_mC_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_chr,
                                                       gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim_chr
                                                      )
gg_cow_Stochasticity_mat_filt_pca_dim_chr <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_chr,
                                                       align = "hv",
                                                       axis = "l",
                                                       nrow = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), ncol = 1)
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_chr.pdf"),
       plot = gg_cow_Stochasticity_mat_filt_pca_dim_chr,
       height = 4*length(gg_cow_list_Stochasticity_mat_filt_pca_dim_chr), width = 4*length(chrName), limitsize = F)


# Overlay loadings
gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim_loadings <- gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_loadings <- gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_Agreement_Stochasticity_mat_filt_pca_dim_loadings <- gg_Agreement_Stochasticity_mat_filt_pca_dim +
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

gg_mC_Stochasticity_mat_filt_pca_dim_loadings <- gg_mC_Stochasticity_mat_filt_pca_dim +
  geom_segment(data = Stochasticity_mat_filt_pca_loadings,
               mapping = aes(x = 0, y = 0, xend = (PC1*3), yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               colour = "grey0") +
  annotate("text", x = (Stochasticity_mat_filt_pca_loadings$PC1*3.2), y = (Stochasticity_mat_filt_pca_loadings$PC2*3.2),
           label = Stochasticity_mat_filt_pca_loadings$variables, colour = "grey0")

gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings <- list(
                                                            gg_Stochasticity_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_Agreement_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_mC_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_Stochasticity_cluster_Stochasticity_mat_filt_pca_dim_loadings,
                                                            gg_Agreement_cluster_Stochasticity_mat_filt_pca_dim_loadings
                                                           )
gg_cow_Stochasticity_mat_filt_pca_dim_loadings <- plot_grid(plotlist = gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings,
                                                            align = "hv",
                                                            axis = "l",
                                                            nrow = 1, ncol = length(gg_cow_list_Stochasticity_mat_filt_pca_dim_loadings))
ggsave(paste0(plotDir, "gg_Stochasticity_mat_filt_pca_dim_loadings.pdf"),
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

con_fk_df_all_filt_tab <- base::merge(x = con_fk_df_all_filt, y = Col_Rep1_IRratio,
                                      by.x = "parent", by.y = "parent")

print(cor.test(con_fk_df_all_filt_tab$fk_kappa_all, con_fk_df_all_filt_tab$IRratio_mean, method = "spearman"))
#-0.1034838
print(cor.test(con_fk_df_all_filt_tab$fk_kappa_all, con_fk_df_all_filt_tab$IRratio_median, method = "spearman"))
#-0.2840318

print(cor.test(con_fk_df_all_filt_tab$mean_stocha_all, con_fk_df_all_filt_tab$IRratio_mean, method = "spearman"))
#-0.1120429
print(cor.test(con_fk_df_all_filt_tab$mean_stocha_all, con_fk_df_all_filt_tab$IRratio_median, method = "spearman"))
#-0.3074874


# Plot relationships and define groups
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


ggTrend_mean_mC_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                   mapping = aes(x = mean_mC_all, y = fk_kappa_all, colour = Agreement_cluster),
                                                   paletteName = "Dark2",
                                                   xvar = mean_mC_all,
                                                   yvar = fk_kappa_all,
                                                   clusterlab = bquote(atop("Agreement &", "m"*.(context)~"cluster")),
                                                   xlab = bquote(.(featName)*" mean m"*.(context)),
                                                   ylab = bquote(.(featName)*" agreement (m"*.(context)*")"),
                                                   xaxtrans = log10_trans(),
                                                   yaxtrans = log10_trans(),
                                                   xbreaks = trans_breaks("log10", function(x) 10^x),
                                                   ybreaks = trans_breaks("log10", function(x) 10^x),
                                                   xlabels = trans_format("log10", math_format(10^.x)),
                                                   ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_mC_all_fk_kappa_all_filt <- ggTrend_mean_mC_all_fk_kappa_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_mean_mC_all_mean_stocha_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
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

ggTrend_mean_stocha_all_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = Agreement_cluster),
                                                       paletteName = "Dark2",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       clusterlab = bquote(atop("Agreement &", "m"*.(context)~"cluster")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" agreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_mean_stocha_all_fk_kappa_all_filt <- ggTrend_mean_stocha_all_fk_kappa_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_fk_kappa_all_mean_stocha_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                                       mapping = aes(x = mean_stocha_all, y = fk_kappa_all, colour = Stochasticity_cluster),
                                                       paletteName = "Set1",
                                                       xvar = mean_stocha_all,
                                                       yvar = fk_kappa_all,
                                                       clusterlab = bquote(atop("Stochasticity &", "m"*.(context)~"cluster")),
                                                       xlab = bquote(.(featName)*" mean stochasticity (m"*.(context)*")"),
                                                       ylab = bquote(.(featName)*" agreement (m"*.(context)*")"),
                                                       xaxtrans = log10_trans(),
                                                       yaxtrans = log10_trans(),
                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
                                                       xlabels = trans_format("log10", math_format(10^.x)),
                                                       ylabels = trans_format("log10", math_format(10^.x)))
ggTrend_fk_kappa_all_mean_stocha_all_filt <- ggTrend_fk_kappa_all_mean_stocha_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list1 <- list(
                     ggTrend_mean_mC_all_fk_kappa_all_filt,
                     ggTrend_mean_mC_all_mean_stocha_all_filt,
                     ggTrend_mean_stocha_all_fk_kappa_all_filt,
                     ggTrend_fk_kappa_all_mean_stocha_all_filt
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
con_fk_df_all_filt_kappa_mC_cluster1 <- con_fk_df_all_filt %>%
  dplyr::filter(Agreement_cluster == "Cluster 1")
con_fk_df_all_filt_kappa_mC_cluster2 <- con_fk_df_all_filt %>%
  dplyr::filter(Agreement_cluster == "Cluster 2")
con_fk_df_all_filt_kappa_mC_cluster3 <- con_fk_df_all_filt %>%
  dplyr::filter(Agreement_cluster == "Cluster 3")
con_fk_df_all_filt_kappa_mC_cluster4 <- con_fk_df_all_filt %>%
  dplyr::filter(Agreement_cluster == "Cluster 4")

write.table(con_fk_df_all_filt_kappa_mC_cluster1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_cluster2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_cluster3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mC_cluster4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_cluster4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all cluster
con_fk_df_all_filt_stocha_mC_cluster1 <- con_fk_df_all_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 1")
con_fk_df_all_filt_stocha_mC_cluster2 <- con_fk_df_all_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 2")
con_fk_df_all_filt_stocha_mC_cluster3 <- con_fk_df_all_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 3")
con_fk_df_all_filt_stocha_mC_cluster4 <- con_fk_df_all_filt %>%
  dplyr::filter(Stochasticity_cluster == "Cluster 4")

write.table(con_fk_df_all_filt_stocha_mC_cluster1,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_cluster2,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_cluster3,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mC_cluster4,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_cluster4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
