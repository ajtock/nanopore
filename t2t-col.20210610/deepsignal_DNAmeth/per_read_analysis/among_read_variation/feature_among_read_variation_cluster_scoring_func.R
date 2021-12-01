#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Group reads overlapping each feature into a constant number of clusters (e.g., 2, by pam or k-means)
# 2. Score among-read variation/agreement (e.g., Fleiss' kappa) for each cluster
# 3. Identify features with high among-read agreement for each cluster,
# with low variance between clusters in among-read agreement scores

# Usage on hydrogen node7:
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_cluster_scoring_func.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 2 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' CEN180"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#k <- 2
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "CEN180"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
k <- as.integer(args[3])
context <- args[4]
NAmax <- as.numeric(args[5])
CPUpc <- as.numeric(args[6])
chrName <- unlist(strsplit(args[7], split = ","))
featName <- args[8]

if(context == "CpG") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
} else if(context == "CHG") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
} else if(context == "CHH") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
}

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
library(irr)
library(dplyr)
library(tidyr)
library(cluster)
library(fpc)
#library(data.table)
#library(segmentSeq)
library(ComplexHeatmap)
#library(RColorBrewer)
library(scales)
#library(circlize)
 
library(ggplot2)
library(cowplot)
#library(ggcorrplot)
library(viridis)
library(ggthemes)
library(tidyquant)
#library(grid)

outDir <- paste0("read_clusters", k, "_", featName, "/")
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

# Read in feature annotation
if(featName == "CEN180") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/CEN180/CEN180_in_", refbase,
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand",
                      "HORlengthsSum", "HORcount", "percentageIdentity")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$score)
} else if(featName == "gene") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/genes/", refbase, "_representative_mRNA",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$score)
} else if(featName == "GYPSY") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$score)
} else {
  stop(print("featName not one of CEN180, gene or GYPSY"))
}

# Load coordinates for mitochondrial insertion on Chr2, in BED format
mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
                       header = F)
colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
mito_ins_GR <- GRanges(seqnames = "Chr2",
                       ranges = IRanges(start = min(mito_ins$start)+1,
                                        end = max(mito_ins$end)),
                       strand = "*")

# Mask out featGR within mitochondrial insertion on Chr2
fOverlaps_feat_mito_ins <- findOverlaps(query = featGR,
                                        subject = mito_ins_GR,
                                        type = "any",
                                        select = "all",
                                        ignore.strand = T)
if(length(fOverlaps_feat_mito_ins) > 0) {
  featGR <- featGR[-unique(queryHits(fOverlaps_feat_mito_ins))]
}

# Read in the raw output .tsv file from Deepsignal methylation model
tab_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                    sampleName, "_MappedOn_", refbase, "_", context, "_raw_", chrName[x], ".tsv"),
             header = F)
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  tab <- dplyr::bind_rows(tab_list)
} else {
  tab <- tab_list[[1]]
}
rm(tab_list); gc()

# For each feat:
# 1. group overlapping reads into k clusters,
# 2. calculate a measure of among-read agreement in methylation state (e.g., Fleiss' kappa)
# for each read cluster,
# 3. calculate mean and sd cluster among-read agreement
  
# Get DNA methylation proportions that overlap each featName
fOverlapsStrand <- function(chr_featGR, chr_tabGR_str) {
  ## Note: findOverlaps() approach does not work where a window does not overlap
  ##       any positions in chr_tabGR, which can occur with smaller genomeBinSize
  # Identify overlapping windows and midpoint coordinates
  fOverlaps_str <- findOverlaps(query = chr_featGR,
                                subject = chr_tabGR_str,
                                type = "any",
                                select = "all",
                                ignore.strand = T)
  fOverlaps_str
}


makeDFx_strand <- function(strand, fOverlaps_str, chr_tabGR_str, chr_featGR, x) {

  chr_tabGR_str_x <- chr_tabGR_str[subjectHits(fOverlaps_str[queryHits(fOverlaps_str) == x])]

  if(length(chr_tabGR_str_x) > 0) {

    chr_tabGR_str_x <- sortSeqlevels(chr_tabGR_str_x)
    chr_tabGR_str_x <- sort(chr_tabGR_str_x, by = ~ read + start)

    df_str_x <- data.frame(pos = start(chr_tabGR_str_x),
                           read = chr_tabGR_str_x$read,
                           call = chr_tabGR_str_x$call)

    pwider_str_x <- as.data.frame(tidyr::pivot_wider(data = df_str_x,
                                                     names_from = read,
#                                                      names_prefix = "read_",
                                                     values_from = call))
    pwider_str_x <- pwider_str_x[ with(data = pwider_str_x, expr = order(pos)), ]
    rownames(pwider_str_x) <- pwider_str_x[,1]
    pwider_str_x <- pwider_str_x[ , -1, drop = F]

    # kappam.fleiss() uses only rows (cytosines) with complete information
    # across all columns (reads)
    # Therefore, remove columns (reads) containing > NAmax proportion NAs to
    # to retain more cytosines in the data.frame for kappa calculation

    mask_cols <- apply(pwider_str_x, MARGIN = 2, FUN = function(col) sum(is.na(col)) >= nrow(pwider_str_x) * NAmax)    
    # Report proportion of columns (reads) to be retained:
    prop_reads_retained_str_x <- sum(!(mask_cols)) / ncol(pwider_str_x)
    # Report number of columns (reads) to be retained:
    num_reads_retained_str_x <- sum(!(mask_cols)) 
    # Conditionally remove columns (reads) containing > NAmax proportion NAs
    if(sum(mask_cols) > 0) {
      pwider_str_x <- pwider_str_x[ , !(mask_cols), drop = F]
    }

    # Identify rows (cytosines) containing any NAs across the retained columns (reads),
    # as these will not be used by kappam.fleiss() in any case
    mask_rows <- apply(pwider_str_x, MARGIN = 1, FUN = function(row) sum(is.na(row)) > 0)
    # Report proportion of rows (cytosines) to be retained:
    prop_Cs_retained_str_x <- sum(!(mask_rows)) / nrow(pwider_str_x) 
    # Report number of rows (cytosines) to be retained:
    num_Cs_retained_str_x <- sum(!(mask_rows))
    # Conditionally remove rows (cytosines) containing any NAs
    if(sum(mask_rows) > 0) {
      pwider_str_x <- pwider_str_x[ !(mask_rows), , drop = F]
    }

    # Define clusters of reads within each window using
    # cluster::pam() (for predefined k) or fpc::pamk() (for dynamic k determination)
    # ("partitioning around medoids with estimation of number of clusters")
    if(nrow(pwider_str_x) >= min_Cs && nrow(pwider_str_x) <= max_Cs &&
       ncol(pwider_str_x) >= min_reads && ncol(pwider_str_x) <= max_reads) {

      set.seed(20000)
      pamk_pwider_str_x <- pam(x = t(pwider_str_x),
                               k = k,
                               metric = "euclidean",
                               do.swap = T,
                               cluster.only = T,
                               diss = F,
                               pamonce = 0)

#       htmp <- Heatmap(t(as.matrix(pwider_str_x)),
#                       col = c("0" = "blue", "1" = "red"),
#                       row_split = paste0("Cluster", pamk_pwider_str_x),
#                       show_column_dend = F, 
#                       cluster_columns = F,
#                       heatmap_legend_param = list(title = context,
#                                                   title_position = "topcenter",
#                                                   title_gp = gpar(font = 2, fontsize = 12),
#                                                   legend_direction = "horizontal",
#                                                   labels_gp = gpar(fontsize = 10)),
#                       column_title = paste0(chrName[i], ":", start(chr_featGR[x]), "-" , end(chr_featGR[x])))
#       pdf(paste0(plotDir,
#                  sampleName, "_MappedOn_", refbase, "_", context,
#                  "_read_clusters", k, "_",
#                  "_NAmax", NAmax, "_", featName, "_", x,
#                  "_", chrName[i], "_", start(chr_featGR[x]), "_", end(chr_featGR[x]), "_", strand,
#                  ".pdf"),
#           height = 10, width = 50)
#       draw(htmp,
#            heatmap_legend_side = "bottom",
#            gap = unit(c(1), "mm"))
#       dev.off()

      # Calculate Fleiss' kappa for each cluster
      fkappa_pwider_str_x_k_list <- lapply(1:k, function(x) {
        kappam.fleiss(pwider_str_x[,which(pamk_pwider_str_x == x)],
                      detail = F)
      })

      fkappa_pwider_str_x_kappa_median <- median( sapply(fkappa_pwider_str_x_k_list, function(x) x$value ) )
      fkappa_pwider_str_x_kappa_mean <- mean( sapply(fkappa_pwider_str_x_k_list, function(x) x$value ) )
      fkappa_pwider_str_x_kappa_sd <- sd( sapply(fkappa_pwider_str_x_k_list, function(x) x$value ) )

      fkappa_pwider_str_x_pval_median <- median( sapply(fkappa_pwider_str_x_k_list, function(x) x$p.value ) )
      fkappa_pwider_str_x_pval_mean <- mean( sapply(fkappa_pwider_str_x_k_list, function(x) x$p.value ) )
      fkappa_pwider_str_x_pval_sd <- sd( sapply(fkappa_pwider_str_x_k_list, function(x) x$p.value ) )

      fkappa_pwider_str_x_zstat_median <- median( sapply(fkappa_pwider_str_x_k_list, function(x) x$statistic ) )
      fkappa_pwider_str_x_zstat_mean <- mean( sapply(fkappa_pwider_str_x_k_list, function(x) x$statistic ) )
      fkappa_pwider_str_x_zstat_sd <- sd( sapply(fkappa_pwider_str_x_k_list, function(x) x$statistic ) )

      fkappa_pwider_str_x_k_reads <- sapply(1:k, function(x) fkappa_pwider_str_x_k_list[[x]]$raters )
      fkappa_pwider_str_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_reads)))
      colnames(fkappa_pwider_str_x_k_reads_df) <- paste0("k", 1:k, "_reads_str")

      fkappa_pwider_str_x_k_Cs <- sapply(1:k, function(x) fkappa_pwider_str_x_k_list[[x]]$subjects )
      fkappa_pwider_str_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_Cs)))
      colnames(fkappa_pwider_str_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_str")

    } else {

      fkappa_pwider_str_x_kappa_median <- NaN 
      fkappa_pwider_str_x_kappa_mean <- NaN 
      fkappa_pwider_str_x_kappa_sd <- NaN 

      fkappa_pwider_str_x_pval_median <- NaN 
      fkappa_pwider_str_x_pval_mean <- NaN 
      fkappa_pwider_str_x_pval_sd <- NaN 

      fkappa_pwider_str_x_zstat_median <- NaN 
      fkappa_pwider_str_x_zstat_mean <- NaN 
      fkappa_pwider_str_x_zstat_sd <- NaN 

      fkappa_pwider_str_x_k_reads <- sapply(1:k, function(x) NaN)
      fkappa_pwider_str_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_reads)))
      colnames(fkappa_pwider_str_x_k_reads_df) <- paste0("k", 1:k, "_reads_str")

      fkappa_pwider_str_x_k_Cs <- sapply(1:k, function(x) NaN)
      fkappa_pwider_str_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_Cs)))
      colnames(fkappa_pwider_str_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_str")

    }

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                  start = start(chr_featGR[x]),
                                  end = end(chr_featGR[x]),
                                  midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),
                                  strand = strand(chr_featGR[x]),
                                  name = chr_featGR[x]$name,
                                  score = chr_featGR[x]$score,

                                  fk_kappa_median_str = fkappa_pwider_str_x_kappa_median,
                                  fk_kappa_mean_str = fkappa_pwider_str_x_kappa_mean,
                                  fk_kappa_sd_str = fkappa_pwider_str_x_kappa_sd,

                                  fk_pval_median_str = fkappa_pwider_str_x_pval_median,
                                  fk_pval_mean_str = fkappa_pwider_str_x_pval_mean,
                                  fk_pval_sd_str = fkappa_pwider_str_x_pval_sd,

                                  fk_zstat_median_str = fkappa_pwider_str_x_zstat_median,
                                  fk_zstat_mean_str = fkappa_pwider_str_x_zstat_mean,
                                  fk_zstat_sd_str = fkappa_pwider_str_x_zstat_sd,

                                  fkappa_pwider_str_x_k_reads_df,
                                  fkappa_pwider_str_x_k_Cs_df
                                 ) 

  } else {

    fkappa_pwider_str_x_k_reads <- sapply(1:k, function(x) NaN)
    fkappa_pwider_str_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_reads)))
    colnames(fkappa_pwider_str_x_k_reads_df) <- paste0("k", 1:k, "_reads_str")

    fkappa_pwider_str_x_k_Cs <- sapply(1:k, function(x) NaN)
    fkappa_pwider_str_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_str_x_k_Cs)))
    colnames(fkappa_pwider_str_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_str")

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                  start = start(chr_featGR[x]),
                                  end = end(chr_featGR[x]),
                                  midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),
                                  strand = strand(chr_featGR[x]),
                                  name = chr_featGR[x]$name,
                                  score = chr_featGR[x]$score,

                                  fk_kappa_median_str = NaN,
                                  fk_kappa_mean_str = NaN,
                                  fk_kappa_sd_str = NaN,

                                  fk_pval_median_str = NaN,
                                  fk_pval_mean_str = NaN,
                                  fk_pval_sd_str = NaN,

                                  fk_zstat_median_str = NaN,
                                  fk_zstat_mean_str = NaN,
                                  fk_zstat_sd_str = NaN,

                                  fkappa_pwider_str_x_k_reads_df,
                                  fkappa_pwider_str_x_k_Cs_df
                                 ) 

  }

fk_df_str_win_x

}


con_fk_df_all <- data.frame()
for(chrIndex in 1:length(chrName)) {

  chr_featGR <- featGR[seqnames(featGR) == chrName[chrIndex]]
  chr_tab <- tab[tab[,1] == chrName[chrIndex],]
  chr_tabGR <- GRanges(seqnames = chrName[chrIndex],
                       ranges = IRanges(start = chr_tab[,2],
                                        width = 1),
                       strand = chr_tab[,3],
                       read = chr_tab[,5],
                       call = chr_tab[,9])
  chr_tabGR_fwd <- chr_tabGR[strand(chr_tabGR) == "+"]
  chr_tabGR_rev <- chr_tabGR[strand(chr_tabGR) == "-"]

  fOverlaps_fwd <- fOverlapsStrand(chr_tabGR_str = chr_tabGR_fwd, chr_featGR = chr_featGR)
  fOverlaps_rev <- fOverlapsStrand(chr_tabGR_str = chr_tabGR_rev, chr_featGR = chr_featGR)

  # Analyse each strand separately
  # fwd
  makeDFx_list_fwd <- mclapply(1:length(chr_featGR), function(x) {
    makeDFx_strand(strand = "fwd",
                   fOverlaps_str = fOverlaps_fwd,
                   chr_tabGR_str = chr_tabGR_fwd,
                   chr_featGR = chr_featGR,
                   x = x)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_fwd <- dplyr::bind_rows(makeDFx_list_fwd, .id = "column_label")
  
  chr_fk_df_fwd <- data.frame(chr_fk_df_fwd,
                              fk_adj_pval_median_str = p.adjust(chr_fk_df_fwd$fk_pval_median_str, method = "BH"),
                              fk_adj_pval_mean_str = p.adjust(chr_fk_df_fwd$fk_pval_median_str, method = "BH"),
                              fk_adj_pval_sd_str = p.adjust(chr_fk_df_fwd$fk_pval_sd_str, method = "BH"))
  
  makeDFx_list_rev <- mclapply(1:length(chr_featGR), function(x) {
    makeDFx_strand(strand = "rev",
                   fOverlaps_str = fOverlaps_rev,
                   chr_tabGR_str = chr_tabGR_rev,
                   chr_featGR = chr_featGR,
                   x = x)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_rev <- dplyr::bind_rows(makeDFx_list_rev, .id = "column_label")
  
  chr_fk_df_rev <- data.frame(chr_fk_df_rev,
                              fk_adj_pval_median_str = p.adjust(chr_fk_df_rev$fk_pval_median_str, method = "BH"),
                              fk_adj_pval_mean_str = p.adjust(chr_fk_df_rev$fk_pval_median_str, method = "BH"),
                              fk_adj_pval_sd_str = p.adjust(chr_fk_df_rev$fk_pval_sd_str, method = "BH"))
  
  stopifnot(identical(chr_fk_df_fwd[,1:8], chr_fk_df_rev[,1:8]))
  
  chr_fk_df_all_mean_list <- lapply(9:ncol(chr_fk_df_fwd), function(x) {
    sapply(1:nrow(chr_fk_df_fwd), function(y) {
      mean(c(chr_fk_df_fwd[y, x], chr_fk_df_rev[y, x]), na.rm = T)
    })
  })
  
  chr_fk_df_all <- data.frame(chr_fk_df_fwd[,1:8],
                              dplyr::bind_cols(chr_fk_df_all_mean_list))
  colnames(chr_fk_df_all) <- sub("_str", "_all", colnames(chr_fk_df_fwd))

  con_fk_df_all <- rbind(con_fk_df_all, chr_fk_df_all)

}

con_fk_df_all_filt <- con_fk_df_all %>%
  dplyr::filter(k1_reads_all >= 10) %>%
  dplyr::filter(k2_reads_all >= 10) %>%
  dplyr::filter(k1_Cs_all >= 2) %>%
  dplyr::filter(k2_Cs_all >= 2) %>%
  dplyr::filter(fk_kappa_sd_all <= 0.5)


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
                         digits = 5))))
}

ggTrend_score_fk_kappa_mean_all <- trendPlot(dataFrame = con_fk_df_all,
                                             xvar = score,
                                             yvar = fk_kappa_mean_all,
                                             xlab = bquote(italic(.(featName))*" divergence"),
                                             ylab = bquote(italic(.(featName))*" Fleiss' kappa (m"*.(context)*")"),
                                             xtrans = log10_trans(),
                                             ytrans = "identity",
                                             xbreaks = trans_breaks("log10", function(x) 10^x),
                                             ybreaks = waiver(),
                                             xlabels = trans_format("log10", math_format(10^.x)),
                                             ylabels = waiver())
ggTrend_score_fk_kappa_mean_all <- ggTrend_score_fk_kappa_mean_all +
  facet_grid(cols = vars(chr), scales = "free_x")

ggTrend_score_fk_kappa_mean_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
                                             xvar = score,
                                             yvar = fk_kappa_mean_all,
                                             xlab = bquote(italic(.(featName))*" divergence"),
                                             ylab = bquote(italic(.(featName))*" Fleiss' kappa (m"*.(context)*")"),
                                             xtrans = log10_trans(),
                                             ytrans = "identity",
                                             xbreaks = trans_breaks("log10", function(x) 10^x),
                                             ybreaks = waiver(),
                                             xlabels = trans_format("log10", math_format(10^.x)),
                                             ylabels = waiver())
ggTrend_score_fk_kappa_mean_all_filt <- ggTrend_score_fk_kappa_mean_all_filt +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_list1 <- list(
                     ggTrend_score_fk_kappa_mean_all,
                     ggTrend_score_fk_kappa_mean_all_filt
                    )

gg_cow1 <- plot_grid(plotlist = gg_cow_list1,
                     labels = "AUTO", label_size = 30,
                     align = "hv",
                     axis = "l",
                     nrow = length(gg_cow_list1), ncol = 1)

ggsave(paste0(plotDir,
              featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
              "_NAmax", NAmax, "_all_trendPlot_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow1,
       height = 5*length(gg_cow_list1), width = 5*length(chrName), limitsize = F)
