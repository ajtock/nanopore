#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Group reads overlapping each feature into a constant number of clusters (e.g., 2, by pam or k-means)
# 2. Score among-read variation/agreement (e.g., Fleiss' kappa) for each cluster
# 3. Identify features with high among-read agreement for each cluster,
# with low variance between clusters in among-read agreement scores

# Usage on hydrogen node7:
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_cluster_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 2 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' CEN180"

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
#library(viridis)
#library(scales)
#library(circlize)
 
#library(ggplot2)
#library(cowplot)
##library(ggcorrplot)
#library(viridis)
#library(ggthemes)
#library(tidyquant)
##library(grid)

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
                    strand = feat$strand)
} else if(featName == "gene") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/genes/", refbase, "_representative_mRNA",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand)
} else if(featName == "GYPSY") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand)
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

## Identify and remove reads whose alignment start and end coordinates are both
## contained wholly within the boundaries of the mitochondrial insertion on Chr2,
## as we cannot be sure that these reads come from the nuclear genome
#  # Get reads that overlap mito_ins_GR
#  tab_mito <- tab[tab[,1] == as.character(seqnames(mito_ins_GR)) &
#                  tab[,2] >= start(mito_ins_GR) &
#                  tab[,2] <= end(mito_ins_GR),]
#  tab_mito_reads <- unique(tab_mito[,5])
# 
#  read_within_mito_ins <- function(DSrawDF, readID, mito_ins_GR) {
#    DSrawDF_read <- DSrawDF[DSrawDF[,5] == readID,]
#    stopifnot(unique(DSrawDF_read[,1]) == as.character(seqnames(mito_ins_GR)))
#    bool <- min(DSrawDF_read[,2], na.rm = T) >= start(mito_ins_GR) &&
#            max(DSrawDF_read[,2], na.rm = T) <= end(mito_ins_GR)
#    return(bool)
#  }
# 
#  tab_mito_reads_bool <- mclapply(tab_mito_reads, function(x) {
#    read_within_mito_ins(DSrawDF = tab,
#                         readID = x,
#                         mito_ins_GR = mito_ins_GR)
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
#
#  tab_within_mito_reads <- tab_mito_reads[unlist(tab_mito_reads_bool)]
#
#  tab <- tab[!(tab[,5] %in% tab_within_mito_reads),]

# For each featName:
# 1. group overlapping reads into k clusters,
# 2. calculate a measure of among-read agreement in methylation state (e.g., Fleiss' kappa)
# for each read cluster,
# 3. calculate mean and sd cluster among-read agreement
print(outDir)
for(i in seq_along(chrName)) {
  # Get DNA methylation proportions that overlap each featName
  chr_featGR <- featGR[seqnames(featGR) == chrName[i]] 
  chr_tab <- tab[tab[,1] == chrName[i],]
  chr_tabGR <- GRanges(seqnames = chrName[i],
                       ranges = IRanges(start = chr_tab[,2],
                                        width = 1),
                       strand = chr_tab[,3],
                       read = chr_tab[,5],
                       call = chr_tab[,9])

  chr_tabGR_fwd <- chr_tabGR[strand(chr_tabGR) == "+"]
  chr_tabGR_rev <- chr_tabGR[strand(chr_tabGR) == "-"]

  ## Note: findOverlaps() approach does not work where a window does not overlap
  ##       any positions in chr_tabGR, which can occur with smaller genomeBinSize
  # Identify overlapping windows and midpoint coordinates
  fOverlaps_fwd <- findOverlaps(query = chr_featGR,
                                subject = chr_tabGR_fwd,
                                type = "any",
                                select = "all",
                                ignore.strand = T)

  fOverlaps_rev <- findOverlaps(query = chr_featGR,
                                subject = chr_tabGR_rev,
                                type = "any",
                                select = "all",
                                ignore.strand = T)

  fk_df_win_list <- mclapply(seq_along(chr_featGR), function(x) {
#  fk_df_win_list <- lapply(seq_along(chr_featGR), function(x) {
#    print(x)

    # Analyse each strand separately
    # fwd
    chr_tabGR_fwd_x <- chr_tabGR_fwd[subjectHits(fOverlaps_fwd[queryHits(fOverlaps_fwd) == x])]
    if(length(chr_tabGR_fwd_x) > 0) {
      chr_tabGR_fwd_x <- sortSeqlevels(chr_tabGR_fwd_x)
      chr_tabGR_fwd_x <- sort(chr_tabGR_fwd_x, by = ~ read + start)

      df_fwd_x <- data.frame(pos = start(chr_tabGR_fwd_x),
                             read = chr_tabGR_fwd_x$read,
                             call = chr_tabGR_fwd_x$call)

#      # tidyr::spread() is deprecated; use tidyr::pivot_wider() instead 
#      spread_fwd_x <- tidyr::spread(data = df_fwd_x,
#                                    key = read,
#                                    value = call)
##                                    sep = "_")
#      spread_fwd_x <- spread_x[ with(data = spread_x, expr = order(pos)), ]
 
      pwider_fwd_x <- as.data.frame(tidyr::pivot_wider(data = df_fwd_x,
                                                       names_from = read,
#                                                       names_prefix = "read_",
                                                       values_from = call))
      pwider_fwd_x <- pwider_fwd_x[ with(data = pwider_fwd_x, expr = order(pos)), ]
      rownames(pwider_fwd_x) <- pwider_fwd_x[,1]
      pwider_fwd_x <- pwider_fwd_x[ , -1, drop = F]

      # kappam.fleiss() uses only rows (cytosines) with complete information
      # across all columns (reads)
      # Therefore, remove columns (reads) containing > NAmax proportion NAs to
      # to retain more cytosines in the data.frame for kappa calculation

      mask_cols <- apply(pwider_fwd_x, MARGIN = 2, FUN = function(col) sum(is.na(col)) >= nrow(pwider_fwd_x) * NAmax)    
      # Report proportion of columns (reads) to be retained:
      prop_reads_retained_fwd_x <- sum(!(mask_cols)) / ncol(pwider_fwd_x)
      # Report number of columns (reads) to be retained:
      num_reads_retained_fwd_x <- sum(!(mask_cols)) 
      # Conditionally remove columns (reads) containing > NAmax proportion NAs
      if(sum(mask_cols) > 0) {
        pwider_fwd_x <- pwider_fwd_x[ , !(mask_cols), drop = F]
      }

      # Identify rows (cytosines) containing any NAs across the retained columns (reads),
      # as these will not be used by kappam.fleiss() in any case
      mask_rows <- apply(pwider_fwd_x, MARGIN = 1, FUN = function(row) sum(is.na(row)) > 0)
      # Report proportion of rows (cytosines) to be retained:
      prop_Cs_retained_fwd_x <- sum(!(mask_rows)) / nrow(pwider_fwd_x) 
      # Report number of rows (cytosines) to be retained:
      num_Cs_retained_fwd_x <- sum(!(mask_rows))
      # Conditionally remove rows (cytosines) containing any NAs
      if(sum(mask_rows) > 0) {
        pwider_fwd_x <- pwider_fwd_x[ !(mask_rows), , drop = F]
      }

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

      # Define clusters of reads within each window using
      # cluster::pam() (for predefined k) or fpc::pamk() (for dynamic k determination)
      # ("partitioning around medoids with estimation of number of clusters")
      if(nrow(pwider_fwd_x) >= min_Cs && nrow(pwider_fwd_x) <= max_Cs &&
         ncol(pwider_fwd_x) >= min_reads && ncol(pwider_fwd_x) <= max_reads) {
        set.seed(20000)
        pamk_pwider_fwd_x <- pam(x = t(pwider_fwd_x),
                                 k = k,
                                 metric = "euclidean",
                                 do.swap = T,
                                 cluster.only = T,
                                 diss = F,
                                 pamonce = 0)
  
#        htmp <- Heatmap(t(as.matrix(pwider_fwd_x)),
#                        col = c("0" = "blue", "1" = "red"),
#                        row_split = paste0("Cluster", pamk_pwider_fwd_x$clustering),
#                        show_column_dend = F, 
#                        cluster_columns = F,
#                        heatmap_legend_param = list(title = context,
#                                                    title_position = "topcenter",
#                                                    title_gp = gpar(font = 2, fontsize = 12),
#                                                    legend_direction = "horizontal",
#                                                    labels_gp = gpar(fontsize = 10)),
#                        column_title = paste0(chrName[i], ":", start(chr_featGR[x]), "-" , end(chr_featGR[x])))
#        pdf(paste0(plotDir,
#                   sampleName, "_MappedOn_", refbase, "_", context,
#                   "_read_clusters", k, "_",
#                   "_NAmax", NAmax, "_", featName, "_", x,
#                   "_", chrName[i], "_", start(chr_featGR[x]), "_", end(chr_featGR[x]), "_fwd",
#                   ".pdf"),
#            height = 10, width = 50)
#        draw(htmp,
#             heatmap_legend_side = "bottom",
#             gap = unit(c(1), "mm"))
#        dev.off()

        # Calculate Fleiss' kappa for each cluster
        fkappa_pwider_fwd_x_k_list <- lapply(1:k, function(x) {
          kappam.fleiss(pwider_fwd_x[,which(pamk_pwider_fwd_x == x)],
                        detail = F)
        })

        fkappa_pwider_fwd_x_kappa_median <- median( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$value ) )
        fkappa_pwider_fwd_x_kappa_mean <- mean( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$value ) )
        fkappa_pwider_fwd_x_kappa_sd <- sd( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$value ) )

        fkappa_pwider_fwd_x_pval_median <- median( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$p.value ) )
        fkappa_pwider_fwd_x_pval_mean <- mean( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$p.value ) )
        fkappa_pwider_fwd_x_pval_sd <- sd( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$p.value ) )

        fkappa_pwider_fwd_x_zstat_median <- median( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$statistic ) )
        fkappa_pwider_fwd_x_zstat_mean <- mean( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$statistic ) )
        fkappa_pwider_fwd_x_zstat_sd <- sd( sapply(fkappa_pwider_fwd_x_k_list, function(x) x$statistic ) )

        fkappa_pwider_fwd_x_k_reads <- sapply(1:k, function(x) fkappa_pwider_fwd_x_k_list[[x]]$raters )
        fkappa_pwider_fwd_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_reads)))
        colnames(fkappa_pwider_fwd_x_k_reads_df) <- paste0("k", 1:k, "_reads_fwd")

        fkappa_pwider_fwd_x_k_Cs <- sapply(1:k, function(x) fkappa_pwider_fwd_x_k_list[[x]]$subjects )
        fkappa_pwider_fwd_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_Cs)))
        colnames(fkappa_pwider_fwd_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_fwd")

      } else {

        fkappa_pwider_fwd_x_kappa_median <- NaN 
        fkappa_pwider_fwd_x_kappa_mean <- NaN 
        fkappa_pwider_fwd_x_kappa_sd <- NaN 

        fkappa_pwider_fwd_x_pval_median <- NaN 
        fkappa_pwider_fwd_x_pval_mean <- NaN 
        fkappa_pwider_fwd_x_pval_sd <- NaN 

        fkappa_pwider_fwd_x_zstat_median <- NaN 
        fkappa_pwider_fwd_x_zstat_mean <- NaN 
        fkappa_pwider_fwd_x_zstat_sd <- NaN 

        fkappa_pwider_fwd_x_k_reads <- sapply(1:k, function(x) NaN)
        fkappa_pwider_fwd_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_reads)))
        colnames(fkappa_pwider_fwd_x_k_reads_df) <- paste0("k", 1:k, "_reads_fwd")

        fkappa_pwider_fwd_x_k_Cs <- sapply(1:k, function(x) NaN)
        fkappa_pwider_fwd_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_Cs)))
        colnames(fkappa_pwider_fwd_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_fwd")

      }

      fk_df_fwd_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                    start = start(chr_featGR[x]),
                                    end = end(chr_featGR[x]),
                                    midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),

                                    fk_kappa_median_fwd = fkappa_pwider_fwd_x_kappa_median,
                                    fk_kappa_mean_fwd = fkappa_pwider_fwd_x_kappa_mean,
                                    fk_kappa_sd_fwd = fkappa_pwider_fwd_x_kappa_sd,

                                    fk_pval_median_fwd = fkappa_pwider_fwd_x_pval_median,
                                    fk_pval_mean_fwd = fkappa_pwider_fwd_x_pval_mean,
                                    fk_pval_sd_fwd = fkappa_pwider_fwd_x_pval_sd,

                                    fk_zstat_median_fwd = fkappa_pwider_fwd_x_zstat_median,
                                    fk_zstat_mean_fwd = fkappa_pwider_fwd_x_zstat_mean,
                                    fk_zstat_sd_fwd = fkappa_pwider_fwd_x_zstat_sd,

                                    fkappa_pwider_fwd_x_k_reads_df,
                                    fkappa_pwider_fwd_x_k_Cs_df
                                   ) 

    } else {

      fkappa_pwider_fwd_x_k_reads <- sapply(1:k, function(x) NaN)
      fkappa_pwider_fwd_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_reads)))
      colnames(fkappa_pwider_fwd_x_k_reads_df) <- paste0("k", 1:k, "_reads_fwd")

      fkappa_pwider_fwd_x_k_Cs <- sapply(1:k, function(x) NaN)
      fkappa_pwider_fwd_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_fwd_x_k_Cs)))
      colnames(fkappa_pwider_fwd_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_fwd")

      fk_df_fwd_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                    start = start(chr_featGR[x]),
                                    end = end(chr_featGR[x]),
                                    midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),

                                    fk_kappa_median_fwd = NaN,
                                    fk_kappa_mean_fwd = NaN,
                                    fk_kappa_sd_fwd = NaN,

                                    fk_pval_median_fwd = NaN,
                                    fk_pval_mean_fwd = NaN,
                                    fk_pval_sd_fwd = NaN,

                                    fk_zstat_median_fwd = NaN,
                                    fk_zstat_mean_fwd = NaN,
                                    fk_zstat_sd_fwd = NaN,

                                    fkappa_pwider_fwd_x_k_reads_df,
                                    fkappa_pwider_fwd_x_k_Cs_df
                                   ) 

    }


    # Analyse each strand separately
    # rev
    chr_tabGR_rev_x <- chr_tabGR_rev[subjectHits(fOverlaps_rev[queryHits(fOverlaps_rev) == x])]
    if(length(chr_tabGR_rev_x) > 0) {
      chr_tabGR_rev_x <- sortSeqlevels(chr_tabGR_rev_x)
      chr_tabGR_rev_x <- sort(chr_tabGR_rev_x, by = ~ read + start)

      df_rev_x <- data.frame(pos = start(chr_tabGR_rev_x),
                             read = chr_tabGR_rev_x$read,
                             call = chr_tabGR_rev_x$call)

#      # tidyr::spread() is deprecated; use tidyr::pivot_wider() instead 
#      spread_rev_x <- tidyr::spread(data = df_rev_x,
#                                    key = read,
#                                    value = call)
##                                    sep = "_")
#      spread_rev_x <- spread_x[ with(data = spread_x, expr = order(pos)), ]
 
      pwider_rev_x <- as.data.frame(tidyr::pivot_wider(data = df_rev_x,
                                                       names_from = read,
#                                                       names_prefix = "read_",
                                                       values_from = call))
      pwider_rev_x <- pwider_rev_x[ with(data = pwider_rev_x, expr = order(pos)), ]
      rownames(pwider_rev_x) <- pwider_rev_x[,1]
      pwider_rev_x <- pwider_rev_x[ , -1, drop = F]

      # kappam.fleiss() uses only rows (cytosines) with complete information
      # across all columns (reads)
      # Therefore, remove columns (reads) containing > NAmax proportion NAs to
      # to retain more cytosines in the data.frame for kappa calculation

      mask_cols <- apply(pwider_rev_x, MARGIN = 2, FUN = function(col) sum(is.na(col)) >= nrow(pwider_rev_x) * NAmax)    
      # Report proportion of columns (reads) to be retained:
      prop_reads_retained_rev_x <- sum(!(mask_cols)) / ncol(pwider_rev_x)
      # Report number of columns (reads) to be retained:
      num_reads_retained_rev_x <- sum(!(mask_cols)) 
      # Conditionally remove columns (reads) containing > NAmax proportion NAs
      if(sum(mask_cols) > 0) {
        pwider_rev_x <- pwider_rev_x[ , !(mask_cols), drop = F]
      }

      # Identify rows (cytosines) containing any NAs across the retained columns (reads),
      # as these will not be used by kappam.fleiss() in any case
      mask_rows <- apply(pwider_rev_x, MARGIN = 1, FUN = function(row) sum(is.na(row)) > 0)
      # Report proportion of rows (cytosines) to be retained:
      prop_Cs_retained_rev_x <- sum(!(mask_rows)) / nrow(pwider_rev_x) 
      # Report number of rows (cytosines) to be retained:
      num_Cs_retained_rev_x <- sum(!(mask_rows))
      # Conditionally remove rows (cytosines) containing any NAs
      if(sum(mask_rows) > 0) {
        pwider_rev_x <- pwider_rev_x[ !(mask_rows), , drop = F]
      }

      # Define clusters of reads within each window using
      # cluster::pam() (for predefined constant k) or fpc::pamk() (for dynamic k determination)
      # ("partitioning around medoids with estimation of number of clusters")
      if(nrow(pwider_rev_x) >= 10 && nrow(pwider_rev_x) <= 50 &&
         ncol(pwider_rev_x) >= 10 && ncol(pwider_rev_x) <= 50) {
        set.seed(20000)
        pamk_pwider_rev_x <- pam(x = t(pwider_rev_x),
                                 k = k,
                                 metric = "euclidean",
                                 do.swap = T,
                                 cluster.only = T,
                                 diss = F,
                                 pamonce = 0)
  
#        htmp <- Heatmap(t(as.matrix(pwider_rev_x)),
#                        col = c("0" = "blue", "1" = "red"),
#                        row_split = paste0("Cluster", pamk_pwider_rev_x$clustering),
#                        show_column_dend = F, 
#                        cluster_columns = F,
#                        heatmap_legend_param = list(title = context,
#                                                    title_position = "topcenter",
#                                                    title_gp = gpar(font = 2, fontsize = 12),
#                                                    legend_direction = "horizontal",
#                                                    labels_gp = gpar(fontsize = 10)),
#                        column_title = paste0(chrName[i], ":", start(chr_featGR[x]), "-" , end(chr_featGR[x])))
#        pdf(paste0(plotDir,
#                   sampleName, "_MappedOn_", refbase, "_", context,
#                   "_read_clusters", k, "_",
#                   "_NAmax", NAmax, "_", featName, "_", x,
#                   "_", chrName[i], "_", start(chr_featGR[x]), "_", end(chr_featGR[x]), "_rev",
#                   ".pdf"),
#            height = 10, width = 50)
#        draw(htmp,
#             heatmap_legend_side = "bottom",
#             gap = unit(c(1), "mm"))
#        dev.off()

        # Calculate Fleiss' kappa for each cluster
        fkappa_pwider_rev_x_k_list <- lapply(1:k, function(x) {
          kappam.fleiss(pwider_rev_x[,which(pamk_pwider_rev_x == x)],
                        detail = F)
        })

        fkappa_pwider_rev_x_kappa_median <- median( sapply(fkappa_pwider_rev_x_k_list, function(x) x$value ) )
        fkappa_pwider_rev_x_kappa_mean <- mean( sapply(fkappa_pwider_rev_x_k_list, function(x) x$value ) )
        fkappa_pwider_rev_x_kappa_sd <- sd( sapply(fkappa_pwider_rev_x_k_list, function(x) x$value ) )

        fkappa_pwider_rev_x_pval_median <- median( sapply(fkappa_pwider_rev_x_k_list, function(x) x$p.value ) )
        fkappa_pwider_rev_x_pval_mean <- mean( sapply(fkappa_pwider_rev_x_k_list, function(x) x$p.value ) )
        fkappa_pwider_rev_x_pval_sd <- sd( sapply(fkappa_pwider_rev_x_k_list, function(x) x$p.value ) )

        fkappa_pwider_rev_x_zstat_median <- median( sapply(fkappa_pwider_rev_x_k_list, function(x) x$statistic ) )
        fkappa_pwider_rev_x_zstat_mean <- mean( sapply(fkappa_pwider_rev_x_k_list, function(x) x$statistic ) )
        fkappa_pwider_rev_x_zstat_sd <- sd( sapply(fkappa_pwider_rev_x_k_list, function(x) x$statistic ) )

        fkappa_pwider_rev_x_k_reads <- sapply(1:k, function(x) fkappa_pwider_rev_x_k_list[[x]]$raters )
        fkappa_pwider_rev_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_reads)))
        colnames(fkappa_pwider_rev_x_k_reads_df) <- paste0("k", 1:k, "_reads_rev")

        fkappa_pwider_rev_x_k_Cs <- sapply(1:k, function(x) fkappa_pwider_rev_x_k_list[[x]]$subjects )
        fkappa_pwider_rev_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_Cs)))
        colnames(fkappa_pwider_rev_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_rev")

      } else {

        fkappa_pwider_rev_x_kappa_median <- NaN 
        fkappa_pwider_rev_x_kappa_mean <- NaN 
        fkappa_pwider_rev_x_kappa_sd <- NaN 

        fkappa_pwider_rev_x_pval_median <- NaN 
        fkappa_pwider_rev_x_pval_mean <- NaN 
        fkappa_pwider_rev_x_pval_sd <- NaN 

        fkappa_pwider_rev_x_zstat_median <- NaN 
        fkappa_pwider_rev_x_zstat_mean <- NaN 
        fkappa_pwider_rev_x_zstat_sd <- NaN 

        fkappa_pwider_rev_x_k_reads <- sapply(1:k, function(x) NaN)
        fkappa_pwider_rev_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_reads)))
        colnames(fkappa_pwider_rev_x_k_reads_df) <- paste0("k", 1:k, "_reads_rev")

        fkappa_pwider_rev_x_k_Cs <- sapply(1:k, function(x) NaN)
        fkappa_pwider_rev_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_Cs)))
        colnames(fkappa_pwider_rev_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_rev")

      }

      fk_df_rev_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                    start = start(chr_featGR[x]),
                                    end = end(chr_featGR[x]),
                                    midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),

                                    fk_kappa_median_rev = fkappa_pwider_rev_x_kappa_median,
                                    fk_kappa_mean_rev = fkappa_pwider_rev_x_kappa_mean,
                                    fk_kappa_sd_rev = fkappa_pwider_rev_x_kappa_sd,

                                    fk_pval_median_rev = fkappa_pwider_rev_x_pval_median,
                                    fk_pval_mean_rev = fkappa_pwider_rev_x_pval_mean,
                                    fk_pval_sd_rev = fkappa_pwider_rev_x_pval_sd,

                                    fk_zstat_median_rev = fkappa_pwider_rev_x_zstat_median,
                                    fk_zstat_mean_rev = fkappa_pwider_rev_x_zstat_mean,
                                    fk_zstat_sd_rev = fkappa_pwider_rev_x_zstat_sd,

                                    fkappa_pwider_rev_x_k_reads_df,
                                    fkappa_pwider_rev_x_k_Cs_df
                                   ) 

    } else {

      fkappa_pwider_rev_x_k_reads <- sapply(1:k, function(x) NaN)
      fkappa_pwider_rev_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_reads)))
      colnames(fkappa_pwider_rev_x_k_reads_df) <- paste0("k", 1:k, "_reads_rev")

      fkappa_pwider_rev_x_k_Cs <- sapply(1:k, function(x) NaN)
      fkappa_pwider_rev_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_rev_x_k_Cs)))
      colnames(fkappa_pwider_rev_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_rev")

      fk_df_rev_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                                    start = start(chr_featGR[x]),
                                    end = end(chr_featGR[x]),
                                    midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),

                                    fk_kappa_median_rev = NaN,
                                    fk_kappa_mean_rev = NaN,
                                    fk_kappa_sd_rev = NaN,

                                    fk_pval_median_rev = NaN,
                                    fk_pval_mean_rev = NaN,
                                    fk_pval_sd_rev = NaN,

                                    fk_zstat_median_rev = NaN,
                                    fk_zstat_mean_rev = NaN,
                                    fk_zstat_sd_rev = NaN,

                                    fkappa_pwider_rev_x_k_reads_df,
                                    fkappa_pwider_rev_x_k_Cs_df
                                   ) 

    }

    fkappa_pwider_all_x_k_reads <- sapply(1:k, function(x) {
      mean(c(fkappa_pwider_fwd_x_k_list[[x]]$raters, fkappa_pwider_rev_x_k_list[[x]]$raters), na.rm = T) })
    fkappa_pwider_all_x_k_reads_df <- as.data.frame(t(matrix(fkappa_pwider_all_x_k_reads)))
    colnames(fkappa_pwider_all_x_k_reads_df) <- paste0("k", 1:k, "_reads_all")

    fkappa_pwider_all_x_k_Cs <- sapply(1:k, function(x) {
      mean(c(fkappa_pwider_fwd_x_k_list[[x]]$raters, fkappa_pwider_rev_x_k_list[[x]]$subjects), na.rm = T) })
    fkappa_pwider_all_x_k_Cs_df <- as.data.frame(t(matrix(fkappa_pwider_all_x_k_Cs)))
    colnames(fkappa_pwider_all_x_k_Cs_df) <- paste0("k", 1:k, "_Cs_all")


    # Make data.frame with relevant info for genomic window
    fk_df_win_x <- data.frame(chr = seqnames(chr_featGR[x]),
                              start = start(chr_featGR[x]),
                              end = end(chr_featGR[x]),
                              midpoint = round((start(chr_featGR[x])+end(chr_featGR[x]))/2),
                              
                              fk_df_fwd_win_x[,5:ncol(fk_df_fwd_win_x)],
                              fk_df_rev_win_x[,5:ncol(fk_df_rev_win_x)],

                              fk_kappa_median_all = mean(c(fk_df_fwd_win_x$fk_kappa_median_fwd, fk_df_rev_win_x$fk_kappa_median_rev), na.rm = T),
                              fk_kappa_mean_all = mean(c(fk_df_fwd_win_x$fk_kappa_mean_fwd, fk_df_rev_win_x$fk_kappa_mean_rev), na.rm = T),
                              fk_kappa_sd_all = mean(c(fk_df_fwd_win_x$fk_kappa_sd_fwd, fk_df_rev_win_x$fk_kappa_sd_rev), na.rm = T),

                              fk_pval_median_all = mean(c(fk_df_fwd_win_x$fk_pval_median_fwd, fk_df_rev_win_x$fk_pval_median_rev), na.rm = T),
                              fk_pval_mean_all = mean(c(fk_df_fwd_win_x$fk_pval_mean_fwd, fk_df_rev_win_x$fk_pval_mean_rev), na.rm = T),
                              fk_pval_sd_all = mean(c(fk_df_fwd_win_x$fk_pval_sd_fwd, fk_df_rev_win_x$fk_pval_sd_rev), na.rm = T),

                              fk_zstat_median_all = mean(c(fk_df_fwd_win_x$fk_zstat_median_fwd, fk_df_rev_win_x$fk_zstat_median_rev), na.rm = T),
                              fk_zstat_mean_all = mean(c(fk_df_fwd_win_x$fk_zstat_mean_fwd, fk_df_rev_win_x$fk_zstat_mean_rev), na.rm = T),
                              fk_zstat_sd_all = mean(c(fk_df_fwd_win_x$fk_zstat_sd_fwd, fk_df_rev_win_x$fk_zstat_sd_rev), na.rm = T),

                              fkappa_pwider_all_x_k_reads_df,
                              fkappa_pwider_all_x_k_Cs_df
                              
                             )

    fk_df_win_x
#  })
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
                               
  fk_df <- dplyr::bind_rows(fk_df_win_list, .id = "column_label")

  fk_df <- data.frame(fk_df,
                      fk_adj_pval_median_fwd = p.adjust(fk_df$fk_pval_median_fwd, method = "BH"),
                      fk_adj_pval_mean_fwd = p.adjust(fk_df$fk_pval_median_fwd, method = "BH"),
                      fk_adj_pval_sd_fwd = p.adjust(fk_df$fk_pval_sd_fwd, method = "BH"),

                      fk_adj_pval_median_rev = p.adjust(fk_df$fk_pval_median_rev, method = "BH"),
                      fk_adj_pval_mean_rev = p.adjust(fk_df$fk_pval_median_rev, method = "BH"),
                      fk_adj_pval_sd_rev = p.adjust(fk_df$fk_pval_sd_rev, method = "BH"),

                      fk_adj_pval_median_all = p.adjust(fk_df$fk_pval_median_all, method = "BH"),
                      fk_adj_pval_mean_all = p.adjust(fk_df$fk_pval_median_all, method = "BH"),
                      fk_adj_pval_sd_all = p.adjust(fk_df$fk_pval_sd_all, method = "BH"))

  write.table(fk_df,
              file = paste0(outDir,
                            sampleName, "_MappedOn_", refbase, "_", context,
                            "_read_clusters", k, "_", featName,  genomeBinName, "_genomeStepSize", genomeStepName,
                            "_NAmax", NAmax,
                            "_min_reads", min_reads, "_max_reads", max_reads,
                            "_min_Cs", min_Cs, "_max_Cs", max_Cs,
                            "_per_read_cluster_var_df_", chrName[i], ".tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
}

#  chrPlot <- function(dataFrame, xvar, yvar, xlab, ylab, colour) {
#    xvar <- enquo(xvar)
#    yvar <- enquo(yvar)
#    ggplot(data = dataFrame,
#           mapping = aes(x = !!xvar,
#                         y = !!yvar)) +
##    geom_line(size = 1, colour = colour) +
#    geom_ma(ma_fun = SMA, n = 10, colour = colour, linetype = 1, size = 1) +
##    geom_smooth(colour = colour, fill = colour, alpha = 0.6,
##                method = "gam", formula = y ~ s(x, bs = "cs")) +
#    scale_x_continuous(
#                       labels = function(x) x/1e6) +
#    labs(x = xlab,
#         y = ylab) +
#    theme_bw() +
#    theme(
#          axis.ticks = element_line(size = 0.5, colour = "black"),
#          axis.ticks.length = unit(0.25, "cm"),
#          axis.text.x = element_text(size = 16, colour = "black"),
#          axis.text.y = element_text(size = 16, colour = "black"),
#          axis.title = element_text(size = 18, colour = "black"),
#          axis.line = element_line(size = 0.5, colour = "black"),
##          panel.grid = element_blank(),
#          panel.background = element_blank(),
##          panel.border = element_rect(size = 1.0, colour = "black"),
#          plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"))
#  }
#
#  gg_fk_kappa_all <- chrPlot(dataFrame = fk_df,
#                             xvar = midpoint,
#                             yvar = fk_kappa_all,
#                             xlab = paste0(chrName[i], " (Mb; ", genomeBinNamePlot, " windows, ", genomeStepNamePlot, " step)"),
#                             ylab = bquote("Fleiss' kappa per-read m"*.(context)),
#                             colour = "red") 
#  ggsave(paste0(plotDir,
#                sampleName, "_MappedOn_", refbase, "_", context,
#                "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#                "_NAmax", NAmax, "_Fleiss_kappa_all_", chrName[i],
#                ".pdf"),
#         plot = gg_fk_kappa_all,
#         height = 5, width = 30, limitsize = F)
#
#ggsave(paste0(plotDir,
#              "DNAmethAllContexts_",
#              gsub("_MappedOn_t2t-col.20210610_.+", "", ChIPNames)[1], "_",
#              "_avgProfiles_around",
#              "_CEN180_ranLoc_CENAthila_nonCENAthila_nonCENGypsy_in_t2t-col.20210610_",
#              paste0(chrName[i], collapse = "_"), ".pdf"),
#       plot = ggObjGA_combined,
#       height = 6.5, width = 7*5, limitsize = FALSE)
#
#
#  
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_Fleiss_kappa_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$fk_kappa_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Fleiss' kappa per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
##  par(new = T)
##  plot(x = fk_df$midpoint, y = fk_df$fk_num_Cs_all, type = "l", lwd = 0.5, col = "skyblue",
##       xaxt = "n", yaxt = "n",
##       xlab = "", ylab = "",
##       main = "")
##  p <- par('usr')
##  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.0), xpd = NA, srt = -90, col = "skyblue",
##       labels = bquote(.(context)*" sites per window"))
##  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#
## CENH3_in_bodies vs HORlengthsSum
#ggTrend5 <- ggplot(data = CEN180,
#                   mapping = aes(x = CENH3_in_bodies,
#                                 y = HORlengthsSum+1)) +
#  geom_hex(bins = 60) +
#  scale_y_continuous(trans = log10_trans(),
#                     breaks = trans_breaks("log10", function(x) 10^x),
#                     labels = trans_format("log10", math_format(10^.x))) +
#  annotation_logticks(sides = "l") +
#  scale_fill_viridis() +
#  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
#              method = "gam", formula = y ~ s(x, bs = "cs")) +
#  labs(x = "CENH3",
#       y = "Repetitiveness") +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 0.5, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 16, colour = "black"),
#        axis.text.y = element_text(size = 16, colour = "black"),
#        axis.title = element_text(size = 18, colour = "black"),
#        panel.grid = element_blank(),
#        panel.background = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 18)) +
#  ggtitle(bquote(italic(r[s]) ~ "=" ~
#                 .(round(cor.test(CEN180$CENH3_in_bodies, CEN180$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
#                         digits = 2)) *
#                 ";" ~ italic(P) ~ "=" ~
#                 .(round(min(0.5, cor.test(CEN180$CENH3_in_bodies, CEN180$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(CEN180)[1]/100) )),
#                         digits = 5)) ~
#                 "(CEN180 in" ~ .(paste0(chrName[i], collapse = ",")) * ")"))
#
#
#
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_Fleiss_kappa_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$fk_kappa_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Fleiss' kappa per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
##  par(new = T)
##  plot(x = fk_df$midpoint, y = fk_df$fk_num_Cs_all, type = "l", lwd = 0.5, col = "skyblue",
##       xaxt = "n", yaxt = "n",
##       xlab = "", ylab = "",
##       main = "")
##  p <- par('usr')
##  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.0), xpd = NA, srt = -90, col = "skyblue",
##       labels = bquote(.(context)*" sites per window"))
##  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_mean_stocha_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$mean_stocha_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Mean stoch. per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_sd_stocha_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$sd_stocha_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("SD stoch. per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_mean_mean_acf_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$mean_mean_acf_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Mean mean ACF per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_mean_min_acf_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$mean_min_acf_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Mean min. ACF per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_mean_max_acf_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$mean_max_acf_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("Mean max. ACF per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()
#
#  pdf(paste0(plotDir,
#             sampleName, "_MappedOn_", refbase, "_", context,
#             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#             "_NAmax", NAmax, "_pam_clusters_", chrName[i],
#             ".pdf"), height = 5, width = 30)
#  par(mfrow = c(1, 1))
#  par(mar = c(4.1, 4.1, 3.1, 4.1))
#  par(mgp = c(3, 1, 0))
#  plot(x = fk_df$midpoint, y = fk_df$pamk_nc_all, type = "l", lwd = 1.5, col = "red",
#       yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
#        text = bquote("PAM clusters per-read m"*.(context)))
#  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
#  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName[i], " (", genomeBinName, " window, ", genomeStepName, " step)"))
#  dev.off()


#  par(new = T)
#  plot(x = fk_df$midpoint, y = -log10(fk_df$fk_adj_pval_all+1e-10), type = "l", lwd = 1.5, col = "blue",
#       ylim = c(0,
#                pmax(-log10(0.05), max(-log10(fk_df$fk_adj_pval_all+1e-10), na.rm = T))),
#       xlab = "", ylab = "",
#       main = "")
#  p <- par('usr')
#  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3.0), xpd = NA, srt = -90, col = "blue",
#       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
#  abline(h = -log10(0.05), lty = 5, lwd = 1, col = "blue")


#  par(new = T)
#  plot(x = fk_df$midpoint, y = fk_df$fk_prop_reads_all, type = "l", col = "blue")
#  par(new = T)
#  plot(x = fk_df$midpoint, y = fk_df$fk_num_reads_all, type = "l", col = "navy")
#  par(new = T)
#  plot(x = fk_df$midpoint, y = fk_df$fk_prop_Cs_all, type = "l", col = "green")
#  par(new = T)
#  plot(x = fk_df$midpoint, y = fk_df$fk_num_Cs_all, type = "l", col = "darkgreen")
#  dev.off()

