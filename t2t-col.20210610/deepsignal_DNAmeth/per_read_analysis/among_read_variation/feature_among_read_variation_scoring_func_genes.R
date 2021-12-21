#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Score among-read variation/agreement (e.g., Fleiss' kappa) for each feature
# 2. Examine relationships between feature among-read agreement and other metrics

# Usage on hydrogen node7:
# csmit -m 260G -c 1 "/applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_genes.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' gene"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]

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

outDir <- paste0(featName, "/")
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
                    score = feat$HORlengthsSum)
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
# 1. calculate a measure of among-read agreement in methylation state (e.g., Fleiss' kappa)

# Get reads that overlap each featName
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

# Function to calculate among-read agreement for a given feature x
makeDFx_strand <- function(fOverlaps_str, chr_tabGR_str, chr_featGR, featNum) {

  chr_tabGR_str_x <- chr_tabGR_str[subjectHits(fOverlaps_str[queryHits(fOverlaps_str) == featNum])]

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

    mean_mC_pwider_str_x <- mean(as.matrix(pwider_str_x), na.rm = T)

    # Calculate Fleiss' kappa
    if(nrow(pwider_str_x) >= min_Cs && nrow(pwider_str_x) <= max_Cs &&
       ncol(pwider_str_x) >= min_reads && ncol(pwider_str_x) <= max_reads) {

      # Calculate Fleiss' kappa
      fkappa_pwider_str_x <- kappam.fleiss(pwider_str_x, detail = F)

      # Sanity checks
      stopifnot(fkappa_pwider_str_x$raters == num_reads_retained_str_x)
      stopifnot(fkappa_pwider_str_x$subjects == num_Cs_retained_str_x)

      fkappa_pwider_str_x_kappa <- fkappa_pwider_str_x$value
      fkappa_pwider_str_x_pval <- fkappa_pwider_str_x$p.value
      fkappa_pwider_str_x_zstat <- fkappa_pwider_str_x$statistic
      fkappa_pwider_str_x_reads <- fkappa_pwider_str_x$raters
      fkappa_pwider_str_x_Cs <- fkappa_pwider_str_x$subjects

      # Calculate absolute differences between methylation statuses of neighbouring Cs within each read 
      absdiff_pwider_str_x <- abs(diff(as.matrix(pwider_str_x)))
      # Calculate the mean absolute difference for each read
      colMeans_absdiff_pwider_str_x <- colMeans(absdiff_pwider_str_x, na.rm = T)
      # Across all reads overlapping a given feature, calculate the mean of mean absolute differences
      mean_stocha_pwider_str_x <- mean(colMeans_absdiff_pwider_str_x, na.rm = T)
      # Across all reads overlapping a given feature, calculate the sd of mean absolute differences
      sd_stocha_pwider_str_x <- sd(colMeans_absdiff_pwider_str_x, na.rm = T)

      # Calculate autocorrelations between methylation statuses of neighbouring Cs within each read
      acf_pwider_str_x_list <- apply(pwider_str_x, MARGIN = 2,
                                     FUN = function(col) acf(col, lag.max = 10, plot = F, na.action = na.pass))
      mean_min_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf) != "NaN") {
          min(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1])
        } else {
          NA
        }
      }), na.rm = T)
      mean_max_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf) != "NaN") {
          max(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1])
        } else {
          NA
        }
      }), na.rm = T)
      mean_mean_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf) != "NaN") {
          mean(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1])
        } else {
          NA
        }
      }), na.rm = T)

    } else {

      fkappa_pwider_str_x_kappa <- NaN
      fkappa_pwider_str_x_pval <- NaN
      fkappa_pwider_str_x_zstat <- NaN
      fkappa_pwider_str_x_reads <- NaN
      fkappa_pwider_str_x_Cs <- NaN

      mean_stocha_pwider_str_x <- NaN
      sd_stocha_pwider_str_x <- NaN
      mean_min_acf_pwider_str_x <- NaN
      mean_max_acf_pwider_str_x <- NaN
      mean_mean_acf_pwider_str_x <- NaN

    }

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[featNum]),
                                  start = start(chr_featGR[featNum]),
                                  end = end(chr_featGR[featNum]),
                                  midpoint = round((start(chr_featGR[featNum])+end(chr_featGR[featNum]))/2),
                                  strand = strand(chr_featGR[featNum]),
                                  name = chr_featGR[featNum]$name,
                                  score = chr_featGR[featNum]$score,

                                  mean_mC_str = mean_mC_pwider_str_x,

                                  fk_kappa_str = fkappa_pwider_str_x_kappa,
                                  fk_pval_str = fkappa_pwider_str_x_pval,
                                  fk_zstat_str = fkappa_pwider_str_x_zstat,
                                  fk_reads_str = fkappa_pwider_str_x_reads,
                                  fk_Cs_str = fkappa_pwider_str_x_Cs,

                                  mean_stocha_str = mean_stocha_pwider_str_x,
                                  sd_stocha_str = sd_stocha_pwider_str_x,
                                  mean_min_acf_str = mean_min_acf_pwider_str_x,
                                  mean_max_acf_str = mean_max_acf_pwider_str_x,
                                  mean_mean_acf_str = mean_mean_acf_pwider_str_x
                                 ) 

  } else {

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[featNum]),
                                  start = start(chr_featGR[featNum]),
                                  end = end(chr_featGR[featNum]),
                                  midpoint = round((start(chr_featGR[featNum])+end(chr_featGR[featNum]))/2),
                                  strand = strand(chr_featGR[featNum]),
                                  name = chr_featGR[featNum]$name,
                                  score = chr_featGR[featNum]$score,

                                  mean_mC_str = mean_mC_pwider_str_x,

                                  fk_kappa_str = NaN,
                                  fk_pval_str = NaN,
                                  fk_zstat_str = NaN,
                                  fk_reads_str = NaN,
                                  fk_Cs_str = NaN,

                                  mean_stocha_str = NaN,
                                  sd_stocha_str = NaN,
                                  mean_min_acf_str = NaN,
                                  mean_max_acf_str = NaN,
                                  mean_mean_acf_str = NaN
                                 ) 

  }

fk_df_str_win_x

}


con_fk_df_all <- data.frame()
for(chrIndex in 1:length(chrName)) {

  print(chrName[chrIndex])

  chr_featGR <- featGR[seqnames(featGR) == chrName[chrIndex]]
  chr_tab <- tab[tab[,1] == chrName[chrIndex],]
  chr_tabGR <- GRanges(seqnames = chrName[chrIndex],
                       ranges = IRanges(start = chr_tab[,2]+1,
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
    makeDFx_strand(fOverlaps_str = fOverlaps_fwd,
                   chr_tabGR_str = chr_tabGR_fwd,
                   chr_featGR = chr_featGR,
                   featNum = x)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_fwd <- dplyr::bind_rows(makeDFx_list_fwd, .id = "column_label")
  
  chr_fk_df_fwd <- data.frame(chr_fk_df_fwd,
                              fk_adj_pval_str = p.adjust(chr_fk_df_fwd$fk_pval_str, method = "BH"))

  # rev  
  makeDFx_list_rev <- mclapply(1:length(chr_featGR), function(x) {
    makeDFx_strand(fOverlaps_str = fOverlaps_rev,
                   chr_tabGR_str = chr_tabGR_rev,
                   chr_featGR = chr_featGR,
                   featNum = x)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_rev <- dplyr::bind_rows(makeDFx_list_rev, .id = "column_label")
  
  chr_fk_df_rev <- data.frame(chr_fk_df_rev,
                              fk_adj_pval_str = p.adjust(chr_fk_df_rev$fk_pval_str, method = "BH"))
  
  stopifnot(identical(chr_fk_df_fwd[,1:8], chr_fk_df_rev[,1:8]))
 
  # Take mean of fwd and rev equivalent columns 
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

if(context == "CpG") {
  min_Cs <- 10
  max_Cs <- Inf
  min_reads <- 10
  max_reads <- Inf
} else if(context == "CHG") {
  min_Cs <- 15
  max_Cs <- Inf
  min_reads <- 15
  max_reads <- Inf
} else if(context == "CHH") {
  min_Cs <- 20
  max_Cs <- Inf
  min_reads <- 20
  max_reads <- Inf
}


con_fk_df_all_filt <- con_fk_df_all %>%
  dplyr::filter(fk_reads_all >= min_reads) %>%
  dplyr::filter(fk_Cs_all >= min_Cs)


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

#ggTrend_score_fk_kappa_all <- trendPlot(dataFrame = con_fk_df_all,
#                                        xvar = score,
#                                        yvar = fk_kappa_all,
#                                        xlab = bquote(.(featName)*" repetitiveness"),
#                                        ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
#                                        xaxtrans = log10_trans(),
#                                        yaxtrans = "identity",
#                                        xbreaks = trans_breaks("log10", function(x) 10^x),
#                                        ybreaks = waiver(),
#                                        xlabels = trans_format("log10", math_format(10^.x)),
#                                        ylabels = waiver())
#ggTrend_score_fk_kappa_all <- ggTrend_score_fk_kappa_all +
#  facet_grid(cols = vars(chr), scales = "free_x")

#ggTrend_score_fk_kappa_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
#                                             xvar = score,
#                                             yvar = fk_kappa_all,
#                                             xlab = bquote(.(featName)*" repetitiveness"),
#                                             ylab = bquote(.(featName)*" Fleiss' kappa (m"*.(context)*")"),
#                                             xaxtrans = log10_trans(),
#                                             yaxtrans = "identity",
#                                             xbreaks = trans_breaks("log10", function(x) 10^x),
#                                             ybreaks = waiver(),
#                                             xlabels = trans_format("log10", math_format(10^.x)),
#                                             ylabels = waiver())
#ggTrend_score_fk_kappa_all_filt <- ggTrend_score_fk_kappa_all_filt +
#  facet_grid(cols = vars(chr), scales = "free_x")

#ggTrend_score_fk_reads_all <- trendPlot(dataFrame = con_fk_df_all,
#                                        xvar = score,
#                                        yvar = fk_reads_all,
#                                        xlab = bquote(.(featName)*" repetitiveness"),
#                                        ylab = bquote(.(featName)*" read coverage (m"*.(context)*")"),
#                                        xaxtrans = log10_trans(),
#                                        yaxtrans = "identity",
#                                        xbreaks = trans_breaks("log10", function(x) 10^x),
#                                        ybreaks = waiver(),
#                                        xlabels = trans_format("log10", math_format(10^.x)),
#                                        ylabels = waiver())
#ggTrend_score_fk_reads_all <- ggTrend_score_fk_reads_all +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#ggTrend_score_fk_reads_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
#                                             xvar = score,
#                                             yvar = fk_reads_all,
#                                             xlab = bquote(.(featName)*" repetitiveness"),
#                                             ylab = bquote(.(featName)*" read coverage (m"*.(context)*")"),
#                                             xaxtrans = log10_trans(),
#                                             yaxtrans = "identity",
#                                             xbreaks = trans_breaks("log10", function(x) 10^x),
#                                             ybreaks = waiver(),
#                                             xlabels = trans_format("log10", math_format(10^.x)),
#                                             ylabels = waiver())
#ggTrend_score_fk_reads_all_filt <- ggTrend_score_fk_reads_all_filt +
#  facet_grid(cols = vars(chr), scales = "free_x")

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
  fk_kappa_all_high <- 0.55
  fk_kappa_all_mid  <- 0.35
  fk_kappa_all_low  <- 0.04
  mean_stocha_all_high <- 0.28
  mean_stocha_all_mid  <- 0.17
  mean_stocha_all_low  <- 0.08
  mean_min_acf_all_high <- -0.05
  mean_min_acf_all_mid  <- -0.10
  mean_min_acf_all_low  <- -0.15
  mean_mC_all_high  <- 0.50
  mean_mC_all_mid   <- 0.15
  mean_mC_all_low   <- 0.05
} else if(context == "CHH") {
  fk_kappa_all_high <- 0.55
  fk_kappa_all_mid  <- 0.35
  fk_kappa_all_low  <- 0.04
  mean_stocha_all_high <- 0.28
  mean_stocha_all_mid  <- 0.17
  mean_stocha_all_low  <- 0.08
  mean_min_acf_all_high <- -0.05
  mean_min_acf_all_mid  <- -0.10
  mean_min_acf_all_low  <- -0.15
  mean_mC_all_high  <- 0.18
  mean_mC_all_mid   <- 0.06
  mean_mC_all_low   <- 0.02
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

#ggTrend_mean_mC_all_mean_min_acf_all <- trendPlot(dataFrame = con_fk_df_all,
#                                                  mapping = aes(x = mean_mC_all, y = mean_min_acf_all),
#                                                  xvar = mean_mC_all,
#                                                  yvar = mean_min_acf_all,
#                                                  xlab = bquote(.(featName)*" mean m"*.(context)),
#                                                  ylab = bquote(.(featName)*" mean min. ACF (m"*.(context)*")"),
#                                                  xaxtrans = log10_trans(),
#                                                  yaxtrans = log10_trans(),
#                                                  xbreaks = trans_breaks("log10", function(x) 10^x),
#                                                  ybreaks = trans_breaks("log10", function(x) 10^x),
#                                                  xlabels = trans_format("log10", math_format(10^.x)),
#                                                  ylabels = trans_format("log10", math_format(10^.x)))
#ggTrend_mean_mC_all_mean_min_acf_all <- ggTrend_mean_mC_all_mean_min_acf_all +
##  geom_hline(yintercept = mean_min_acf_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
##  geom_hline(yintercept = mean_min_acf_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
##  geom_hline(yintercept = mean_min_acf_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
#  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
#  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
#  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#ggTrend_mean_mC_all_mean_min_acf_all_filt <- trendPlot(dataFrame = con_fk_df_all_filt,
#                                                       mapping = aes(x = mean_mC_all, y = mean_min_acf_all),
#                                                       xvar = mean_mC_all,
#                                                       yvar = mean_min_acf_all,
#                                                       xlab = bquote(.(featName)*" mean m"*.(context)),
#                                                       ylab = bquote(.(featName)*" mean min. ACF (m"*.(context)*")"),
#                                                       xaxtrans = log10_trans(),
#                                                       yaxtrans = log10_trans(),
#                                                       xbreaks = trans_breaks("log10", function(x) 10^x),
#                                                       ybreaks = trans_breaks("log10", function(x) 10^x),
#                                                       xlabels = trans_format("log10", math_format(10^.x)),
#                                                       ylabels = trans_format("log10", math_format(10^.x)))
#ggTrend_mean_mC_all_mean_min_acf_all_filt <- ggTrend_mean_mC_all_mean_min_acf_all_filt +
##  geom_hline(yintercept = mean_min_acf_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
##  geom_hline(yintercept = mean_min_acf_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
##  geom_hline(yintercept = mean_min_acf_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
#  geom_vline(xintercept = mean_mC_all_high, linetype = "dashed", size = 1, colour = "darkorange1") +
#  geom_vline(xintercept = mean_mC_all_mid, linetype = "dashed", size = 1, colour = "magenta1") +
#  geom_vline(xintercept = mean_mC_all_low, linetype = "dashed", size = 1, colour = "dodgerblue1") +
#  facet_grid(cols = vars(chr), scales = "free_x")

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

gg_cow_list1 <- list(
                     ggTrend_mean_mC_all_fk_kappa_all,
                     ggTrend_mean_mC_all_fk_kappa_all_filt,
                     ggTrend_mean_mC_all_mean_stocha_all,
                     ggTrend_mean_mC_all_mean_stocha_all_filt,
#                     ggTrend_mean_mC_all_mean_min_acf_all,
#                     ggTrend_mean_mC_all_mean_min_acf_all_filt,
                     ggTrend_fk_reads_all_fk_kappa_all,
                     ggTrend_fk_reads_all_fk_kappa_all_filt,
                     ggTrend_fk_reads_all_mean_stocha_all,
                     ggTrend_fk_reads_all_mean_stocha_all_filt,
                     ggTrend_fk_Cs_all_fk_kappa_all,
                     ggTrend_fk_Cs_all_fk_kappa_all_filt,
                     ggTrend_fk_Cs_all_mean_stocha_all,
                     ggTrend_fk_Cs_all_mean_stocha_all_filt
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


# Extract feature groups (based on trend plots) to enable enrichment analysis

# Filter by fk_kappa_all and mean_mC_all
con_fk_df_all_filt_kappa_low_mC_low_group1 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_low)

con_fk_df_all_filt_kappa_mid_mC_low_group2 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all >  fk_kappa_all_low) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_low)

con_fk_df_all_filt_kappa_mid_mC_mid_group3 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all <= fk_kappa_all_mid) %>%
  dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_mid)

con_fk_df_all_filt_kappa_high_mC_mid_group4 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all >  fk_kappa_all_mid) %>%
  dplyr::filter(mean_mC_all  >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_mid)

con_fk_df_all_filt_kappa_high_mC_high_group5 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all <=  fk_kappa_all_high) %>%
  dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_high)

con_fk_df_all_filt_kappa_vhigh_mC_high_group6 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all >  fk_kappa_all_high) %>%
  dplyr::filter(mean_mC_all  >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all  <= mean_mC_all_high)

con_fk_df_all_filt_kappa_low_mC_vhigh_group7 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all <= fk_kappa_all_low) %>%
  dplyr::filter(mean_mC_all  > mean_mC_all_high)

con_fk_df_all_filt_kappa_mid_mC_vhigh_group8 <- con_fk_df_all_filt %>%
  dplyr::filter(fk_kappa_all > fk_kappa_all_low) %>%
  dplyr::filter(mean_mC_all  >  mean_mC_all_high)

write.table(con_fk_df_all_filt_kappa_low_mC_low_group1,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mid_mC_low_group2,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mid_mC_mid_group3,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_high_mC_mid_group4,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_high_mC_high_group5,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_vhigh_mC_high_group6,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_low_mC_vhigh_group7,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mid_mC_vhigh_group8,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all
con_fk_df_all_filt_stocha_low_mC_low_group1 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_low)

con_fk_df_all_filt_stocha_mid_mC_low_group2 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_low)

con_fk_df_all_filt_stocha_mid_mC_mid_group3 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

con_fk_df_all_filt_stocha_high_mC_mid_group4 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

con_fk_df_all_filt_stocha_high_mC_high_group5 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <=  mean_stocha_all_high) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_high)

con_fk_df_all_filt_stocha_vhigh_mC_high_group6 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_high) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_high)

con_fk_df_all_filt_stocha_mid_mC_vhigh_group7 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_high)

con_fk_df_all_filt_stocha_high_mC_vhigh_group8 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all > mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_high)

write.table(con_fk_df_all_filt_stocha_low_mC_low_group1,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_low_group2,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_mid_group3,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_mid_group4,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_high_group5,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_vhigh_mC_high_group6,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_vhigh_group7,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_vhigh_group8,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


