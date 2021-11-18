#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Among-read variation/agreement (Fleiss' kappa) [this script]
# 2. Variance of read stochasticity (variance of per-read, per-window proportion of inter-C intervals that represent a methylation state change [variance across all reads that overlap window examined])
# 3. Per-read autocorrelation, but could be hard to derive a general measure across reads overlapping a given region (maybe look at variance of pairwise read correlation coefficients between autocorrelation values, either standardised to e.g. 100 autocorrelation values per read, or use only reads with info across the same Cs) 

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CHG 0.50 1.00 Chr1"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CpG 0.50 1.00 Chr1"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1", split = ","))

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
library(GenomicRanges)
library(irr)
library(dplyr)
library(tidyr)
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

# Identify and remove reads whose alignment start and end coordinates are both
# contained wholly within the boundaries of the mitochondrial insertion on Chr2,
# as we cannot be sure that these reads come from the nuclear genome
if(length(chrName) == 1 && chrName == "Chr2") {

  tab <- tab[tab[,1] == chrName,]
 
  # Get reads that overlap mito_ins_GR
  tab_mito <- tab[tab[,1] == as.character(seqnames(mito_ins_GR)) &
                  tab[,2] >= start(mito_ins_GR) &
                  tab[,2] <= end(mito_ins_GR),]
  tab_mito_reads <- unique(tab_mito[,5])
 
  read_within_mito_ins <- function(DSrawDF, readID, mito_ins_GR) {
    DSrawDF_read <- DSrawDF[DSrawDF[,5] == readID,]
    stopifnot(unique(DSrawDF_read[,1]) == as.character(seqnames(mito_ins_GR)))
    bool <- min(DSrawDF_read[,2], na.rm = T) >= start(mito_ins_GR) &&
            max(DSrawDF_read[,2], na.rm = T) <= end(mito_ins_GR)
    return(bool)
  }
 
  tab_mito_reads_bool <- mclapply(tab_mito_reads, function(x) {
    read_within_mito_ins(DSrawDF = tab,
                         readID = x,
                         mito_ins_GR = mito_ins_GR)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

  tab_within_mito_reads <- tab_mito_reads[unlist(tab_mito_reads_bool)]

  tab <- tab[!(tab[,5] %in% tab_within_mito_reads),]

}

# For each genomeBinSize-bp window with a step of genomeStepSize-bp,
# calculate Fleiss' kappa statistic as a measure of among-read agreement
# in methylation state
print(genomeBinName)
print(genomeStepName)
for(i in seq_along(chrName)) {
  # Define sliding windows of width genomeBinSize bp,
  # with a step of genomeStepSize vp
  ## Note: the active code creates windows of genomeBinSize bp only,
  ## whereas the commented-out code creates windows decreasing from genomeBinSize bp to genomeStepSize bp
  ## at the right-hand end of each chromosome ( from chrLens[i]-genomeBinSize to chrLens[i] ),
  winStarts <- seq(from = 1,
#                   to = chrLens[i],
                   to = chrLens[i]-genomeBinSize,
                   by = genomeStepSize)
#  stopifnot(winStarts[length(winStarts)] == chrLens[i])
  if(chrLens[i] - winStarts[length(winStarts)] >= genomeBinSize) {
    winStarts <- c(winStarts,
                   winStarts[length(winStarts)]+genomeStepSize)
  }
  winEnds <- seq(from = winStarts[1]+genomeBinSize-1,
                 to = chrLens[i],
                 by = genomeStepSize)
  winEnds <- c(winEnds,
               rep(chrLens[i], times = length(winStarts)-length(winEnds)))
  stopifnot(winEnds[length(winEnds)] == chrLens[i])
  stopifnot(length(winStarts) == length(winEnds))

  winGR <- GRanges(seqnames = chrName[i],
                   ranges = IRanges(start = winStarts,
                                    end = winEnds),
                   strand = "*")
  print(winGR)

#  # Define adjacent windows
#  winSeq <- seq(from = 1, to = chrLens[i], by = genomeBinSize)
#  winIR <- IRanges(start = winSeq,
#                   width = genomeBinSize)
#  winIR <- winIR[-length(winIR)]
#  winIR <- append(winIR,
#                  IRanges(start = winSeq[length(winSeq)],
#                          end = chrLens[i]))
#  winGR <- GRanges(seqnames = chrName[i],
#                   ranges = winIR,
#                   strand = "*")
#  print(winGR)

  # Define per-read-window midpoint coordinates as GRanges objects
  # and get corresponding DNA methylation proportions that overlap genomic windows
  chr_tab <- tab[tab[,1] == chrName[i],]
  chr_tab_GR <- GRanges(seqnames = chrName[i],
                        ranges = IRanges(start = chr_tab[,2],
                                         width = 1),
                        strand = chr_tab[,3],
                        read = chr_tab[,5],
                        call = chr_tab[,9])

  chr_tab_GR_fwd <- chr_tab_GR[strand(chr_tab_GR) == "+"]
#  chr_tab_GR_fwd <- sort(chr_tab_GR_fwd, by = ~ read + start)
  chr_tab_GR_rev <- chr_tab_GR[strand(chr_tab_GR) == "-"]
#  chr_tab_GR_rev <- sort(chr_tab_GR_rev, by = ~ read + start)

  ## Note: findOverlaps() approach does not work where a window does not overlap
  ##       any positions in chr_tab_GR, which can occur with smaller genomeBinSize
  # Identify overlapping windows and midpoint coordinates
  fOverlaps_fwd <- findOverlaps(query = winGR,
                                subject = chr_tab_GR_fwd,
                                type = "any",
                                select = "all",
                                ignore.strand = T)

  fOverlaps_rev <- findOverlaps(query = winGR,
                                subject = chr_tab_GR_rev,
                                type = "any",
                                select = "all",
                                ignore.strand = T)

  fk_df_win_list <- mclapply(seq_along(winGR), function(x) {
#  fk_df_win_list <- lapply(seq_along(winGR), function(x) {
#    print(x)

    # Analyse each strand separately
    # fwd
    chr_tab_GR_fwd_x <- chr_tab_GR_fwd[subjectHits(fOverlaps_fwd[queryHits(fOverlaps_fwd) == x])]
    if(length(chr_tab_GR_fwd_x) > 0) {
      chr_tab_GR_fwd_x <- sortSeqlevels(chr_tab_GR_fwd_x)
      chr_tab_GR_fwd_x <- sort(chr_tab_GR_fwd_x, by = ~ read + start)

      df_fwd_x <- data.frame(pos = start(chr_tab_GR_fwd_x),
                             read = chr_tab_GR_fwd_x$read,
                             call = chr_tab_GR_fwd_x$call)

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
      # across all columns (reads) considered
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
#      # Conditionally remove rows (cytosines) containing any NAs
#      if(sum(mask_rows) > 0) {
#        pwider_fwd_x <- pwider_fwd_x[ !(mask_rows), , drop = F]
#      }

      # Calculate Fleiss' kappa
      fkappa_pwider_fwd_x <- kappam.fleiss(pwider_fwd_x, detail = F)

      # Sanity checks
      stopifnot(fkappa_pwider_fwd_x$raters == num_reads_retained_fwd_x)
      stopifnot(fkappa_pwider_fwd_x$subjects == num_Cs_retained_fwd_x)


      # Calculate within-read site-to-site stochasticity in methylation status
      # Get absolute differences in methylation status between Cs for each read
      if(nrow(pwider_fwd_x) > 1) {
        absdiff_pwider_fwd_x <- abs(diff(as.matrix(pwider_fwd_x)))
        # Calculate the mean absolute difference for each read
        colMeans_absdiff_pwider_fwd_x <- colMeans(absdiff_pwider_fwd_x, na.rm = T)
        # Across all reads overlapping a given window, calculate the mean of mean absolute differences
        mean_stocha_pwider_fwd_x <- mean(colMeans_absdiff_pwider_fwd_x, na.rm = T)
        # Across all reads overlapping a given window, calculate the sd of mean absolute differences
        sd_stocha_pwider_fwd_x <- sd(colMeans_absdiff_pwider_fwd_x, na.rm = T)

        # Calculate autocorrelations between methylation statuses of neighbouring Cs within each read
        acf_pwider_fwd_x_list <- apply(pwider_fwd_x, MARGIN = 2,
                                       FUN = function(col) acf(col, lag.max = 10, plot = F, na.action = na.pass))
        mean_min_acf_pwider_fwd_x <- mean(sapply(seq_along(acf_pwider_fwd_x_list), function(col) {
          min(as.vector(acf_pwider_fwd_x_list[[col]]$acf)[-1])
        }))
        mean_max_acf_pwider_fwd_x <- mean(sapply(seq_along(acf_pwider_fwd_x_list), function(col) {
          max(as.vector(acf_pwider_fwd_x_list[[col]]$acf)[-1])
        }))
        mean_mean_acf_pwider_fwd_x <- mean(sapply(seq_along(acf_pwider_fwd_x_list), function(col) {
          mean(as.vector(acf_pwider_fwd_x_list[[col]]$acf)[-1])
        }))
      } else {
        mean_stocha_pwider_fwd_x <- NaN
        sd_stocha_pwider_fwd_x <- NaN
        mean_min_acf_pwider_fwd_x <- NaN
        mean_max_acf_pwider_fwd_x <- NaN
        mean_mean_acf_pwider_fwd_x <- NaN
      }  

#      # Define clusters of reads within each window using fpc::pamk()
#      # ("partitioning around medoids with estimation of number of clusters")
#      set.seed(20000)
#      pamk_pwider_fwd_x <- pamk(t(pwider_fwd_x),
##                                krange = 1:(nrow(t(pwider_fwd_x))-1),
#                                krange = 1:(round(nrow(t(pwider_fwd_x))/2)),
#                                criterion = "asw",
#                                usepam = T,
#                                scaling = F,
#                                alpha = 0.001,
##                                ns = 10,
##                                seed = 100000,
#                                diss = F)

#      htmp <- Heatmap(t(as.matrix(pwider_fwd_x)),
#                      col = c("0" = "blue", "1" = "red"),
#                      row_split = paste0("Cluster", pamk_pwider_fwd_x$pamobject$clustering),
#                      show_column_dend = F, 
#                      cluster_columns = F,
#                      heatmap_legend_param = list(title = context,
#                                                  title_position = "topcenter",
#                                                  title_gp = gpar(font = 2, fontsize = 12),
#                                                  legend_direction = "horizontal",
#                                                  labels_gp = gpar(fontsize = 10)),
#                      column_title = paste0(chrName[i], ":", start(winGR[x]), "-" , end(winGR[x])))
#      pdf(paste0(plotDir,
#                 sampleName, "_MappedOn_", refbase, "_", context,
#                 "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#                 "_NAmax", NAmax, "_", chrName[i], "_", start(winGR[x]), "_", end(winGR[x]), "_fwd",
#                 ".pdf"),
#          height = 10, width = 50)
#      draw(htmp,
#           heatmap_legend_side = "bottom",
#           gap = unit(c(1), "mm"))
#      dev.off()

      fk_df_fwd_win_x <- data.frame(chr = seqnames(winGR[x]),
                                    start = start(winGR[x]),
                                    end = end(winGR[x]),
                                    midpoint = round((start(winGR[x])+end(winGR[x]))/2),

                                    fk_kappa_fwd = fkappa_pwider_fwd_x$value,
                                    fk_pval_fwd = fkappa_pwider_fwd_x$p.value,
                                    fk_zstat_fwd = fkappa_pwider_fwd_x$statistic,
                                    fk_num_reads_fwd = fkappa_pwider_fwd_x$raters,
                                    fk_num_Cs_fwd = fkappa_pwider_fwd_x$subjects,
                                    fk_prop_reads_fwd = prop_reads_retained_fwd_x,
                                    fk_prop_Cs_fwd = prop_Cs_retained_fwd_x,
 
#                                    pamk_nc_fwd = pamk_pwider_fwd_x$nc,

                                    mean_stocha_fwd = mean_stocha_pwider_fwd_x,
                                    sd_stocha_fwd = sd_stocha_pwider_fwd_x,
                                    mean_min_acf_fwd = mean_min_acf_pwider_fwd_x,
                                    mean_max_acf_fwd = mean_max_acf_pwider_fwd_x,
                                    mean_mean_acf_fwd = mean_mean_acf_pwider_fwd_x)

    } else {

      fk_df_fwd_win_x <- data.frame(chr = seqnames(winGR[x]),
                                    start = start(winGR[x]),
                                    end = end(winGR[x]),
                                    midpoint = round((start(winGR[x])+end(winGR[x]))/2),

                                    fk_kappa_fwd = NaN,
                                    fk_pval_fwd = NaN,
                                    fk_zstat_fwd = NaN,
                                    fk_num_reads_fwd = NaN,
                                    fk_num_Cs_fwd = NaN,
                                    fk_prop_reads_fwd = NaN,
                                    fk_prop_Cs_fwd = NaN,
 
#                                    pamk_nc_fwd = NaN,

                                    mean_stocha_fwd = NaN,
                                    sd_stocha_fwd = NaN,
                                    mean_min_acf_fwd = NaN,
                                    mean_max_acf_fwd = NaN,
                                    mean_mean_acf_fwd = NaN)
 
    }


    # Analyse each strand separately
    # rev
    chr_tab_GR_rev_x <- chr_tab_GR_rev[subjectHits(fOverlaps_rev[queryHits(fOverlaps_rev) == x])]
    if(length(chr_tab_GR_rev_x) > 0) {
      chr_tab_GR_rev_x <- sortSeqlevels(chr_tab_GR_rev_x)
      chr_tab_GR_rev_x <- sort(chr_tab_GR_rev_x, by = ~ read + start)

      df_rev_x <- data.frame(pos = start(chr_tab_GR_rev_x),
                             read = chr_tab_GR_rev_x$read,
                             call = chr_tab_GR_rev_x$call)

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
      # across all columns (reads) considered
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
#      # Conditionally remove rows (cytosines) containing any NAs
#      if(sum(mask_rows) > 0) {
#        pwider_rev_x <- pwider_rev_x[ !(mask_rows), , drop = F]
#      }

      # Calculate Fleiss' kappa
      fkappa_pwider_rev_x <- kappam.fleiss(pwider_rev_x, detail = F)

      # Sanity checks
      stopifnot(fkappa_pwider_rev_x$raters == num_reads_retained_rev_x)
      stopifnot(fkappa_pwider_rev_x$subjects == num_Cs_retained_rev_x)


      if(nrow(pwider_rev_x) > 1) {
        absdiff_pwider_rev_x <- abs(diff(as.matrix(pwider_rev_x)))
        # Calculate the mean absolute difference for each read
        colMeans_absdiff_pwider_rev_x <- colMeans(absdiff_pwider_rev_x, na.rm = T)
        # Across all reads overlapping a given window, calculate the mean of mean absolute differences
        mean_stocha_pwider_rev_x <- mean(colMeans_absdiff_pwider_rev_x, na.rm = T)
        # Across all reads overlapping a given window, calculate the sd of mean absolute differences
        sd_stocha_pwider_rev_x <- sd(colMeans_absdiff_pwider_rev_x, na.rm = T)

        # Calculate autocorrelations between methylation statuses of neighbouring Cs within each read
        acf_pwider_rev_x_list <- apply(pwider_rev_x, MARGIN = 2,
                                       FUN = function(col) acf(col, lag.max = 10, plot = F, na.action = na.pass))
        mean_min_acf_pwider_rev_x <- mean(sapply(seq_along(acf_pwider_rev_x_list), function(col) {
          min(as.vector(acf_pwider_rev_x_list[[col]]$acf)[-1])
        }))
        mean_max_acf_pwider_rev_x <- mean(sapply(seq_along(acf_pwider_rev_x_list), function(col) {
          max(as.vector(acf_pwider_rev_x_list[[col]]$acf)[-1])
        }))
        mean_mean_acf_pwider_rev_x <- mean(sapply(seq_along(acf_pwider_rev_x_list), function(col) {
          mean(as.vector(acf_pwider_rev_x_list[[col]]$acf)[-1])
        }))
      } else {
        mean_stocha_pwider_rev_x <- NaN
        sd_stocha_pwider_rev_x <- NaN
        mean_min_acf_pwider_rev_x <- NaN
        mean_max_acf_pwider_rev_x <- NaN
        mean_mean_acf_pwider_rev_x <- NaN
      }  


#      # Define clusters of reads within each window using fpc::pamk()
#      # ("partitioning around medoids with estimation of number of clusters")
#      set.seed(20000)
#      pamk_pwider_rev_x <- pamk(t(pwider_rev_x),
##                                krange = 1:(nrow(t(pwider_rev_x))-1),
#                                krange = 1:(round(nrow(t(pwider_rev_x))/2)),
#                                criterion = "asw",
#                                usepam = T,
#                                scaling = F,
#                                alpha = 0.001,
##                                ns = 10,
##                                seed = 100000,
#                                diss = F)

#      htmp <- Heatmap(t(as.matrix(pwider_rev_x)),
#                      col = c("0" = "blue", "1" = "red"),
#                      row_split = paste0("Cluster", pamk_pwider_rev_x$pamobject$clustering),
#                      show_column_dend = F, 
#                      cluster_columns = F,
#                      heatmap_legend_param = list(title = context,
#                                                  title_position = "topcenter",
#                                                  title_gp = gpar(font = 2, fontsize = 12),
#                                                  legend_direction = "horizontal",
#                                                  labels_gp = gpar(fontsize = 10)),
#                      column_title = paste0(chrName[i], ":", start(winGR[x]), "-" , end(winGR[x])))
#      pdf(paste0(plotDir,
#                 sampleName, "_MappedOn_", refbase, "_", context,
#                 "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#                 "_NAmax", NAmax, "_", chrName[i], "_", start(winGR[x]), "_", end(winGR[x]), "_rev",
#                 ".pdf"),
#          height = 10, width = 50)
#      draw(htmp,
#           heatmap_legend_side = "bottom",
#           gap = unit(c(1), "mm"))
#      dev.off()

      fk_df_rev_win_x <- data.frame(chr = seqnames(winGR[x]),
                                    start = start(winGR[x]),
                                    end = end(winGR[x]),
                                    midpoint = round((start(winGR[x])+end(winGR[x]))/2),

                                    fk_kappa_rev = fkappa_pwider_rev_x$value,
                                    fk_pval_rev = fkappa_pwider_rev_x$p.value,
                                    fk_zstat_rev = fkappa_pwider_rev_x$statistic,
                                    fk_num_reads_rev = fkappa_pwider_rev_x$raters,
                                    fk_num_Cs_rev = fkappa_pwider_rev_x$subjects,
                                    fk_prop_reads_rev = prop_reads_retained_rev_x,
                                    fk_prop_Cs_rev = prop_Cs_retained_rev_x,
 
#                                    pamk_nc_rev = pamk_pwider_rev_x$nc,

                                    mean_stocha_rev = mean_stocha_pwider_rev_x,
                                    sd_stocha_rev = sd_stocha_pwider_rev_x,
                                    mean_min_acf_rev = mean_min_acf_pwider_rev_x,
                                    mean_max_acf_rev = mean_max_acf_pwider_rev_x,
                                    mean_mean_acf_rev = mean_mean_acf_pwider_rev_x)

    } else {

      fk_df_rev_win_x <- data.frame(chr = seqnames(winGR[x]),
                                    start = start(winGR[x]),
                                    end = end(winGR[x]),
                                    midpoint = round((start(winGR[x])+end(winGR[x]))/2),

                                    fk_kappa_rev = NaN,
                                    fk_pval_rev = NaN,
                                    fk_zstat_rev = NaN,
                                    fk_num_reads_rev = NaN,
                                    fk_num_Cs_rev = NaN,
                                    fk_prop_reads_rev = NaN,
                                    fk_prop_Cs_rev = NaN,
 
#                                    pamk_nc_rev = NaN,

                                    mean_stocha_rev = NaN,
                                    sd_stocha_rev = NaN,
                                    mean_min_acf_rev = NaN,
                                    mean_max_acf_rev = NaN,
                                    mean_mean_acf_rev = NaN)

    }
 
    # Make data.frame with relevant info for genomic window
    fk_df_win_x <- data.frame(chr = seqnames(winGR[x]),
                              start = start(winGR[x]),
                              end = end(winGR[x]),
                              midpoint = round((start(winGR[x])+end(winGR[x]))/2),
                              
                              fk_df_fwd_win_x[,5:ncol(fk_df_fwd_win_x)],
                              fk_df_rev_win_x[,5:ncol(fk_df_rev_win_x)],

                              fk_kappa_all = mean(c(fk_df_fwd_win_x$fk_kappa_fwd, fk_df_rev_win_x$fk_kappa_rev), na.rm = T),
                              fk_pval_all = mean(c(fk_df_fwd_win_x$fk_pval_fwd, fk_df_rev_win_x$fk_pval_rev), na.rm = T),
                              fk_zstat_all = mean(c(fk_df_fwd_win_x$fk_zstat_fwd, fk_df_rev_win_x$fk_zstat_rev), na.rm = T),
                              fk_num_reads_all = mean(c(fk_df_fwd_win_x$fk_num_reads_fwd, fk_df_rev_win_x$fk_num_reads_rev), na.rm = T),
                              fk_num_Cs_all = mean(c(fk_df_fwd_win_x$fk_num_Cs_fwd, fk_df_rev_win_x$fk_num_Cs_rev), na.rm = T),
                              fk_prop_reads_all = mean(c(fk_df_fwd_win_x$fk_prop_reads_fwd, fk_df_rev_win_x$fk_prop_reads_rev), na.rm = T),
                              fk_prop_Cs_all = mean(c(fk_df_fwd_win_x$fk_prop_Cs_fwd, fk_df_rev_win_x$fk_prop_Cs_rev), na.rm = T),

#                              pamk_nc_all = mean(c(fk_df_fwd_win_x$pamk_nc_fwd, fk_df_rev_win_x$pamk_nc_rev), na.rm = T)

                              mean_stocha_all = mean(c(fk_df_fwd_win_x$mean_stocha_fwd, fk_df_rev_win_x$mean_stocha_rev), na.rm = T),
                              sd_stocha_all = mean(c(fk_df_fwd_win_x$sd_stocha_fwd, fk_df_rev_win_x$sd_stocha_rev), na.rm = T),
                              mean_min_acf_all = mean(c(fk_df_fwd_win_x$mean_min_acf_fwd, fk_df_rev_win_x$mean_min_acf_rev), na.rm = T),
                              mean_max_acf_all = mean(c(fk_df_fwd_win_x$mean_max_acf_fwd, fk_df_rev_win_x$mean_max_acf_rev), na.rm = T),
                              mean_mean_acf_all = mean(c(fk_df_fwd_win_x$mean_mean_acf_fwd, fk_df_rev_win_x$mean_mean_acf_rev), na.rm = T))

    fk_df_win_x
#  })
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
                               
  fk_df <- dplyr::bind_rows(fk_df_win_list, .id = "column_label")

  fk_df <- data.frame(fk_df,
                      fk_adj_pval_fwd = p.adjust(fk_df$fk_pval_fwd, method = "BH"),
                      fk_adj_pval_rev = p.adjust(fk_df$fk_pval_rev, method = "BH"),
                      fk_adj_pval_all = p.adjust(fk_df$fk_pval_all, method = "BH"))

  write.table(fk_df,
              file = paste0(outDir,
                            sampleName, "_MappedOn_", refbase, "_", context,
                            "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                            "_NAmax", NAmax, "_per_read_var_df_", chrName[i], ".tsv"),
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

