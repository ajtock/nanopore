#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Among-read variation/agreement (Fleiss' kappa) [this script]
# 2. Variance of read stochasticity (variance of per-read, per-window proportion of inter-C intervals that represent a methylation state change [variance across all reads that overlap window examined])
# 3. Per-read autocorrelation, but could be hard to derive a general measure across reads overlapping a given region (maybe look at variance of pairwise read correlation coefficients between autocorrelation values, either standardised to e.g. 100 autocorrelation values per read, or use only reads with info across the same Cs) 

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 2000 CHG 0.50 0.20 Chr1"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 2000 CpG 0.50 0.20 Chr1"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 0.20
#chrName <- unlist(strsplit("Chr4", split = ","))

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
 
if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

if(floor(log10(genomeStepSize)) + 1 < 4) {
  genomeStepName <- paste0(genomeStepSize, "bp")
} else if(floor(log10(genomeStepSize)) + 1 >= 4 &
          floor(log10(genomeStepSize)) + 1 <= 6) {
  genomeStepName <- paste0(genomeStepSize/1e3, "kb")
} else if(floor(log10(genomeStepSize)) + 1 >= 7) {
  genomeStepName <- paste0(genomeStepSize/1e6, "Mb")
}

outDir <- paste0("genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[1:5]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Read in the raw output .tsv file from Deepsignal methylation model
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                         sampleName, "_MappedOn_", refbase, "_", context, "_raw_", chrName, ".tsv"),
                  header = F)
 
# For each genomeBinSize-bp window with a step of genomeStepSize-bp,
# calculate Fleiss' kappa statistic as a measure of among-read variation
# in methylation state
print(genomeBinName)
print(genomeStepName)
for(i in seq_along(chrs)) {
  # Define sliding windows of width genomeBinSize bp,
  # with a step of genomeStepSize vp
  ## Note: the active code creates windows of genomeBinSize bp only,
  ## whereas the commented-out code creates windows decreasing from genomeBinSize bp to genomeStepSize bp
  ## at the right-hand end of each chromosome ( from chrLens[x]-genomeBinSize to chrLens[x] ),
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

  winGR <- GRanges(seqnames = chrs[i],
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
#  winGR <- GRanges(seqnames = chrs[i],
#                   ranges = winIR,
#                   strand = "*")
#  print(winGR)

  # Define per-read-window midpoint coordinates as GRanges objects
  # and get corresponding DNA methylation proportions that overlap genomic windows
  chr_tab <- tab[tab[,1] == chrs[i],]
#  chr_tab <- tab
  chr_tab_GR <- GRanges(seqnames = chrs[i],
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
                                ignore.strand = TRUE)

  fOverlaps_rev <- findOverlaps(query = winGR,
                                subject = chr_tab_GR_rev,
                                type = "any",
                                select = "all",
                                ignore.strand = TRUE)

  fk_df_win_list <- mclapply(seq_along(winGR), function(x) {

    # Analyse each strand separately
    # fwd
    chr_tab_GR_fwd_x <- chr_tab_GR_fwd[subjectHits(fOverlaps_fwd[queryHits(fOverlaps_fwd) == x])]
    chr_tab_GR_fwd_x <- sortSeqlevels(chr_tab_GR_fwd_x)
    chr_tab_GR_fwd_x <- sort(chr_tab_GR_fwd_x, by = ~ read + start)

    df_fwd_x <- data.frame(pos = start(chr_tab_GR_fwd_x),
                           read = chr_tab_GR_fwd_x$read,
                           call = chr_tab_GR_fwd_x$call)

#    # tidyr::spread() is deprecated; use tidyr::pivot_wider() instead 
#    spread_fwd_x <- tidyr::spread(data = df_fwd_x,
#                                  key = read,
#                                  value = call)
##                                  sep = "_")
#    spread_fwd_x <- spread_x[ with(data = spread_x, expr = order(pos)), ]
 
    pwider_fwd_x <- as.data.frame(tidyr::pivot_wider(data = df_fwd_x,
                                                     names_from = read,
#                                                     names_prefix = "read_",
                                                     values_from = call))
    pwider_fwd_x <- pwider_fwd_x[ with(data = pwider_fwd_x, expr = order(pos)), ]
    rownames(pwider_fwd_x) <- pwider_fwd_x[,1]
    pwider_fwd_x <- pwider_fwd_x[,-1]

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

    # Remove rows (cytosines) containing any NAs across the retained columns (reads),
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

    # Calculate Fleiss' kappa
    fkappa_pwider_fwd_x <- kappam.fleiss(pwider_fwd_x, detail = F)

    # Sanity checks
    stopifnot(fkappa_pwider_fwd_x$raters == num_reads_retained_fwd_x)
    stopifnot(fkappa_pwider_fwd_x$subjects == num_Cs_retained_fwd_x)

    # Define clusters of reads within each window using fpc::pamk()
    # ("partitioning around medoids with estimation of number of clusters")
    set.seed(20000)
    pamk_pwider_fwd_x <- pamk(t(pwider_fwd_x),
#                              krange = 1:(nrow(t(pwider_fwd_x))-1),
                              krange = 1:(round(nrow(t(pwider_fwd_x))/2)),
                              criterion = "asw",
                              usepam = FALSE,
                              scaling = FALSE,
                              alpha = 0.001,
#                              ns = 10,
#                              seed = 100000,
                              diss = FALSE)

  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)

#    htmp <- Heatmap(t(as.matrix(pwider_fwd_x)),
#                    col = c("0" = "blue", "1" = "red"),
#                    row_split = paste0("Cluster", pamk_pwider_fwd_x$pamobject$clustering),
#                    show_column_dend = FALSE, 
#                    cluster_columns = FALSE,
#                    heatmap_legend_param = list(title = context,
#                                                title_position = "topcenter",
#                                                title_gp = gpar(font = 2, fontsize = 12),
#                                                legend_direction = "horizontal",
#                                                labels_gp = gpar(fontsize = 10)),
#                    column_title = paste0(chrName, ":", start(winGR[x]), "-" , end(winGR[x])))
#    pdf(paste0(plotDir,
#               sampleName, "_MappedOn_", refbase, "_", context,
#               "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#               "_NAmax", NAmax, "_", chrName, "_", start(winGR[x]), "_", end(winGR[x]), "_fwd",
#               ".pdf"),
#        height = 10, width = 50)
#    draw(htmp,
#         heatmap_legend_side = "bottom",
#         gap = unit(c(1), "mm"))
#    dev.off()


    # Analyse each strand separately
    # rev
    chr_tab_GR_rev_x <- chr_tab_GR_rev[subjectHits(fOverlaps_rev[queryHits(fOverlaps_rev) == x])]
    chr_tab_GR_rev_x <- sortSeqlevels(chr_tab_GR_rev_x)
    chr_tab_GR_rev_x <- sort(chr_tab_GR_rev_x, by = ~ read + start)

    df_rev_x <- data.frame(pos = start(chr_tab_GR_rev_x),
                           read = chr_tab_GR_rev_x$read,
                           call = chr_tab_GR_rev_x$call)

#    # tidyr::spread() is deprecated; use tidyr::pivot_wider() instead 
#    spread_rev_x <- tidyr::spread(data = df_rev_x,
#                                  key = read,
#                                  value = call)
##                                  sep = "_")
#    spread_rev_x <- spread_x[ with(data = spread_x, expr = order(pos)), ]
 
    pwider_rev_x <- as.data.frame(tidyr::pivot_wider(data = df_rev_x,
                                                     names_from = read,
#                                                     names_prefix = "read_",
                                                     values_from = call))
    pwider_rev_x <- pwider_rev_x[ with(data = pwider_rev_x, expr = order(pos)), ]
    rownames(pwider_rev_x) <- pwider_rev_x[,1]
    pwider_rev_x <- pwider_rev_x[,-1]

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

    # Remove rows (cytosines) containing any NAs across the retained columns (reads),
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

    # Calculate Fleiss' kappa
    fkappa_pwider_rev_x <- kappam.fleiss(pwider_rev_x, detail = F)

    # Sanity checks
    stopifnot(fkappa_pwider_rev_x$raters == num_reads_retained_rev_x)
    stopifnot(fkappa_pwider_rev_x$subjects == num_Cs_retained_rev_x)

    # Define clusters of reads within each window using fpc::pamk()
    # ("partitioning around medoids with estimation of number of clusters")
    set.seed(20000)
    pamk_pwider_rev_x <- pamk(t(pwider_rev_x),
                              krange = 1:10,
                              criterion = "multiasw",
                              usepam = FALSE,
                              scaling = FALSE,
                              alpha = 0.001,
                              ns = 2,
                              seed = 100000,
                              diss = FALSE)
   
#    htmp <- Heatmap(t(as.matrix(pwider_rev_x)),
#                    col = c("0" = "blue", "1" = "red"),
#                    row_split = paste0("Cluster", pamk_pwider_rev_x$pamobject$clustering),
#                    show_column_dend = FALSE, 
#                    cluster_columns = FALSE,
#                    heatmap_legend_param = list(title = context,
#                                                title_position = "topcenter",
#                                                title_gp = gpar(font = 2, fontsize = 12),
#                                                legend_direction = "horizontal",
#                                                labels_gp = gpar(fontsize = 10)),
#                    column_title = paste0(chrName, ":", start(winGR[x]), "-" , end(winGR[x])))
#    pdf(paste0(plotDir,
#               sampleName, "_MappedOn_", refbase, "_", context,
#               "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
#               "_NAmax", NAmax, "_", chrName, "_", start(winGR[x]), "_", end(winGR[x]), "_rev",
#               ".pdf"),
#        height = 10, width = 50)
#    draw(htmp,
#         heatmap_legend_side = "bottom",
#         gap = unit(c(1), "mm"))
#    dev.off()


    # Make data.frame with relevant info for genomic window
    fk_df_win_x <- data.frame(chr = seqnames(winGR[x]),
                              start = start(winGR[x]),
                              end = end(winGR[x]),
                              midpoint = round(start(winGR[x])+end(winGR[x])/2),
                              fk_kappa_fwd = fkappa_pwider_fwd_x$value,
                              fk_pval_fwd = fkappa_pwider_fwd_x$p.value,
                              fk_zstat_fwd = fkappa_pwider_fwd_x$statistic,
                              fk_num_reads_fwd = fkappa_pwider_fwd_x$raters,
                              fk_num_Cs_fwd = fkappa_pwider_fwd_x$subjects,
                              fk_prop_reads_fwd = prop_reads_retained_fwd_x,
                              fk_prop_Cs_fwd = prop_Cs_retained_fwd_x,
                              fk_kappa_rev = fkappa_pwider_rev_x$value,
                              fk_pval_rev = fkappa_pwider_rev_x$p.value,
                              fk_zstat_rev = fkappa_pwider_rev_x$statistic,
                              fk_num_reads_rev = fkappa_pwider_rev_x$raters,
                              fk_num_Cs_rev = fkappa_pwider_rev_x$subjects,
                              fk_prop_reads_rev = prop_reads_retained_rev_x,
                              fk_prop_Cs_rev = prop_Cs_retained_rev_x,
                              fk_kappa_all = mean(c(fkappa_pwider_fwd_x$value, fkappa_pwider_rev_x$value)),
                              fk_pval_all = mean(c(fkappa_pwider_fwd_x$p.value, fkappa_pwider_rev_x$p.value)),
                              fk_zstat_all = mean(c(fkappa_pwider_fwd_x$statistic, fkappa_pwider_rev_x$statistic)),
                              fk_num_reads_all = mean(c(fkappa_pwider_fwd_x$raters, fkappa_pwider_rev_x$raters)),
                              fk_num_Cs_all = mean(c(fkappa_pwider_fwd_x$subjects, fkappa_pwider_rev_x$subjects)),
                              fk_prop_reads_all = mean(c(prop_reads_retained_fwd_x, prop_reads_retained_rev_x)),
                              fk_prop_Cs_all = mean(c(prop_Cs_retained_fwd_x, prop_Cs_retained_rev_x)),
                              pamk_pwider_fwd_nc = pamk_pwider_fwd_x$nc,
                              pamk_pwider_rev_nc = pamk_pwider_rev_x$nc,
                              pamk_pwider_all_nc = mean(c(pamk_pwider_fwd_x$nc, pamk_pwider_rev_x$nc), na.rm = T))

    fk_df_win_x
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
                               
  fk_df <- dplyr::bind_rows(fk_df_win_list, .id = "column_label")

  fk_df <- data.frame(fk_df,
                      fk_adj_pval_fwd = p.adjust(fk_df$fk_pval_fwd, method = "BH"),
                      fk_adj_pval_rev = p.adjust(fk_df$fk_pval_rev, method = "BH"),
                      fk_adj_pval_all = p.adjust(fk_df$fk_pval_all, method = "BH"))

  pdf(paste0(plotDir,
             sampleName, "_MappedOn_", refbase, "_", context,
             "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
             "_NAmax", NAmax, "_", chrName,
             ".pdf"), height = 5, width = 30)
  par(mfrow = c(1, 1))
  par(mar = c(4.1, 4.1, 3.1, 4.1))
  par(mgp = c(3, 1, 0))

  plot(x = fk_df$midpoint, y = fk_df$fk_kappa_all, type = "l", lwd = 1.5, col = "red",
       yaxt = "n",
       xlab = "", ylab = "",
       main = "")
  mtext(side = 2, line = 2.25, cex = 1.5, col = "red",
        text = bquote("Fleiss' kappa per-read m"*.(context)))
#        text = bquote("Fleiss' kappa per-read m"*.(context)))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)

#  par(new = T)
#  plot(x = fk_df$midpoint, y = fk_df$fk_num_Cs_all, type = "l", lwd = 0.5, col = "skyblue",
#       xaxt = "n", yaxt = "n",
#       xlab = "", ylab = "",
#       main = "")
#  p <- par('usr')
#  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.0), xpd = NA, srt = -90, col = "skyblue",
#       labels = bquote(.(context)*" sites per window"))
#  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)

  mtext(side = 1, line = 2.25, cex = 1.5, text = paste0(chrName, " (", genomeBinName, " window, ", genomeStepName, " step)"))

  dev.off()

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

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 2000
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 0.20
#chrName <- unlist(strsplit("Chr4", split = ","))

#write.table(per_read_DNAmeth_DF,
#            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
#                          "_raw_readBinSize", readBinCs, "Cs_per_readWin_midpoint.tsv"),
#            quote = F, sep = "\t", row.names = F, col.names = T)
}

##  # Convert fOverlaps into list object equivalent to that
##  # generated by segmentSeq::getOverlaps(), in which each
##  # list element corresponds to a vector of sequentially numbered indices of
##  # read midpoint coordinates that overlap a given genomic window
##  fOverlapsList <- mclapply(seq_along(unique(queryHits(fOverlaps))),
##                            function(x) {
##                              subjectHits(fOverlaps[queryHits(fOverlaps) == x])
##                            }, mc.cores = detectCores(), mc.preschedule = TRUE)
#  fOverlapsList <- getOverlaps(coordinates = winGR,
#                               segments = chr_tab_GR,
#                               overlapType = "overlapping",
#                               whichOverlaps = TRUE,
#                               ignoreStrand = TRUE)
#
#  # Get per-read-window methylation proportion values overlapping each genomic window
#  win_mProp_list <- lapply(fOverlapsList, function(x) {
#                      data.frame(matrix(data = chr_tab[,3][x], nrow = 1))
#                    })
#  # Convert into matrix in which each column corresponds to a genomic window
#  win_mProp_matrix <- t(as.matrix(x = bind_rows(win_mProp_list)))
#  colnames(win_mProp_matrix) <- round(start(winGR)/1e6, digits = 1)
#  # Remove columns where fewer than 2 rows are not NA
#  win_mProp_matrix <- win_mProp_matrix[,which(colSums(is.na(win_mProp_matrix)) < nrow(win_mProp_matrix) - 1)]  
#
## Generate a character vector of the unique readIDs in the file
#readIDs <- unique(tab$V5)
#
## Loop within parallelised loop to calculate the per-read methylation proportion
## across sequential adjacent readBinCs-Cs-containing windows along each read
#per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
#  y <- tab[tab[,5] == x,]              # Get rows (cytosine positions) for read x
#  y <- y[order(y$V2, decreasing = F),] # Order the rows by ascending position in the chromosome
# 
#  # Define window start coordinates within read
#  winStarts <- seq(from = 1,
#                   to = nrow(y),
#                   by = readBinCs)
#  # Remove the last winStart value if is the same as the number of rows (total number of Cs in read)
#  if(winStarts[length(winStarts)] == nrow(y)) {
#    winStarts <- winStarts[-length(winStarts)]
#  }
#  # Remove the last winStart value if there are fewer than readBinCs Cs from
#  # this value to the last C in the read (the last row), so that the last window
#  # always has as much or more methylation-state information than the other windows
#  tryCatch(
#    {
#      if(length(winStarts) > 1 && nrow(y) - winStarts[length(winStarts)] + 1 < readBinCs) {
#        winStarts <- winStarts[-length(winStarts)]
#      }
#    },
#    error=function(cond) {
#      message(paste(x, "read is problematic for winStarts"))
#      message("Here's the original error message:")
#      message(cond)
#      # Choose a return value in case of error
#      return(NA)
#    }
#  )
#
#  # Define window end coordinates within read
#  tryCatch(
#    {
#      if(nrow(y) >= readBinCs) {
#        winEnds <- seq(from = readBinCs,
#                       to = nrow(y),
#                       by = readBinCs)
#        if(winEnds[length(winEnds)] != nrow(y)) {
#          winEnds <- c(winEnds, nrow(y))
#        }
#        # Remove the penultimate winEnd value if there are fewer than readBinCs Cs from
#        # this value to the last C in the read (the last row), so that the last window
#        # always has as much methylation-state information as, or more than, the other windows
#        if(nrow(y) - winEnds[(length(winEnds) - 1)] < readBinCs) {
#          winEnds <- winEnds[-(length(winEnds) - 1)]
#        }
#      } else {
#        winEnds <- nrow(y)
#      }
#    },
#    error=function(cond) {
#      message(paste(x, "read is problematic for winEnds"))
#      message("Here's the original error message:")
#      message(cond)
#      # Choose a return value in case of error
#      return(NA)
#    }
#  )
#
#  tryCatch(
#    {
#      stopifnot(length(winStarts) == length(winEnds))
#    },
#    error = function(cond) {
#      message(paste(x, "read: length of winStarts is not equal to length of winEnds"))
#      message("Here's the original error message:")
#      message(cond)
#      # Choose a return value in case of error
#      return(NA)
#    }
#  )
#
#  methDat <- NULL
#  for(i in 1:length(winStarts)) {
#    readwin <- y[winStarts[i] : winEnds[i],]
#    # Proceed only if readwin contains rows (cytosine positions)
#    if(dim(readwin)[1] > 0) {
#      midpoint <- ( min(readwin[,2]) + max(readwin[,2]) ) / 2
#      per_readwin_methyl_mean <- mean(readwin[,9])
#      per_readwin_methyl_sd <- sd(readwin[,9])
#      methDat_mean <- data.frame(chr = readwin[,1][1],
#                                 midpoint = midpoint,
#                                 per_readwin_methyl_mean = per_readwin_methyl_mean,
#                                 per_readwin_methyl_sd = per_readwin_methyl_sd,
#                                 start = min(readwin[,2]),
#                                 end = max(readwin[,2]),
#                                 read = readwin[,5][1],
#                                 stringsAsFactors = F)
#
#      methDat <- rbind(methDat, methDat_mean)
#    }
#  }
# 
#  # Andy to Matt: Removed return(methDat) as return() works only within a function
#  #               (this may be the source of the issue you mentioned);
#  #               Replaced with methDat :
#  methDat
#  # Due to the large numbers of reads that are forked to each CPU, these are
#  # RAM-heavy combined processes, which calls for using fewer CPUs than are
#  # available on a given node (e.g., 30% of the available CPUs)
#}, mc.cores = detectCores()*CPUpc, mc.preschedule = T))
#
#print("Done")
#
#per_read_DNAmeth_DF <- per_read_DNAmeth_DF[
#                         order(per_read_DNAmeth_DF[,1], per_read_DNAmeth_DF[,2]),
#                       ]
#
#write.table(per_read_DNAmeth_DF,
#            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
#                          "_raw_readBinSize", readBinCs, "Cs_per_readWin_midpoint.tsv"),
#            quote = F, sep = "\t", row.names = F, col.names = T)
