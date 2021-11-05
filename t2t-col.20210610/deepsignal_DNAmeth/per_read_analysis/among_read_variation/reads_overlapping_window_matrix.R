#!/applications/R/R-4.0.0/bin/Rscript

# Make matrices in which columns correspond to cytosines within a
# genomeBinSize window and rows correspond to ONT reads with
# DeepSignal-called methylation status
 
# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript reads_overlapping_window_matrix.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 200000 200000 CHG 0.50 1.00 Chr1"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript reads_overlapping_window_matrix.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 200000 200000 CpG 0.50 1.00 Chr1"

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 200000
#genomeStepSize <- 200000
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
library(dplyr)
 
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
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

CEN <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.centromeres"), header = T)
CEN <- CEN[which(CEN$chr %in% chrName),]
CENstart <- CEN$start
CENend <- CEN$end
CENGR <- GRanges(seqnames = CEN$chr,
                 ranges = IRanges(start = CEN$start,
                                  end = CEN$end),
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

# For each genomeBinSize-bp window with a step of genomeStepSize-bp,
# calculate Fleiss' kappa statistic as a measure of among-read variation
# in methylation state
print(genomeBinName)
print(genomeStepName)
for(i in seq_along(chrName)) {
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

  winGR <- GRanges(seqnames = chrName[i],
                   ranges = IRanges(start = winStarts,
                                    end = winEnds),
                   strand = "*")
  print(winGR)

  # Define GRanges object containing each cytosine and methylation status call
  chr_tab <- tab[tab[,1] == chrName[i],]
  chr_tab_GR <- GRanges(seqnames = chrName[i],
                        ranges = IRanges(start = chr_tab[,2],
                                         width = 1),
                        strand = chr_tab[,3],
                        read = chr_tab[,5],
                        call = chr_tab[,9])

  chr_tab_GR_fwd <- chr_tab_GR[strand(chr_tab_GR) == "+"]
  chr_tab_GR_rev <- chr_tab_GR[strand(chr_tab_GR) == "-"]

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

  # Process each strand separately
  fwd_df_win_list <- mclapply(seq_along(winGR), function(x) {
#  fwd_df_win_list <- lapply(seq_along(winGR), function(x) {
#    print(x)

    chr_tab_GR_fwd_x <- chr_tab_GR_fwd[subjectHits(fOverlaps_fwd[queryHits(fOverlaps_fwd) == x])]

    if(length(chr_tab_GR_fwd_x) > 0) {
      chr_tab_GR_fwd_x <- sortSeqlevels(chr_tab_GR_fwd_x)
      chr_tab_GR_fwd_x <- sort(chr_tab_GR_fwd_x, by = ~ read + start)

      df_fwd_x <- data.frame(pos = start(chr_tab_GR_fwd_x),
                             read = chr_tab_GR_fwd_x$read,
                             call = chr_tab_GR_fwd_x$call)

      pwider_fwd_x <- as.data.frame(tidyr::pivot_wider(data = df_fwd_x,
                                                       names_from = read,
                                                       values_from = call))
      pwider_fwd_x <- t(pwider_fwd_x[ with(data = pwider_fwd_x, expr = order(pos)), ])
      colnames(pwider_fwd_x) <- as.vector(pwider_fwd_x[1, ])
#      pwider_fwd_x <- cbind(rownames(pwider_fwd_x), pwider_fwd_x)
      pwider_fwd_x <- pwider_fwd_x[-1, , drop = F]
#      colnames(pwider_fwd_x)[1] <- "readID"
#      rownames(pwider_fwd_x) <- 1:nrow(pwider_fwd_x)
    } else {
      pwider_fwd_x <- NULL
    }
    
    pwider_fwd_x

#  })
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

  names(fwd_df_win_list) <- paste0(chrName[i], "_", start(winGR), "_", end(winGR))

  save(fwd_df_win_list,
       file = paste0(outDir,
                     sampleName, "_MappedOn_", refbase,
                     "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "_matrix_list_",
                     "_per_read_methylation_calls_fwd_strand_", context, "_", chrName[i], ".RData"))

  rm(fwd_df_win_list); gc()
#  load(paste0(outDir,
#              sampleName, "_MappedOn_", refbase,
#              "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "_matrix_list_",
#              "_per_read_methylation_calls_fwd_strand_", context, "_", chrName[i], ".RData"))

  # Process each strand separately
  rev_df_win_list <- mclapply(seq_along(winGR), function(x) {
#  rev_df_win_list <- lapply(seq_along(winGR), function(x) {
#    print(x)

    chr_tab_GR_rev_x <- chr_tab_GR_rev[subjectHits(fOverlaps_rev[queryHits(fOverlaps_rev) == x])]

    if(length(chr_tab_GR_rev_x) > 0) {
      chr_tab_GR_rev_x <- sortSeqlevels(chr_tab_GR_rev_x)
      chr_tab_GR_rev_x <- sort(chr_tab_GR_rev_x, by = ~ read + start)

      df_rev_x <- data.frame(pos = start(chr_tab_GR_rev_x),
                             read = chr_tab_GR_rev_x$read,
                             call = chr_tab_GR_rev_x$call)

      pwider_rev_x <- as.data.frame(tidyr::pivot_wider(data = df_rev_x,
                                                       names_from = read,
                                                       values_from = call))
      pwider_rev_x <- t(pwider_rev_x[ with(data = pwider_rev_x, expr = order(pos)), ])
      colnames(pwider_rev_x) <- as.vector(pwider_rev_x[1, ])
#      pwider_rev_x <- cbind(rownames(pwider_rev_x), pwider_rev_x)
      pwider_rev_x <- pwider_rev_x[-1, , drop = F]
#      colnames(pwider_rev_x)[1] <- "readID"
#      rownames(pwider_rev_x) <- 1:nrow(pwider_rev_x)
    } else {
      pwider_rev_x <- NULL
    }
    
    pwider_rev_x

#  })
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

  names(rev_df_win_list) <- paste0(chrName[i], "_", start(winGR), "_", end(winGR))

  save(rev_df_win_list,
       file = paste0(outDir,
                     sampleName, "_MappedOn_", refbase,
                     "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "_matrix_list_",
                     "_per_read_methylation_calls_rev_strand_", context, "_", chrName[i], ".RData"))

  rm(rev_df_win_list); gc()
#  load(paste0(outDir,
#              sampleName, "_MappedOn_", refbase,
#              "_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "_matrix_list_",
#              "_per_read_methylation_calls_rev_strand_", context, "_", chrName[i], ".RData"))

}
