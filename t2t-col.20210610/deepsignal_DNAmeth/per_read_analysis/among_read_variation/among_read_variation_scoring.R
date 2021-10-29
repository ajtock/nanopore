#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 1 CHG 0.30 Chr1"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript among_read_variation_scoring.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 1 CpG 0.20 Chr1" 

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#readBinCs <- 20
#genomeBinSize <- 1000
#genomeStepSize <- 1
#context <- "CpG"
#CPUpc <- 0.20
#chrName <- unlist(strsplit("Chr4", split = ","))

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
genomeBinSize <- as.integer(args[3])
genomeStepSize <- as.integer(args[4])
context <- args[5]
CPUpc <- as.numeric(args[6])
chrName <- unlist(strsplit(args[7], split = ","))

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
library(irr)
library(dplyr)
library(tidyr)
#library(data.table)
#library(segmentSeq)
#library(ComplexHeatmap)
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
#  # Define sliding windows of width genomeBinSize bp,
#  # with a step of genomeStepSize vp
#  ## Note: the active code creates windows of genomeBinSize bp only,
#  ## whereas the commented-out code creates windows decreasing from genomeBinSize bp to genomeStepSize bp
#  ## at the right-hand end of each chromosome ( from chrLens[x]-genomeBinSize to chrLens[x] ),
#  winStarts <- seq(from = 1,
##                   to = chrLens[i],
#                   to = chrLens[i]-genomeBinSize,
#                   by = genomeStepSize)
##  stopifnot(winStarts[length(winStarts)] == chrLens[i])
#  if(chrLens[i] - winStarts[length(winStarts)] >= genomeBinSize) {
#    winStarts <- c(winStarts,
#                   winStarts[length(winStarts)]+genomeStepSize)
#  }
#  winEnds <- seq(from = winStarts[1]+genomeBinSize-1,
#                 to = chrLens[i],
#                 by = genomeStepSize)
#  stopifnot(winEnds[length(winEnds)] == chrLens[i])
#  winEnds <- c(winEnds,
#               rep(chrLens[i], times = length(winStarts)-length(winEnds)))
#  stopifnot(length(winStarts) == length(winEnds))
#
#  winGR <- GRanges(seqnames = chrs[i],
#                   ranges = IRanges(start = winStarts,
#                                    end = winEnds),
#                   strand = "*")
#  print(winGR)

  # Define adjacent windows
  winSeq <- seq(from = 1, to = chrLens[i], by = genomeBinSize)
  winIR <- IRanges(start = winSeq,
                   width = genomeBinSize)
  winIR <- winIR[-length(winIR)]
  winIR <- append(winIR,
                  IRanges(start = winSeq[length(winSeq)],
                          end = chrLens[i]))
  winGR <- GRanges(seqnames = chrs[i],
                   ranges = winIR,
                   strand = "*")
  print(winGR)

  # Define per-read-window midpoint coordinates as GRanges objects
  # and get corresponding DNA methylation proportions that overlap genomic windows
#  chr_tab <- tab[tab[,1] == chrs[i],]
  chr_tab <- tab
  chr_tab_GR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chr_tab[,2],
                                         width = 1),
                        strand = chr_tab[,3],
                        read = chr_tab[,5],
                        call = chr_tab[,9])
  chr_tab_GR <- sort(chr_tab_GR, by = ~ seqnames + read + start + end + strand)

  chr_tab_GR_fwd <- chr_tab_GR[strand(chr_tab_GR) == "+"]
  chr_tab_GR_fwd <- sort(chr_tab_GR_fwd, by = ~ read + start)
  chr_tab_GR_rev <- chr_tab_GR[strand(chr_tab_GR) == "-"]
  chr_tab_GR_rev <- sort(chr_tab_GR_rev, by = ~ read + start)

  ## Note: findOverlaps() approach does not work where a window does not overlap
  ##       any positions in chr_tab_GR, which can occur with smaller genomeBinSize
  # Identify overlapping windows and midpoint coordinates
  fOverlaps <- findOverlaps(query = winGR,
                            subject = chr_tab_GR,
                            type = "any",
                            select = "all",
                            ignore.strand = TRUE)

  for(x in seq_along(winGR)) {
    chr_tab_GR_x <- chr_tab_GR[subjectHits(fOverlaps[queryHits(fOverlaps) == x])]
    chr_tab_GR_x <- sortSeqlevels(chr_tab_GR_x)
    chr_tab_GR_x <- sort(chr_tab_GR_x, by = ~ seqnames + read + start + end + strand)

    df_x <- data.frame(pos = start(chr_tab_GR_x),
                       read = chr_tab_GR_x$read,
                       call = chr_tab_GR_x$call)

#    spread_x <- tidyr::spread(data = df_x,
#                              key = read,
#                              value = call)
##                              sep = "_")
#    spread_x <- spread_x[ with(data = spread_x, expr = order(pos)), ]
   
    pwider_x <- as.data.frame(tidyr::pivot_wider(data = df_x,
                                                 names_from = read,
#                                                 names_prefix = "read_",
                                                 values_from = call))
    pwider_x <- pwider_x[ with(data = pwider_x, expr = order(pos)), ]
    rownames(pwider_x) <- pwider_x[,1]
    pwider_x <- pwider_x[,-1]
    
#    rownames(pwider_x) <- rownames(spread_x)
#    stopifnot(all.equal(pwider_x, spread_x, check.attributes=F))
    which(colSums(is.na(pwider_x)) == 0)
    which(rowSums(is.na(pwider_x)) == 0)

    tmp <- pwider_x[ , -which(colSums(is.na(pwider_x)) > nrow(pwider_x)*0.5) ]
    which(rowSums(is.na(tmp)) == 0)
    
    kappam.fleiss(tmp)


#  # Convert fOverlaps into list object equivalent to that
#  # generated by segmentSeq::getOverlaps(), in which each
#  # list element corresponds to a vector of sequentially numbered indices of
#  # read midpoint coordinates that overlap a given genomic window
#  fOverlapsList <- mclapply(seq_along(unique(queryHits(fOverlaps))),
#                            function(x) {
#                              subjectHits(fOverlaps[queryHits(fOverlaps) == x])
#                            }, mc.cores = detectCores(), mc.preschedule = TRUE)
  fOverlapsList <- getOverlaps(coordinates = winGR,
                               segments = chr_tab_GR,
                               overlapType = "overlapping",
                               whichOverlaps = TRUE,
                               ignoreStrand = TRUE)

  # Get per-read-window methylation proportion values overlapping each genomic window
  win_mProp_list <- lapply(fOverlapsList, function(x) {
                      data.frame(matrix(data = chr_tab[,3][x], nrow = 1))
                    })
  # Convert into matrix in which each column corresponds to a genomic window
  win_mProp_matrix <- t(as.matrix(x = bind_rows(win_mProp_list)))
  colnames(win_mProp_matrix) <- round(start(winGR)/1e6, digits = 1)
  # Remove columns where fewer than 2 rows are not NA
  win_mProp_matrix <- win_mProp_matrix[,which(colSums(is.na(win_mProp_matrix)) < nrow(win_mProp_matrix) - 1)]  

# Generate a character vector of the unique readIDs in the file
readIDs <- unique(tab$V5)

# Loop within parallelised loop to calculate the per-read methylation proportion
# across sequential adjacent readBinCs-Cs-containing windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
  y <- tab[tab[,5] == x,]              # Get rows (cytosine positions) for read x
  y <- y[order(y$V2, decreasing = F),] # Order the rows by ascending position in the chromosome
 
  # Define window start coordinates within read
  winStarts <- seq(from = 1,
                   to = nrow(y),
                   by = readBinCs)
  # Remove the last winStart value if is the same as the number of rows (total number of Cs in read)
  if(winStarts[length(winStarts)] == nrow(y)) {
    winStarts <- winStarts[-length(winStarts)]
  }
  # Remove the last winStart value if there are fewer than readBinCs Cs from
  # this value to the last C in the read (the last row), so that the last window
  # always has as much or more methylation-state information than the other windows
  tryCatch(
    {
      if(length(winStarts) > 1 && nrow(y) - winStarts[length(winStarts)] + 1 < readBinCs) {
        winStarts <- winStarts[-length(winStarts)]
      }
    },
    error=function(cond) {
      message(paste(x, "read is problematic for winStarts"))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )

  # Define window end coordinates within read
  tryCatch(
    {
      if(nrow(y) >= readBinCs) {
        winEnds <- seq(from = readBinCs,
                       to = nrow(y),
                       by = readBinCs)
        if(winEnds[length(winEnds)] != nrow(y)) {
          winEnds <- c(winEnds, nrow(y))
        }
        # Remove the penultimate winEnd value if there are fewer than readBinCs Cs from
        # this value to the last C in the read (the last row), so that the last window
        # always has as much methylation-state information as, or more than, the other windows
        if(nrow(y) - winEnds[(length(winEnds) - 1)] < readBinCs) {
          winEnds <- winEnds[-(length(winEnds) - 1)]
        }
      } else {
        winEnds <- nrow(y)
      }
    },
    error=function(cond) {
      message(paste(x, "read is problematic for winEnds"))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )

  tryCatch(
    {
      stopifnot(length(winStarts) == length(winEnds))
    },
    error = function(cond) {
      message(paste(x, "read: length of winStarts is not equal to length of winEnds"))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )

  methDat <- NULL
  for(i in 1:length(winStarts)) {
    readwin <- y[winStarts[i] : winEnds[i],]
    # Proceed only if readwin contains rows (cytosine positions)
    if(dim(readwin)[1] > 0) {
      midpoint <- ( min(readwin[,2]) + max(readwin[,2]) ) / 2
      per_readwin_methyl_mean <- mean(readwin[,9])
      per_readwin_methyl_sd <- sd(readwin[,9])
      methDat_mean <- data.frame(chr = readwin[,1][1],
                                 midpoint = midpoint,
                                 per_readwin_methyl_mean = per_readwin_methyl_mean,
                                 per_readwin_methyl_sd = per_readwin_methyl_sd,
                                 start = min(readwin[,2]),
                                 end = max(readwin[,2]),
                                 read = readwin[,5][1],
                                 stringsAsFactors = F)

      methDat <- rbind(methDat, methDat_mean)
    }
  }
 
  # Andy to Matt: Removed return(methDat) as return() works only within a function
  #               (this may be the source of the issue you mentioned);
  #               Replaced with methDat :
  methDat
  # Due to the large numbers of reads that are forked to each CPU, these are
  # RAM-heavy combined processes, which calls for using fewer CPUs than are
  # available on a given node (e.g., 30% of the available CPUs)
}, mc.cores = detectCores()*CPUpc, mc.preschedule = T))

print("Done")

per_read_DNAmeth_DF <- per_read_DNAmeth_DF[
                         order(per_read_DNAmeth_DF[,1], per_read_DNAmeth_DF[,2]),
                       ]

write.table(per_read_DNAmeth_DF,
            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
                          "_raw_readBinSize", readBinCs, "Cs_per_readWin_midpoint.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
