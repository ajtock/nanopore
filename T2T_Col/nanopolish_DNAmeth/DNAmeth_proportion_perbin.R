#!/applications/R/R-4.0.0/bin/Rscript

# This R script is called by ../Snakefile
# Calculate mean DNA methylation proportions in adjacent windows
# for generating chromosome-scale plots

# Usage:
# /applications/R/R-4.0.0/bin/Rscript ./DNAmeth_proportion_perbin.R WT_nanopolishDNAmeth_95_10kb T2T_Col 10000 101

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
genomeBinSize <- as.integer(args[3])
# smoothN must be an odd number
# see profile smoothing below
smoothN <- as.numeric(args[4])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(segmentSeq)

# Load DNA methylation proportion table
CG <- read.table(gzfile(paste0(sampleName,
                               "_MappedOn_", refbase, "_CpG.tsv")),
                 colClasses = c(rep(NA, 2), rep("NULL", 4), NA, "NULL"),
                 header = T)

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/T2T_Col/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1][1:5])
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# For each genomeBinSize-bp adjacent window, calculate the mean of the genomeBinSize
# per-base DNA methylation proportions within that window
print(genomeBinName)
cumMethDat <- NULL
for(i in seq_along(chrs)) {
  # Define adjacent windows
  winSeq <- seq(from = 1, to = chrLens[i], by = genomeBinSize)
  winCum <- winSeq + sumchr[i]
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

  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  # CG
  chrCG <- CG[CG[,1] == chrs[i],]
  chrCG_GR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = chrCG[,2],
                                       width = 1),
                      strand = "*")
  CGoverlaps <- getOverlaps(coordinates = winGR,
                            segments = chrCG_GR,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE,
                            ignoreStrand = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chrCG[,3][x])))
  # Combine in data frame
  methDat <- data.frame(chr = as.character(chrs[i]),
                        window = as.integer(start(winGR)),
                        cumwindow = as.integer(start(winGR) + sumchr[i]),
                        mCG  = as.numeric(CGwinVals))
  cumMethDat <- rbind(cumMethDat, methDat)
}
write.table(cumMethDat,
            file = paste0(
                          "DNAmeth_", sampleName, "_MappedOn_", refbase,
                          "_dedup_binSize", genomeBinName, "_unsmoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

chrMethDat <- lapply(seq_along(chrs), function(x) {
  cumMethDat[cumMethDat$chr == chrs[x],]
})

# Calculate moving average (MA) of current window,
# (smoothN/2)-0.5 previous windows (where smoothN is odd),
# and
# (smoothN/2)-0.5 subsequent windows (where smoothN is odd)
# The higher smoothN is, the greater the smoothing
# Use modulo operation (%%) to confirm that the remainder of
# integer division smoothN %/% 2 is 1 (i.e., smoothN is odd)
stopifnot(smoothN %% 2 == 1)
flank <- smoothN %/% 2
# Define MA filter coefficients
f <- rep(x = 1/smoothN, times = smoothN)

filt_methDat <- NULL
for(i in seq_along(chrs)) {
  # mCG
  filt_chr_mCG <- stats::filter(chrMethDat[[i]]$mCG,
                                filter = f,
                                sides = 2)
  filt_chr_mCG[1:flank] <- filt_chr_mCG[flank+1]
  filt_chr_mCG[(length(filt_chr_mCG)-flank+1):length(filt_chr_mCG)] <- filt_chr_mCG[(length(filt_chr_mCG)-flank)]
  # Combine in data frame
  filt_chrMethDat <- data.frame(chr = chrMethDat[[i]]$chr,
                                window = as.integer(chrMethDat[[i]]$window),
                                cumwindow = as.integer(chrMethDat[[i]]$cumwindow),
                                filt_mCG  = as.numeric(filt_chr_mCG))
  filt_methDat <- rbind(filt_methDat, filt_chrMethDat)
}
write.table(filt_methDat,
            file = paste0(
                          "DNAmeth_", sampleName, "_MappedOn_", refbase,
                          "_dedup_binSize", genomeBinName, "_smoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
