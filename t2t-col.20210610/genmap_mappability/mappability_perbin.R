#!/applications/R/R-4.0.0/bin/Rscript

# Calculate mean mappability in adjacent windows
# for generating chromosome-scale plots

# Usage:
# csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./mappability_perbin.R 'map_K40_E2' t2t-col.20210610 10000 101"

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

# Load mappability track in per-base format
# (containing expanded mappability values where
# they are equal at consecutive bases;
# these are collapsed in bedGraph format)
MAP <- read.table(paste0(sampleName, "/t2t-col.20210610.genmap.sorted.pb"), header = F)

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/", refbase, ".fa.fai"), header = F)
fai <- fai[1:5,]
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])
} else {
  chrs <- fai[,1]
}
chrLens <- fai[,2]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# For each genomeBinSize-bp adjacent window, calculate the mean of the genomeBinSize
# per-base mappability within that window
print(genomeBinName)
cumMAPDat <- NULL
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

  # Define mappability coordinates as GRanges objects
  # and calculate mean mappability in each window
  # MAP
  chrMAP <- MAP[MAP[,1] == chrs[i],]
  chrMAP_GR <- GRanges(seqnames = chrs[i],
                       ranges = IRanges(start = chrMAP[,2],
                                        width = 1),
                       strand = "*")
  MAPoverlaps <- getOverlaps(coordinates = winGR,
                             segments = chrMAP_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  MAPwinVals <- sapply(MAPoverlaps, function(x) mean(as.numeric(chrMAP[,3][x])))
  # Combine in data frame
  MAPDat <- data.frame(chr = as.character(chrs[i]),
                       window = as.integer(start(winGR)),
                       cumwindow = as.integer(start(winGR) + sumchr[i]),
                       mappability  = as.numeric(MAPwinVals))
  cumMAPDat <- rbind(cumMAPDat, MAPDat)
}
write.table(cumMAPDat,
            file = paste0(sampleName, "/", sampleName, "_genmap_MappedOn_", refbase,
                          "_binSize", genomeBinName, "_unsmoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

chrMAPDat <- lapply(seq_along(chrs), function(x) {
  cumMAPDat[cumMAPDat$chr == chrs[x],]
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

filt_MAPDat <- NULL
for(i in seq_along(chrs)) {
  # MAP
  filt_chr_MAP <- stats::filter(chrMAPDat[[i]]$mappability,
                                filter = f,
                                sides = 2)
  filt_chr_MAP[1:flank] <- filt_chr_MAP[flank+1]
  filt_chr_MAP[(length(filt_chr_MAP)-flank+1):length(filt_chr_MAP)] <- filt_chr_MAP[(length(filt_chr_MAP)-flank)]
  # Combine in data frame
  filt_chrMAPDat <- data.frame(chr = chrMAPDat[[i]]$chr,
                               window = as.integer(chrMAPDat[[i]]$window),
                               cumwindow = as.integer(chrMAPDat[[i]]$cumwindow),
                               filt_mappability  = as.numeric(filt_chr_MAP))
  filt_MAPDat <- rbind(filt_MAPDat, filt_chrMAPDat)
}
write.table(filt_MAPDat,
            file = paste0(sampleName, "/", sampleName, "_genmap_MappedOn_", refbase,
                          "_binSize", genomeBinName, "_smoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
