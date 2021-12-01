#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.02.2021

# Calculate windowed TE frequencies along each t2t-col.20210610 chromosome

# Usage:
# /applications/R/R-4.0.0/bin/Rscript gene_frequency_chrProfiles.R 10000 101

#genomeBinSize <- 10000
#maPeriod <- 101

args <- commandArgs(trailingOnly = T)
genomeBinSize <- as.integer(args[1])
maPeriod <- as.integer(args[2])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14841110,3823792,13597188,4203902,11784131)
CENend <- c(17559778,6045243,15733925,6977949,14551809)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

genes <- read.table(paste0("t2t-col.20210610_representative_mRNA_Chr1_Chr2_Chr3_Chr4_Chr5.bed"),
                     header = F)
colnames(genes) <- c("chr", "start0based", "end", "name", "score", "strand")

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = genomeBinSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = genomeBinSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

featureProfile <- NULL
for(i in 1:length(chrs)) {
  # Count features within adjacent windows
  featuresChr <- genes[genes$chr == chrs[i],]
  windowsChrGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  if(dim(featuresChr)[1] > 0) {
    featuresChrGR <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = featuresChr$start+1,
                                              end = featuresChr$end),
                             strand = featuresChr$strand)
  } else {
    featuresChrGR <- GRanges()
  }
  winfeatures <- countOverlaps(query = windowsChrGR,
                               subject = featuresChrGR,
                               type = "any",
                               ignore.strand = T)
  profileChr <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(windowsChrGR)),
                           cumwindow = as.integer(start(windowsChrGR)+sumchr[i]),
                           features = as.integer(winfeatures),
                           stringsAsFactors = F)
  featureProfile <- rbind(featureProfile, profileChr)
}
write.table(featureProfile,
            file = paste0("t2t-col.20210610_gene_frequency_per_", genomeBinName, "_unsmoothed.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")

# Calculate moving average of current window,
#### (maPeriod/2) previous windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 previous windows (where maPeriod is odd),
# and
#### (maPeriod/2) subsequent windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 subsequent windows (where maPeriod is odd)
# (the higher maPeriod is, the greater the smoothing)
stopifnot(maPeriod %% 2 != 0)
flank <- (maPeriod/2)-0.5
# Define MA filter coefficients
f <- rep(1/maPeriod, maPeriod)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  featureProfile[featureProfile$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_chrProfilek <- stats::filter(x = chrProfiles[[x]]$features,
                                    filter = f,
                                    sides = 2)
  filt_chrProfilek[1:flank] <- filt_chrProfilek[flank+1]
  filt_chrProfilek[(length(filt_chrProfilek)-flank+1):length(filt_chrProfilek)] <- filt_chrProfilek[(length(filt_chrProfilek)-flank)]
  filt_featureProfile <- data.frame(chr = as.character(chrProfiles[[x]]$chr),
                                    window = as.integer(chrProfiles[[x]]$window),
                                    cumwindow = as.integer(chrProfiles[[x]]$cumwindow),
                                    filt_features = as.numeric(filt_chrProfilek),
                                    stringsAsFactors = F)
  return(filt_featureProfile)
}, mc.cores = length(chrProfiles))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_featureProfile <- do.call(rbind, filt_chrProfiles)
write.table(filt_featureProfile,
            file = paste0("t2t-col.20210610_gene_frequency_per_", genomeBinName, "_smoothed.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")
