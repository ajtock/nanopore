#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.02.2021

# Calculate windowed TE frequencies along each T2T_Col chromosome

# Usage:
# /applications/R/R-4.0.0/bin/Rscript TEsuperfam_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 101

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#genomeBinSize <- 10000
#maPeriod <- 101

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
genomeBinSize <- as.integer(args[2])
maPeriod <- as.integer(args[3])

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
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Get TE superfamNames
superfamNames <- system(paste0("ls *", paste0(chrName, collapse = "_"), ".bed"), intern = T)
superfamNames <- gsub("T2T_Col_TEs_", "", superfamNames)
superfamNames <- gsub(paste0("_", paste0(chrName, collapse = "_"), ".bed"), "", superfamNames)

superfamNames <- sort(unique(superfamNames))

superfamDFlist <- lapply(seq_along(superfamNames), function(x) {
  tmp <- read.table(paste0("T2T_Col_TEs_", superfamNames[x], "_",
                           paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
  colnames(tmp) <- c("chr", "start0based", "end", "name", "featureID", "strand")
  tmp 
})

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

for(k in seq_along(superfamNames)) {
  superfamkProfile <- NULL
  for(i in 1:length(chrs)) {
    # Count features within adjacent windows
    featuresChr <- superfamDFlist[[k]][superfamDFlist[[k]]$chr == chrs[i],]
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
                             features = as.integer(winfeatures),
                             stringsAsFactors = F)
    superfamkProfile <- rbind(superfamkProfile, profileChr)
  }
  write.table(superfamkProfile,
              file = paste0("T2T_Col_TEs_", superfamNames[k],
                            "_frequency_per_", genomeBinName, "_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), "_unsmoothed.tsv"),
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
    superfamkProfile[superfamkProfile$chr == chrs[x],]
  }, mc.cores = length(chrs))
  
  filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
    filt_superfamkProfile <- NULL
    for(k in seq_along(superfamNames)) {
      filt_chrProfilek <- stats::filter(x = chrProfiles[[x]][,k+2],
                                        filter = f,
                                        sides = 2)
      filt_chrProfilek[1:flank] <- filt_chrProfilek[flank+1]
      filt_chrProfilek[(length(filt_chrProfilek)-flank+1):length(filt_chrProfilek)] <- filt_chrProfilek[(length(filt_chrProfilek)-flank)]
      filt_superfamkProfile <- data.frame(chr = as.character(chrProfiles[[x]]$chr),
                                          window = as.integer(chrProfiles[[x]]$window),
                                          filt_features = as.numeric(filt_chrProfilek),
                                          stringsAsFactors = F)
      filt_superfamkProfile <- cbind(filt_superfamkProfile, filt_superfamkProfile[,3])
    }
    filt_superfamkProfile <- data.frame(filt_superfamkProfile[,1:2],
                                        filt_superfamkProfile,
                                        stringsAsFactors = F)
    colnames(filt_superfamkProfile) <- c("chr", "window",
                                         paste0("filt_quantile", 1:quantiles))
    return(filt_superfamkProfile)
  }, mc.cores = length(chrProfiles))
  
  # Combine list of 1 data.frame per chromosome into one data.frame
  filt_superfamkProfile <- do.call(rbind, filt_chrProfiles)
  write.table(filt_superfamkProfile,
              file = paste0(outDir,
                            "CEN180_frequency_per_", genomeBinName,
                            "_", quantiles, "quantiles_",
                            "_by_", orderingFactor,
                            "_of_CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), "_smoothed.tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
