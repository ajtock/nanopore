#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion.R
# csmit -m 300G -c 48 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CHG"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CpG"

# Calculate the per-read methylation proportion in windows along each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#winSize <- 1000
#context <- "CpG"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
winSize <- as.integer(args[3])
context <- args[4]

if(floor(log10(winSize)) + 1 < 4) {
  winName <- paste0(winSize, "bp")
} else if(floor(log10(winSize)) + 1 >= 4 &
          floor(log10(winSize)) + 1 <= 6) {
  winName <- paste0(winSize/1e3, "kb")
} else if(floor(log10(winSize)) + 1 >= 7) {
  winName <- paste0(winSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(parallel)
 
# Read in the raw output .tsv file from Deepsignal methylation model
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                         sampleName, "_MappedOn_", refbase, "_", context, "_raw.tsv"))
 
# Generate a character vector of the unique readIDs in the file
#readID_list <- NULL
readIDs <- unique(tab$V5)

# Loop within parallelised loop to calculate the per-read methylation mean
# across sequential adjacent winSize-nt windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
  y <- tab[tab[,5] == x,]                # Get rows (cytosine positions) for read x
  Chr <- c(y$V1[1])                    # Get the chromosome that read x aligned to
  y <- y[order(y$V2),]                 # Order the rows by ascending position in the chromosome
   
  # Define window start and end coordinates within read
  winStarts <- seq(from = min(y[,2]),
                   to = max(y[,2]),
                   by = winSize)
  if(winStarts[length(winStarts)] == max(y[,2])) {
    winStarts <- winStarts[-length(winStarts)]
  }
  if(max(y[,2]) - min(y[,2]) + 1 >= winSize) {
    winEnds <- seq(from = winStarts[1] + winSize - 1,
                   to = max(y[,2]),
                   by = winSize)
    if(winEnds[length(winEnds)] != max(y[,2])) {
      winEnds <- c(winEnds, max(y[,2]))
    }
  } else {
    winEnds <- max(y[,2])
  }

  filt_methDat <- NULL
  for(i in 1:length(winStarts)) {
    chunk <- y[y[,2] >= winStarts[i] &
               y[,2] <= winEnds[i],]
    # Proceed only if chunk contains rows (cytosine positions)
    if(dim(chunk)[1] > 0) {
      midpoint <- c((max(chunk[,2])+min(chunk[,2]))/2)
      per_read_methyl_freq <- mean(chunk[,9])
      methDat_freq <- data.frame(chr = Chr,
                                 midpoint = midpoint,
                                 per_read_methyl_freq = per_read_methyl_freq,
                                 stringsAsFactors = F)
      filt_methDat <- rbind(filt_methDat, methDat_freq) 
    }
  }
 
  # Andy to Matt: Removed return(filt_methDat) as return() works only within a function (this may be the source of the issue you mentioned);
  #               Replaced with filt_methDat :
  filt_methDat
}, mc.cores = detectCores()))
print("I'm done!")

per_read_DNAmeth_DF <- per_read_DNAmeth_DF[
                         order(per_read_DNAmeth_DF[,1], per_read_DNAmeth_DF[,2]),
                       ]

write.table(per_read_DNAmeth_DF,
            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
                          "_raw_winSize", winName, "_per_read_midpoint.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
