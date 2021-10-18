#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winByCs_30pc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10 CHG"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winByCs_30pc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10 CpG"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#noOfCs <- 10
#context <- "CpG"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
noOfCs <- as.integer(args[3])
context <- args[4]

options(stringsAsFactors = F)
library(parallel)
 
# Read in the raw output .tsv file from Deepsignal methylation model
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                         sampleName, "_MappedOn_", refbase, "_", context, "_raw_tail100000000.tsv"),
                  header = F)
 
# Generate a character vector of the unique readIDs in the file
readIDs <- unique(tab$V5)

# Loop within parallelised loop to calculate the per-read methylation proportion
# across sequential adjacent noOfCs-Cs-containing windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
#mclapply(readIDs, function(x) {
#  y <- tab[tab[,5] == x,c(1,2)]
  y <- tab[tab[,5] == x,]              # Get rows (cytosine positions) for read x
  y <- y[order(y$V2, decreasing = F),] # Order the rows by ascending position in the chromosome
 
  # Define window start coordinates within read
  winStarts <- seq(from = 1,
                   to = nrow(y),
                   by = noOfCs)
  # Remove the last winStart value if is the same as the number of rows (total number of Cs in read)
  if(winStarts[length(winStarts)] == nrow(y)) {
    winStarts <- winStarts[-length(winStarts)]
  }
  # Remove the last winStart value if there are fewer than noOfCs Cs from
  # this value to the last C in the read (the last row), so that the last window
  # always has as much or more methylation-state information than the other windows
  tryCatch(
    {
      if(nrow(y) - winStarts[length(winStarts)] + 1 < noOfCs) {
        winStarts <- winStarts[-length(winStarts)]
      }
    },
    error=function(cond) {
      message(paste(x, "read is problematic"))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )

  # Define window end coordinates within read
  if(nrow(y) - winStarts[1] + 1 >= noOfCs) {
    winEnds <- seq(from = winStarts[1] + noOfCs - 1,
                   to = nrow(y),
                   by = noOfCs)
    # Remove the last winEnd value if there are fewer than noOfCs Cs from
    # this value to the last C in the read (the last row), so that the last window
    # always has as much or more methylation-state information than the other windows
    if(nrow(y) - winEnds[length(winEnds)] + 1 < noOfCs) {
      winEnds <- winEnds[-length(winEnds)]
    }
    if(winEnds[length(winEnds)] != nrow(y)) {
      winEnds <- c(winEnds, nrow(y))
    }
  } else {
    winEnds <- nrow(y)
  }

  methDat <- NULL
  for(i in 1:length(winStarts)) {
    chunk <- y[winStarts[i] : winEnds[i],]
    # Proceed only if chunk contains rows (cytosine positions)
    if(dim(chunk)[1] > 0) {
      midpoint <- ( min(chunk[,2]) + max(chunk[,2]) ) / 2
      per_chunk_methyl_mean <- mean(chunk[,9])
      methDat_mean <- data.frame(chr = chunk[,1][1],
                                 midpoint = midpoint,
                                 per_chunk_methyl_mean = per_chunk_methyl_mean,
                                 start = min(chunk[,2]),
                                 end = max(chunk[,2]),
                                 read = chunk[,5][1],
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
#}, mc.cores = detectCores()*0.30, mc.preschedule = T)
}, mc.cores = detectCores()*0.30, mc.preschedule = T))

print("I'm done!")

per_read_DNAmeth_DF <- per_read_DNAmeth_DF[
                         order(per_read_DNAmeth_DF[,1], per_read_DNAmeth_DF[,2]),
                       ]

write.table(per_read_DNAmeth_DF,
            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
                          "_raw_readWinSize", noOfCs, "Cs_per_readWin_midpoint.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
