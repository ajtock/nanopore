#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# chmod +x per_read_methylation_proportion_winByCs.R
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winByCs_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 20 CHG 0.20"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winByCs_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 20 CpG 0.20"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#noOfCs <- 20
#context <- "CpG"
#CPUpc <- 0.20

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
noOfCs <- as.integer(args[3])
context <- args[4]
CPUpc <- as.numeric(args[5])

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
 
# Read in the raw output .tsv file from Deepsignal methylation model
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                         sampleName, "_MappedOn_", refbase, "_", context, "_raw.tsv"),
                  header = F)
 
# Generate a character vector of the unique readIDs in the file
readIDs <- unique(tab$V5)

# Loop within parallelised loop to calculate the per-read methylation proportion
# across sequential adjacent noOfCs-Cs-containing windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
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
      if(nrow(y) >= noOfCs) {
        winEnds <- seq(from = noOfCs,
                       to = nrow(y),
                       by = noOfCs)
        if(winEnds[length(winEnds)] != nrow(y)) {
          winEnds <- c(winEnds, nrow(y))
        }
        # Remove the penultimate winEnd value if there are fewer than noOfCs Cs from
        # this value to the last C in the read (the last row), so that the last window
        # always has as much methylation-state information as, or more than, the other windows
        if(nrow(y) - winEnds[(length(winEnds) - 1)] < noOfCs) {
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
      methDat_mean <- data.frame(chr = readwin[,1][1],
                                 midpoint = midpoint,
                                 per_readwin_methyl_mean = per_readwin_methyl_mean,
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
                          "_raw_readBinSize", noOfCs, "Cs_per_readWin_midpoint.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
