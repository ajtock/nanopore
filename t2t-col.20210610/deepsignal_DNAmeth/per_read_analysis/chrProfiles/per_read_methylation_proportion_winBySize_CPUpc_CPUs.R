#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winBySize_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CHG 0.30"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winBySize_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CpG 0.20"

# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#readBinSize <- 1000
#context <- "CpG"
#CPUpc <- 0.20

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
readBinSize <- as.integer(args[3])
context <- args[4]
CPUpc <- as.numeric(args[5])

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomeRanges)
library(dplyr)
 
if(floor(log10(readBinSize)) + 1 < 4) {
  readBinName <- paste0(readBinSize, "bp")
} else if(floor(log10(readBinSize)) + 1 >= 4 &
          floor(log10(readBinSize)) + 1 <= 6) {
  readBinName <- paste0(readBinSize/1e3, "kb")
} else if(floor(log10(readBinSize)) + 1 >= 7) {
  readBinName <- paste0(readBinSize/1e6, "Mb")
}

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai[,1][1:5]
chrLens <- fai[,2][1:5]

# Load coordinates for mitochondrial insertion on Chr2, in BED format
mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
                       header = F)
colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
mito_ins_GR <- GRanges(seqnames = "Chr2",
                       ranges = IRanges(start = min(mito_ins$start)+1,
                                        end = max(mito_ins$end)),
                       strand = "*")

# Read in the raw output .tsv file from Deepsignal methylation model
tab_list <- mclapply(seq_along(chrs), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                    sampleName, "_MappedOn_", refbase, "_", context, "_raw_", chrs[x], ".tsv"),
             header = F)
}, mc.cores = length(chrs), mc.preschedule = F)

if(length(chrs) > 1) {
  tab <- dplyr::bind_rows(tab_list)
} else {
  tab <- tab_list[[1]]
}
rm(tab_list); gc()

# Get reads that overlap mito_ins_GR
tab_mito <- tab[tab[,1] == as.character(seqnames(mito_ins_GR)) &
                tab[,2] >= start(mito_ins_GR) &
                tab[,2] <= end(mito_ins_GR),]
tab_mito_reads <- unique(tab_mito[,5])

read_within_mito_ins <- function(DSrawDF, readID, mito_ins_GR) {
  DSrawDF_read <- DSrawDF[DSrawDF[,5] == readID,]
  stopifnot(unique(DSrawDF_read[,1]) == as.character(seqnames(mito_ins_GR)))
  bool <- min(DSrawDF_read[,2], na.rm = T) >= start(mito_ins_GR) &&
          max(DSrawDF_read[,2], na.rm = T) <= end(mito_ins_GR)
  return(bool)
}

tab_mito_reads_bool <- mclapply(tab_mito_reads, function(x) {
  read_within_mito_ins(DSrawDF = tab,
                       readID = x,
                       mito_ins_GR = mito_ins_GR)
}, mc.cores = detectCores(), mc.preschedule = T)

tab_within_mito_reads <- tab_mito_reads[unlist(tab_mito_reads_bool)]

tab <- tab[!(tab[,5] %in% tab_within_mito_reads),]

# Generate a character vector of the unique readIDs in the file
readIDs <- unique(tab$V5)

# Loop within parallelised loop to calculate the per-read methylation proportion
# across sequential adjacent noOfCs-Cs-containing windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
  y <- tab[tab[,5] == x,]              # Get rows (cytosine positions) for read x
  y <- y[order(y$V2, decreasing = F),] # Order the rows by ascending position in the chromosome
 
  # Define window start and end coordinates within read
  winStarts <- seq(from = min(y[,2]),
                   to = max(y[,2]),
                   by = readBinSize)
  if(winStarts[length(winStarts)] == max(y[,2])) {
    winStarts <- winStarts[-length(winStarts)]
  }
  if(max(y[,2]) - min(y[,2]) + 1 >= readBinSize) {
    winEnds <- seq(from = winStarts[1] + readBinSize - 1,
                   to = max(y[,2]),
                   by = readBinSize)
    if(winEnds[length(winEnds)] != max(y[,2])) {
      winEnds <- c(winEnds, max(y[,2]))
    }
  } else {
    winEnds <- max(y[,2])
  }

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
      if(length(winStarts) > 1 && nrow(y) - winStarts[length(winStarts)] + 1 < noOfCs) {
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
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))

print("Done")

per_read_DNAmeth_DF <- per_read_DNAmeth_DF[
                         order(per_read_DNAmeth_DF[,1], per_read_DNAmeth_DF[,2]),
                       ]

write.table(per_read_DNAmeth_DF,
            file = paste0(sampleName, "_MappedOn_", refbase, "_", context,
                          "_raw_readBinSize", readBinName, "_per_readWin_midpoint.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
