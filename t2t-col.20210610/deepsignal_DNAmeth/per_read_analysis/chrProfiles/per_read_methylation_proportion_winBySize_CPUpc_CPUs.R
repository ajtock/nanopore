#!/applications/R/R-4.0.0/bin/Rscript

# Usage on hydrogen node7:
# csmit -m 300G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winBySize_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CHG 0.30 Chr2"
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript per_read_methylation_proportion_winBySize_CPUpc_CPUs.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 CpG 0.20 Chr2"

# Divide each read into adjacent segments each consisting of a given physical size,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#readBinSize <- 1000
#context <- "CpG"
#CPUpc <- 0.20
#chrName <- unlist(strsplit("Chr2", split = ","))

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
readBinSize <- as.integer(args[3])
context <- args[4]
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
library(dplyr)
 
if(floor(log10(readBinSize)) + 1 < 4) {
  readBinName <- paste0(readBinSize, "bp")
} else if(floor(log10(readBinSize)) + 1 >= 4 &
          floor(log10(readBinSize)) + 1 <= 6) {
  readBinName <- paste0(readBinSize/1e3, "kb")
} else if(floor(log10(readBinSize)) + 1 >= 7) {
  readBinName <- paste0(readBinSize/1e6, "Mb")
}

outDir <- paste0("winBySize/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

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

# Identify and remove reads whose alignment start and end coordinates are both
# contained wholly within the boundaries of the mitochondrial insertion on Chr2,
# as we cannot be sure that these reads come from the nuclear genome
if(length(chrName) == 1 && chrName == "Chr2") {

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
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

  tab_within_mito_reads <- tab_mito_reads[unlist(tab_mito_reads_bool)]

  tab <- tab[!(tab[,5] %in% tab_within_mito_reads),]

}


# Generate a character vector of the unique readIDs in the file
readIDs <- unique(tab[,5])

# Loop within parallelised loop to calculate the per-read methylation proportion
# across sequential adjacent noOfCs-Cs-containing windows along each read
per_read_DNAmeth_DF <- do.call(rbind, mclapply(readIDs, function(x) {
  y <- tab[tab[,5] == x,]              # Get rows (cytosine positions) for read x
  y <- y[order(y$V2, decreasing = F),] # Order the rows by ascending position in the chromosome
 
  # Define window start coordinates within read
  winStarts <- seq(from = min(y[,2]),
                   to = max(y[,2]),
                   by = readBinSize)
  # Remove the last winStart value if it is equal to the read end coordinate
  if(winStarts[length(winStarts)] == max(y[,2])) {
    winStarts <- winStarts[-length(winStarts)]
  }
  # Remove the last winStart value if there are fewer than readBinSize bases from
  # this value to the read end coordinate (the last row in y), so that the last window
  # always has as many or more bases than the other windows
  tryCatch(
    {
      if(( length(winStarts) > 1 ) && ( max(y[,2]) - winStarts[length(winStarts)] + 1 < readBinSize )) {
        winStarts <- winStarts[-length(winStarts)]
      }
    },
    error = function(cond) {
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
      if(max(y[,2]) - min(y[,2]) + 1 >= readBinSize) {
        winEnds <- seq(from = winStarts[1] + readBinSize - 1,
                       to = max(y[,2]),
                       by = readBinSize)
        if(winEnds[length(winEnds)] != max(y[,2])) {
          winEnds <- c(winEnds, max(y[,2]))
        } 
        # Remove the penultimate winEnd value if there are fewer than readBinSize bases from
        # this value to the read end coordinate (the last row in y), so that the last window
        # always has as many or more bases than the other windows
        if(max(y[,2]) - winEnds[(length(winEnds) - 1)] < readBinSize) {
          winEnds <- winEnds[-(length(winEnds) - 1)]
        }
      } else {
        winEnds <- max(y[,2])
      }
    },
    error = function(cond) {
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
    readwin <- y[y[,2] >= winStarts[i] &
                 y[,2] <= winEnds[i],]
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
            file = paste0(outDir, sampleName, "_MappedOn_", refbase, "_", context,
                          "_raw_readBinSize", readBinName, "_per_readWin_midpoint_",
                          chrName, ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
