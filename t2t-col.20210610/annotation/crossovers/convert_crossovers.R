#!/applications/R/R-4.0.0/bin/Rscript

# Create table of crossovers using racon.cos.bed (which has 1-based start coordinates)
# and generate random loci of the same number and width distribution
# Write as BED files 

# Usage:
# ./convert_crossovers.R 'Chr1,Chr2,Chr3,Chr4,Chr5' t2t-col.20210610

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#refbase <- "t2t-col.20210610"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
refbase <- args[2]

options(stringsAsFactors = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(data.table)

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")

CEN <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.centromeres"), header = T)
CEN <- CEN[which(fai$V1 %in% chrName),]
# Not sure Dan defined pericentromeres at 10-kb resolution (appear to be at 100-kb resolution)
periCEN <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.Dan_pericentromeres"), header = T)
periCEN <- periCEN[which(fai$V1 %in% chrName),]

CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CEN$start,
                                  end = CEN$end),
                 strand = "*")
CENGR <- CENGR[which(seqnames(CENGR)@values %in% chrName)]

nonCENGR <- GRanges(seqnames = rep(chrs, 2),
                    ranges = IRanges(start = c(rep(1, length(chrs)),
                                               CEN$end+1),
                                     end = c(CEN$start-1,
                                             chrLens)),
                    strand = "*")
nonCENGR <- nonCENGR[which(seqnames(nonCENGR)@values %in% chrName)]

periCENGR <- GRanges(seqnames = chrs,
                     ranges = IRanges(start = periCEN$start,
                                      end = periCEN$end),
                     strand = "*")
periCENGR <- periCENGR[which(seqnames(periCENGR)@values %in% chrName)]

armGR <- GRanges(seqnames = rep(chrs, 2),
                 ranges = IRanges(start = c(rep(1, length(chrs)),
                                            periCEN$end+1),
                                  end = c(periCEN$start-1,
                                          chrLens)),
                 strand = "*")
armGR <- armGR[which(seqnames(armGR)@values %in% chrName)]

# Load table of crossover coordinates (which has 1-based start coordinates)
crossovers <- read.table("racon.cos.bed", header = F)
crossovers <- crossovers[,1:3]
colnames(crossovers) <- c("chr", "start", "end")
crossovers <- crossovers[which(crossovers$chr %in% chrName),]
print(dim(crossovers))
#[1] 2080    6

crossoversGR <- GRanges(seqnames = crossovers$chr,
                        ranges = IRanges(start = crossovers$start,
                                         end = crossovers$end),
                        strand = "*")
crossoversGR <- sortSeqlevels(crossoversGR)
crossoversGR <- sort(crossoversGR)

crossovers_bed <- data.frame(chr = as.character(seqnames(crossoversGR)),
                             start = as.integer(start(crossoversGR)-1),
                             end = as.integer(end(crossoversGR)),
                             name = as.integer(1:length(crossoversGR)),
                             score = rep("NA", length(crossoversGR)),
                             strand = as.character(strand(crossoversGR)))
write.table(crossovers_bed,
            file = paste0("t2t-col.20210610_crossovers_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as crossoversGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as crossoversGR
chrs <- chrs[chrs %in% chrName]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  crossoversChrGR <- crossoversGR[seqnames(crossoversGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(crossoversChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(crossoversChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(crossoversChrGR)),
                         strand = strand(crossoversChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0("t2t-col.20210610_crossovers_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
