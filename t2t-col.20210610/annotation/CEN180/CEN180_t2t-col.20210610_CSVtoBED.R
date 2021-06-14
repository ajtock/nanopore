#!/applications/R/R-4.0.0/bin/Rscript

# Convert CEN180 sequence coordinates identified in t2t-col.20210610
# from CSV (1-based start coordinates) into BED format (0-based start coordinates)
# Convert coordinates corresponding to the intervening sequences between CEN180
# sequences from CSV into BED format

# Usage:
# ./CEN180_t2t-col.20210610_CSVtoBED.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")

CENstart <- c(14841110,3823792,13597188,4203902,11784131)[which(fai$V1 %in% chrName)]
CENend <- c(17559778,6045243,15733925,6977949,14551809)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of CEN180 sequence coordinates
# with weighted SNVs vs consensus
tab <- read.csv("T2T0610cen180sV.1.0.csv", header = T)
tab <- tab[,-which(colnames(tab) == "strand")]
colnames(tab)[which(colnames(tab) == "strand.1")] <- "strand"
tab$strand <- gsub(pattern = "cen180plus", replacement = "+",
                   x = tab$strand,
                   ignore.case = TRUE)
tab$strand <- gsub(pattern = "cen180minus", replacement = "-",
                   x = tab$strand,
                   ignore.case = TRUE)
tab$chromosome <- gsub(pattern = "^", replacement = "Chr",
                       x = tab$chromosome)
CEN180GR <- GRanges(seqnames = tab$chromosome,
                    ranges = IRanges(start = tab$start,
                                     end = tab$end),
                    strand = tab$strand,
                    index = as.integer(rownames(tab)),
                    weightedSNV = tab$weightedSNV,
                    HORlengthsSum = tab$HORlengthsSum,
                    HORcount = tab$HORcount,
                    percentageIdentity = tab$percentageIdentity)

# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
# Sorting data.frame by multiple columns (including that correspodning to strand)
# would be more complicated
# Necessary to include sorting by strand because intervening CEN180 sequences on the
# opposite strand would otherwise result in non-detection of tandem repeats on the same strand
CEN180GR <- sort(CEN180GR)
CEN180GR <- CEN180GR[seqnames(CEN180GR) %in% chrName]
CEN180_bed <- data.frame(chr = as.character(seqnames(CEN180GR)),
                         start = start(CEN180GR)-1,
                         end = end(CEN180GR),
                         name = CEN180GR$index,
                         score = CEN180GR$weightedSNV,
                         strand = as.character(strand(CEN180GR)),
                         HORlengthsSum = CEN180GR$HORlengthsSum,
                         HORcount = CEN180GR$HORcount,
                         percentageIdentity = CEN180GR$percentageIdentity)
write.table(CEN180_bed,
            file = paste0("CEN180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as CEN180_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as CEN180GR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(CEN180ChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(CEN180ChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(CEN180ChrGR)),
                         strand = strand(CEN180ChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
stopifnot(identical(width(ranLocGR), width(CEN180GR)))
stopifnot(identical(as.character(seqnames(ranLocGR)), as.character(seqnames(CEN180GR))))
stopifnot(identical(strand(ranLocGR), strand(CEN180GR)))
ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = start(ranLocGR)-1,
                         end = end(ranLocGR),
                         name = 1:length(ranLocGR),
                         score = "NA",
                         strand = strand(ranLocGR))
write.table(ranLoc_bed,
            file = paste0("CEN180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as CEN180GR
CENranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  # Contract CENChrGR so that CENrandom loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(CENChrGR) <- end(CENChrGR)-max(width(CEN180ChrGR))-2000
  start(CENChrGR) <- start(CENChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  CENranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {           
                                                                start(CENChrGR[x]) : end(CENChrGR[x])          
                                                              })),
                                         n = length(CEN180ChrGR))
  CENranLocChrGR <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = CENranLocChrStart,
                                             width = width(CEN180ChrGR)),
                            strand = strand(CEN180ChrGR))
  CENranLocGR <- append(CENranLocGR, CENranLocChrGR)
}
stopifnot(identical(width(CENranLocGR), width(CEN180GR)))
stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CEN180GR))))
stopifnot(identical(strand(CENranLocGR), strand(CEN180GR)))
CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = start(CENranLocGR)-1,
                            end = end(CENranLocGR),
                            name = 1:length(CENranLocGR),
                            score = "NA",
                            strand = strand(CENranLocGR))
write.table(CENranLoc_bed,
            file = paste0("CEN180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define adjacent centromeric loci each with widths equal to the median width of
# CEN180 sequences in (a) specified chromosome(s)
# Extract those that do not overlap CEN180 sequences
CENgapGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  CENgapChrStart <- seq(from = start(CENChrGR),
                        to = end(CENChrGR),
                        by = median(width(CEN180ChrGR)))  
  CENgapChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = CENgapChrStart,
                                          width = median(width(CEN180ChrGR))),
                         strand = "*")
  CENgapGR <- append(CENgapGR, CENgapChrGR)
}
CENgapGR_CEN180GR_ol <- findOverlaps(query = CEN180GR,
                                     subject = CENgapGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
CENgapGR <- CENgapGR[-subjectHits(CENgapGR_CEN180GR_ol)]

# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
CENgapGR <- sort(CENgapGR)
CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
CENgap_bed <- data.frame(chr = as.character(seqnames(CENgapGR)),
                         start = start(CENgapGR)-1,
                         end = end(CENgapGR),
                         name = 1:length(CENgapGR),
                         score = "NA",
                         strand = strand(CENgapGR))
write.table(CENgap_bed,
            file = paste0("CENgap_180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
