#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric gaps between CEN180 sequences (either flanked by CEN180 SNVs or not)
# identified in t2t-col.20210610 from CSV (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./CENgapSNVboolean_CSVtoBED_AlexBousios_v080721.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(RColorBrewer)
library(scales)
revSpectralScale11 <- rev(brewer.pal(11, "Spectral"))
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}

CENgapSNVTDir <- "CENgapSNVT/"
CENgapSNVFDir <- "CENgapSNVF/"
system(paste0("[ -d ", CENgapSNVTDir, " ] || mkdir -p ", CENgapSNVTDir))
system(paste0("[ -d ", CENgapSNVFDir, " ] || mkdir -p ", CENgapSNVFDir))

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

## Load table of centromeric gap and Athila sequence coordinates
tab <- read.csv("CENgapSNVboolean_in_t2t-col.20210610.csv",
                header = T)
colnames(tab)[1] <- "gap_name"
colnames(tab)[2] <- "chr"
tab$chr <- gsub(pattern = "^", replacement = "Chr", x = tab$chr)
colnames(tab)[3] <- "gap_start"
colnames(tab)[4] <- "gap_end"
colnames(tab)[5] <- "variant"
tab$variant <- gsub(pattern = "yes", replacement = as.logical(1), x = tab$variant)
tab$variant <- gsub(pattern = "no", replacement = as.logical(0), x = tab$variant)

CENgapSNVT <- tab[tab$variant == "TRUE",]
CENgapSNVF <- tab[tab$variant == "FALSE",]
stopifnot(nrow(tab) == nrow(CENgapSNVT) + nrow(CENgapSNVF))

# Convert CENgapSNVT into GRanges and then BED format
CENgapSNVTGR <- GRanges(seqnames = CENgapSNVT$chr,
                        ranges = IRanges(start = CENgapSNVT$gap_start,
                                         end = CENgapSNVT$gap_end),
                        strand = "*",
                        name = CENgapSNVT$gap_name)
CENgapSNVTGR <- CENgapSNVTGR[seqnames(CENgapSNVTGR) %in% chrName]
CENgapSNVT_bed <- data.frame(chr = as.character(seqnames(CENgapSNVTGR)),
                             start = as.integer(start(CENgapSNVTGR)-1),
                             end = as.integer(end(CENgapSNVTGR)),
                             name = as.character(CENgapSNVTGR$name),
                             score = ".",
                             strand = as.character(strand(CENgapSNVTGR)))
write.table(CENgapSNVT_bed,
            file = paste0(CENgapSNVTDir, "CENgapSNVT_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENgapSNVF into GRanges and then BED format
CENgapSNVFGR <- GRanges(seqnames = CENgapSNVF$chr,
                        ranges = IRanges(start = CENgapSNVF$gap_start,
                                         end = CENgapSNVF$gap_end),
                        strand = "*",
                        name = CENgapSNVF$gap_name)
CENgapSNVFGR <- CENgapSNVFGR[seqnames(CENgapSNVFGR) %in% chrName]
CENgapSNVF_bed <- data.frame(chr = as.character(seqnames(CENgapSNVFGR)),
                             start = as.integer(start(CENgapSNVFGR)-1),
                             end = as.integer(end(CENgapSNVFGR)),
                             name = as.character(CENgapSNVFGR$name),
                             score = ".",
                             strand = as.character(strand(CENgapSNVFGR)))
write.table(CENgapSNVF_bed,
            file = paste0(CENgapSNVFDir, "CENgapSNVF_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
