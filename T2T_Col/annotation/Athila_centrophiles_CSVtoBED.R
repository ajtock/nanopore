#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from CSV (1-based start coordinates) into BED format (0-based start coordinates)

# Usage:
# ./Athila_centrophiles_CSVtoBED.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

AthilaDir <- "CENAthila/"
gapDir <- "CENgap/"
soloLTRDir <- "CENsoloLTR/"
system(paste0("[ -d ", AthilaDir, " ] || mkdir -p ", AthilaDir))
system(paste0("[ -d ", gapDir, " ] || mkdir -p ", gapDir))
system(paste0("[ -d ", soloLTRDir, " ] || mkdir -p ", soloLTRDir))

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of centromeric gap and Athila sequence coordinates
tab <- read.csv("TEs/t2t.largegaps.csv", header = T)
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
# Get rows corresponding to gaps that each contain a paired LTR element
Athila_rows <- tab[tab$TE_class == "Athila",]
print(dim(Athila_rows))
#[1] 50 22
gap_rows <- tab[tab$TE_class %in% c("Athila", "LTR_Gypsy", "LTR_Gypsy/EnSpm"),]
print(dim(gap_rows))
#[1] 52 22
soloLTR_rows <- tab[tab$TE_class == "soloLTR",]
print(dim(soloLTR_rows))
#[1] [1] 14 22

# Convert Athila into GRanges
AthilaGR <- GRanges(seqnames = Athila_rows$chr,
                    ranges = IRanges(start = Athila_rows$ltr1_start_t2t,
                                     end = Athila_rows$ltr2_end_t2t),
                    strand = "*",
                    name = Athila_rows$name,
                    ltrpid = Athila_rows$all.ltrpid,
                    class = Athila_rows$Athila_class)
# Sort AthilaGR object (by chromosome, strand, start coordinate and, finally, end coordinate)
AthilaGR <- sort(AthilaGR)
AthilaGR <- AthilaGR[seqnames(AthilaGR) %in% chrName]
Athila_bed <- data.frame(chr = as.character(seqnames(AthilaGR)),
                         start = as.integer(start(AthilaGR)-1),
                         end = as.integer(end(AthilaGR)),
                         name = as.character(AthilaGR$name),
                         score = as.numeric(AthilaGR$ltrpid),
                         strand = as.character(strand(AthilaGR)))
write.table(Athila_bed,
            file = paste0(AthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert gap into GRanges
gapGR <- GRanges(seqnames = gap_rows$chr,
                 ranges = IRanges(start = gap_rows$start,
                                  end = gap_rows$stop),
                 strand = "*",
                 name = gap_rows$name,
                 ltrpid = gap_rows$all.ltrpid,
                 class = gap_rows$TE_class)
# Sort gapGR object (by chromosome, strand, start coordinate and, finally, end coordinate)
gapGR <- sort(gapGR)
gapGR <- gapGR[seqnames(gapGR) %in% chrName]
gap_bed <- data.frame(chr = as.character(seqnames(gapGR)),
                      start = as.integer(start(gapGR)-1),
                      end = as.integer(end(gapGR)),
                      name = as.character(gapGR$name),
                      score = as.numeric(gapGR$ltrpid),
                      strand = as.character(strand(gapGR)))
write.table(gap_bed,
            file = paste0(gapDir, "CENgap_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert soloLTR into GRanges
soloLTRGR <- GRanges(seqnames = soloLTR_rows$chr,
                     ranges = IRanges(start = soloLTR_rows$start,
                                      end = soloLTR_rows$stop),
                     strand = "*",
                     name = soloLTR_rows$name,
                     ltrpid = soloLTR_rows$all.ltrpid,
                     class = soloLTR_rows$TE_class)
# Sort soloLTRGR object (by chromosome, strand, start coordinate and, finally, end coordinate)
soloLTRGR <- sort(soloLTRGR)
soloLTRGR <- soloLTRGR[seqnames(soloLTRGR) %in% chrName]
soloLTR_bed <- data.frame(chr = as.character(seqnames(soloLTRGR)),
                          start = as.integer(start(soloLTRGR)-1),
                          end = as.integer(end(soloLTRGR)),
                          name = as.character(soloLTRGR$name),
                          score = as.numeric(soloLTRGR$ltrpid),
                          strand = as.character(strand(soloLTRGR)))
write.table(soloLTR_bed,
            file = paste0(soloLTRDir, "CENsoloLTR_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
