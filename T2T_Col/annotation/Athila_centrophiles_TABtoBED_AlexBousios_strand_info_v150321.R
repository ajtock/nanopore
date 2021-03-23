#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_TABtoBED_AlexBousios_strand_info.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

AthilaDir <- "CENAthila/"
system(paste0("[ -d ", AthilaDir, " ] || mkdir -p ", AthilaDir))

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
tab <- read.table("TEs/t2t_Athila_master.txt_coord.corrected_50Athila_strand", header = F)
colnames(tab) <- c("catinfo", "chr", "start", "end", "strand")
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
# Get Athila name and class from "catinfo" column
tab$name <- sub("-[A-z].+", replacement = "", x = tab$catinfo)
tab$name <- sub("^Athila.*-", replacement = "", x = tab$name)
tab$class <- sub(pattern = "-.+", replacement = "", x = tab$catinfo)

# Convert Athila into GRanges
AthilaGR <- GRanges(seqnames = tab$chr,
                    ranges = IRanges(start = tab$start,
                                     end = tab$end),
                    strand = tab$strand,
                    name = tab$name,
                    class = tab$class)
# Sort AthilaGR object (by chromosome, strand, start coordinate and, finally, end coordinate)
AthilaGR <- sortSeqlevels(AthilaGR)
AthilaGR <- sort(AthilaGR)
AthilaGR <- AthilaGR[seqnames(AthilaGR) %in% chrName]
Athila_bed <- data.frame(chr = as.character(seqnames(AthilaGR)),
                         start = as.integer(start(AthilaGR)-1),
                         end = as.integer(end(AthilaGR)),
                         name = as.character(AthilaGR$name),
                         score = as.character(AthilaGR$class),
                         strand = as.character(strand(AthilaGR)))
write.table(Athila_bed,
            file = paste0(AthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
