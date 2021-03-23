#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_TABtoBED_AlexBousios_strand_info_v230321.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
tab <- read.table("TEs/t2t_Athila_master.txt_coord.corrected_50Athila_strand.bed.clean.matrix", header = T, na.strings = "na")
colnames(tab)[1] <- "chr"
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)

# Convert Athila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
# (with the exception of 4 elements without LTR coordinates)
AthilaGR <- GRanges(seqnames = tab$chr,
                    ranges = IRanges(start = tab$genome_left_coord_FL,
                                     end = tab$genome_right_coord_FL),
                    strand = tab$direction,
                    name = tab$name,
                    class = tab$class1)
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
