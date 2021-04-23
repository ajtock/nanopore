#!/applications/R/R-4.0.0/bin/Rscript

# Convert non-centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_noncentrophiles_TABtoBED_AlexBousios_strand_info_v210421.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

AthilaDir <- "nonCENAthila/"
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

# Load table of non-centromeric Athila
tab <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/Athilas_nonCentromeric.bin",
                  header = T, sep = "\t")
colnames(tab)[1] <- "chr"
tab$phylo <- toupper(tab$phylo)
Athila <- tab

# Convert Athila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
AthilaGR <- GRanges(seqnames = Athila$chr,
                    ranges = IRanges(start = Athila$genome_left_coord_FL,
                                     end = Athila$genome_right_coord_FL),
                    strand = Athila$direction,
                    name = Athila$TE_ID,
                    class = Athila$phylo)
AthilaGR <- AthilaGR[seqnames(AthilaGR) %in% chrName]
Athila_bed <- data.frame(chr = as.character(seqnames(AthilaGR)),
                         start = as.integer(start(AthilaGR)-1),
                         end = as.integer(end(AthilaGR)),
                         name = as.character(AthilaGR$name),
                         score = as.character(AthilaGR$class),
                         strand = as.character(strand(AthilaGR)))
write.table(Athila_bed,
            file = paste0(AthilaDir, "nonCENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
