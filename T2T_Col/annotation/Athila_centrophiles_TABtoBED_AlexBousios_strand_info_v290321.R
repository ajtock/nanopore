#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_TABtoBED_AlexBousios_strand_info_v290321.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

AthilaDir <- "CENAthila/"
soloLTRDir <- "CENsoloLTR/" 
gapDir <- "CENgap/"
system(paste0("[ -d ", AthilaDir, " ] || mkdir -p ", AthilaDir))
system(paste0("[ -d ", soloLTRDir, " ] || mkdir -p ", soloLTRDir))
system(paste0("[ -d ", gapDir, " ] || mkdir -p ", gapDir))

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
tab <- read.csv("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/t2t.largegaps--Athilafinal.csv",
                header = T, na.strings = "na")
colnames(tab)[1] <- "chr"
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
# Redefine Athila6A and Athila6B assignments based on Alex's analysis;
# Alex: "Note that you need to include in 6B some elements that currently are in 6A.
#        In the file that I gave you these are the Athila6v2 elements
#        (6 in total of which the original 6B belongs to)."
tab$class_v290321 <- tab$class2
tab$class_v290321[grep("_v2", tab$class_v290321)] <- "Athila6B"
tab$class_v290321[grep("6A", tab$class_v290321)] <- "Athila6A"
tab$class_v290321 <- toupper(tab$class_v290321)
tab$class_v290321 <- gsub("ATHILA1", "ATHILA", tab$class_v290321)

# Get row indices for elements rejected or identified as soloLTRs by Alex;
# Alex: "I have rejected 4 [now 5] elements that are a mosaic of stuff or non-TE sequence that just add noise.
#        These elements do not have coordinates for LTRs/internal domain (search for 'na').
#        They are ATGP2, 4_8 (ATHILA), and 5_37 and 5_42 (both ATHILA6A)."
reject_rowIndices <- which(is.na(tab$genome_5LTR_start))

Athila <- tab[-reject_rowIndices,]
soloLTR <- tab[tab$class == "soloLTR",]
gap <- tab

# Convert Athila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
# (with the exception of 4 elements without LTR coordinates)
AthilaGR <- GRanges(seqnames = Athila$chr,
                    ranges = IRanges(start = Athila$genome_left_coord_FL,
                                     end = Athila$genome_right_coord_FL),
                    strand = Athila$direction,
                    name = Athila$name,
                    class = Athila$class_v290321)
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

# Convert soloLTR into GRanges and then BED format
soloLTRGR <- GRanges(seqnames = soloLTR$chr,
                     ranges = IRanges(start = soloLTR$genome_left_coord_FL,
                                      end = soloLTR$genome_right_coord_FL),
                     strand = soloLTR$direction,
                     name = soloLTR$name,
                     class = soloLTR$class_v290321)
soloLTRGR <- soloLTRGR[seqnames(soloLTRGR) %in% chrName]
if(length(soloLTRGR) > 0) {
  soloLTR_bed <- data.frame(chr = as.character(seqnames(soloLTRGR)),
                            start = as.integer(start(soloLTRGR)-1),
                            end = as.integer(end(soloLTRGR)),
                            name = as.character(soloLTRGR$name),
                            score = as.character(soloLTRGR$class),
                            strand = as.character(strand(soloLTRGR)))
  write.table(soloLTR_bed,
              file = paste0(soloLTRDir, "CENsoloLTR_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert gap into GRanges and then BED format
gapGR <- GRanges(seqnames = gap$chr,
                 ranges = IRanges(start = gap$start,
                                  end = gap$stop),
                 strand = gap$direction,
                 name = gap$name,
                 class = gap$class_v290321)
gapGR <- gapGR[seqnames(gapGR) %in% chrName]
gap_bed <- data.frame(chr = as.character(seqnames(gapGR)),
                      start = as.integer(start(gapGR)-1),
                      end = as.integer(end(gapGR)),
                      name = as.character(gapGR$name),
                      score = as.character(gapGR$class),
                      strand = as.character(strand(gapGR)))
write.table(gap_bed,
            file = paste0(gapDir, "CENgap_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
