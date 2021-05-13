#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_noncentrophiles_TABtoBED_AlexBousios_strand_info_v270421_5LTR.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

CENAthilaDir <- "CENAthila/"
nonCENAthilaDir <- "nonCENAthila/"
system(paste0("[ -d ", CENAthilaDir, " ] || mkdir -p ", CENAthilaDir))
system(paste0("[ -d ", nonCENAthilaDir, " ] || mkdir -p ", nonCENAthilaDir))

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

## Load table of centromeric gap and Athila sequence coordinates
tab <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/Athilas_fulllength_inCentr+outCentr_v270421.bin",
                  header = T, na.strings = "na")
colnames(tab)[1] <- "chr"
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
tab$phylo <- toupper(tab$phylo)

CENAthila <- tab[tab$position == "in_centr",]
nonCENAthila <- tab[tab$position == "out_centr",]

CENAthila_forward <- CENAthila[CENAthila$direction == "+",]
CENAthila_reverse <- CENAthila[CENAthila$direction == "-",]

nonCENAthila_forward <- nonCENAthila[nonCENAthila$direction == "+",]
nonCENAthila_reverse <- nonCENAthila[nonCENAthila$direction == "-",]

# Convert CENAthila into GRanges
# Use genome_5LTR_start and genome_5LTR_end as
# element boundaries so that start and end coordinates
# correspond to 5' LTR start and 5' LTR end coordinates
CENAthilaGR <- sort(c(
                      GRanges(seqnames = CENAthila_forward$chr,
                              ranges = IRanges(start = CENAthila_forward$genome_5LTR_start,
                                               end = CENAthila_forward$genome_5LTR_end),
                              strand = CENAthila_forward$direction,
                              name = CENAthila_forward$gap_name,
                              class = CENAthila_forward$phylo),
                      GRanges(seqnames = CENAthila_reverse$chr,
                              ranges = IRanges(start = CENAthila_reverse$genome_5LTR_end,
                                               end = CENAthila_reverse$genome_5LTR_start),
                              strand = CENAthila_reverse$direction,
                              name = CENAthila_reverse$gap_name,
                              class = CENAthila_reverse$phylo)
                   ), ignore.strand = TRUE)
CENAthilaGR <- CENAthilaGR[seqnames(CENAthilaGR) %in% chrName]
CENAthila_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                         start = as.integer(start(CENAthilaGR)-1),
                         end = as.integer(end(CENAthilaGR)),
                         name = as.character(CENAthilaGR$name),
                         score = as.character(CENAthilaGR$class),
                         strand = as.character(strand(CENAthilaGR)))
write.table(CENAthila_bed,
            file = paste0(CENAthilaDir, "CENAthila5LTR_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert nonCENAthila into GRanges
# Use genome_5LTR_start and genome_5LTR_end as
# element boundaries so that start and end coordinates
# correspond to 5' LTR start and 5' LTR end coordinates
nonCENAthilaGR <- sort(c(
                         GRanges(seqnames = nonCENAthila_forward$chr,
                                 ranges = IRanges(start = nonCENAthila_forward$genome_5LTR_start,
                                                  end = nonCENAthila_forward$genome_5LTR_end),
                                 strand = nonCENAthila_forward$direction,
                                 name = nonCENAthila_forward$gap_name,
                                 class = nonCENAthila_forward$phylo),
                         GRanges(seqnames = nonCENAthila_reverse$chr,
                                 ranges = IRanges(start = nonCENAthila_reverse$genome_5LTR_end,
                                                  end = nonCENAthila_reverse$genome_5LTR_start),
                                 strand = nonCENAthila_reverse$direction,
                                 name = nonCENAthila_reverse$gap_name,
                                 class = nonCENAthila_reverse$phylo)
                      ), ignore.strand = TRUE)
nonCENAthilaGR <- nonCENAthilaGR[seqnames(nonCENAthilaGR) %in% chrName]
nonCENAthila_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                               start = as.integer(start(nonCENAthilaGR)-1),
                               end = as.integer(end(nonCENAthilaGR)),
                               name = as.character(nonCENAthilaGR$name),
                               score = as.character(nonCENAthilaGR$class),
                               strand = as.character(strand(nonCENAthilaGR)))
write.table(nonCENAthila_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila5LTR_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
