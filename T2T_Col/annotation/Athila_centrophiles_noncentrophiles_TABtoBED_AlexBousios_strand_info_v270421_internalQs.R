#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_noncentrophiles_TABtoBED_AlexBousios_strand_info_v270421_internalQs.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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

# Define genomic coordinates of internal quarters (between genome_5LTR_end and genome_3LTR_start)
CENAthila_forward <- data.frame(CENAthila_forward,
                                length_internal = ((CENAthila_forward$genome_3LTR_start - 1) -
                                                   (CENAthila_forward$genome_5LTR_end + 1) + 1))
CENAthila_forward <- data.frame(CENAthila_forward,
                                genome_internalQ1_start = (CENAthila_forward$genome_5LTR_end + 1),
                                genome_internalQ1_end = ((CENAthila_forward$genome_5LTR_end + 1) +
                                                         (floor(CENAthila_forward$length_internal / 4) - 1)))
CENAthila_forward <- data.frame(CENAthila_forward,
                                genome_internalQ2_start = (CENAthila_forward$genome_internalQ1_end + 1),
                                genome_internalQ2_end = ((CENAthila_forward$genome_internalQ1_end + 1) +
                                                         (floor(CENAthila_forward$length_internal / 4) - 1)))
CENAthila_forward <- data.frame(CENAthila_forward,
                                genome_internalQ3_start = (CENAthila_forward$genome_internalQ2_end + 1),
                                genome_internalQ3_end = ((CENAthila_forward$genome_internalQ2_end + 1) +
                                                         (floor(CENAthila_forward$length_internal / 4) - 1)))
CENAthila_forward <- data.frame(CENAthila_forward,
                                genome_internalQ4_start = (CENAthila_forward$genome_internalQ3_end + 1),
                                genome_internalQ4_end = (CENAthila_forward$genome_3LTR_start - 1))

# Define genomic coordinates of internal quarters (between genome_3LTR_start and genome_5LTR_end)
CENAthila_reverse <- data.frame(CENAthila_reverse,
                                length_internal = ((CENAthila_reverse$genome_5LTR_end - 1) -
                                                   (CENAthila_reverse$genome_3LTR_start + 1) + 1))
CENAthila_reverse <- data.frame(CENAthila_reverse,
                                genome_internalQ1_start = ((CENAthila_reverse$genome_5LTR_end - 1) -
                                                           (floor(CENAthila_reverse$length_internal / 4) - 1)),
                                genome_internalQ1_end = (CENAthila_reverse$genome_5LTR_end - 1))
CENAthila_reverse <- data.frame(CENAthila_reverse,
                                genome_internalQ2_start = ((CENAthila_reverse$genome_internalQ1_start - 1) -
                                                           (floor(CENAthila_reverse$length_internal / 4) - 1)),
                                genome_internalQ2_end = (CENAthila_reverse$genome_internalQ1_start - 1))
CENAthila_reverse <- data.frame(CENAthila_reverse,
                                genome_internalQ3_start = ((CENAthila_reverse$genome_internalQ2_start - 1) -
                                                           (floor(CENAthila_reverse$length_internal / 4) - 1)),
                                genome_internalQ3_end = (CENAthila_reverse$genome_internalQ2_start - 1))
CENAthila_reverse <- data.frame(CENAthila_reverse,
                                genome_internalQ4_start = (CENAthila_reverse$genome_3LTR_start + 1),
                                genome_internalQ4_end = CENAthila_reverse$genome_internalQ3_start - 1)

CENAthila <- rbind(CENAthila_forward, CENAthila_reverse)
CENAthila <- CENAthila[with(CENAthila, order(chr, genome_left_coord_FL, genome_right_coord_FL)),]

# Define genomic coordinates of internal quarters (between genome_5LTR_end and genome_3LTR_start)
nonCENAthila_forward <- data.frame(nonCENAthila_forward,
                                   length_internal = ((nonCENAthila_forward$genome_3LTR_start - 1) -
                                                      (nonCENAthila_forward$genome_5LTR_end + 1) + 1))
nonCENAthila_forward <- data.frame(nonCENAthila_forward,
                                   genome_internalQ1_start = (nonCENAthila_forward$genome_5LTR_end + 1),
                                   genome_internalQ1_end = ((nonCENAthila_forward$genome_5LTR_end + 1) +
                                                            (floor(nonCENAthila_forward$length_internal / 4) - 1)))
nonCENAthila_forward <- data.frame(nonCENAthila_forward,
                                   genome_internalQ2_start = (nonCENAthila_forward$genome_internalQ1_end + 1),
                                   genome_internalQ2_end = ((nonCENAthila_forward$genome_internalQ1_end + 1) +
                                                            (floor(nonCENAthila_forward$length_internal / 4) - 1)))
nonCENAthila_forward <- data.frame(nonCENAthila_forward,
                                   genome_internalQ3_start = (nonCENAthila_forward$genome_internalQ2_end + 1),
                                   genome_internalQ3_end = ((nonCENAthila_forward$genome_internalQ2_end + 1) +
                                                            (floor(nonCENAthila_forward$length_internal / 4) - 1)))
nonCENAthila_forward <- data.frame(nonCENAthila_forward,
                                   genome_internalQ4_start = (nonCENAthila_forward$genome_internalQ3_end + 1),
                                   genome_internalQ4_end = (nonCENAthila_forward$genome_3LTR_start - 1))
 
# Define genomic coordinates of internal quarters (between genome_3LTR_start and genome_5LTR_end)
nonCENAthila_reverse <- data.frame(nonCENAthila_reverse,
                                   length_internal = ((nonCENAthila_reverse$genome_5LTR_end - 1) -
                                                      (nonCENAthila_reverse$genome_3LTR_start + 1) + 1))
nonCENAthila_reverse <- data.frame(nonCENAthila_reverse,
                                   genome_internalQ1_start = ((nonCENAthila_reverse$genome_5LTR_end - 1) -
                                                              (floor(nonCENAthila_reverse$length_internal / 4) - 1)),
                                   genome_internalQ1_end = (nonCENAthila_reverse$genome_5LTR_end - 1))
nonCENAthila_reverse <- data.frame(nonCENAthila_reverse,
                                   genome_internalQ2_start = ((nonCENAthila_reverse$genome_internalQ1_start - 1) -
                                                              (floor(nonCENAthila_reverse$length_internal / 4) - 1)),
                                   genome_internalQ2_end = (nonCENAthila_reverse$genome_internalQ1_start - 1))
nonCENAthila_reverse <- data.frame(nonCENAthila_reverse,
                                   genome_internalQ3_start = ((nonCENAthila_reverse$genome_internalQ2_start - 1) -
                                                              (floor(nonCENAthila_reverse$length_internal / 4) - 1)),
                                   genome_internalQ3_end = (nonCENAthila_reverse$genome_internalQ2_start - 1))
nonCENAthila_reverse <- data.frame(nonCENAthila_reverse,
                                   genome_internalQ4_start = (nonCENAthila_reverse$genome_3LTR_start + 1),
                                   genome_internalQ4_end = nonCENAthila_reverse$genome_internalQ3_start - 1)

nonCENAthila <- rbind(nonCENAthila_forward, nonCENAthila_reverse)
nonCENAthila <- nonCENAthila[with(nonCENAthila, order(chr, genome_left_coord_FL, genome_right_coord_FL)),]

Athila <- rbind(CENAthila, nonCENAthila)
write.table(Athila,
            file = "/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/Athilas_fulllength_inCentr+outCentr_v270421_with_internalQs.tsv",
            sep = "\t", quote = F, row.names = F) 

# Convert CENAthilaIQ1 into GRanges
CENAthilaIQ1_GR <- GRanges(seqnames = CENAthila$chr,
                           ranges = IRanges(start = CENAthila$genome_internalQ1_start,
                                            end = CENAthila$genome_internalQ1_end),
                           strand = CENAthila$direction,
                           name = CENAthila$gap_name,
                           class = CENAthila$phylo)
CENAthilaIQ1_GR <- CENAthilaIQ1_GR[seqnames(CENAthilaIQ1_GR) %in% chrName]
CENAthilaIQ1_bed <- data.frame(chr = as.character(seqnames(CENAthilaIQ1_GR)),
                               start = as.integer(start(CENAthilaIQ1_GR)-1),
                               end = as.integer(end(CENAthilaIQ1_GR)),
                               name = as.character(CENAthilaIQ1_GR$name),
                               score = as.character(CENAthilaIQ1_GR$class),
                               strand = as.character(strand(CENAthilaIQ1_GR)))
write.table(CENAthilaIQ1_bed,
            file = paste0(CENAthilaDir, "CENAthilaIQ1_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENAthilaIQ2 into GRanges
CENAthilaIQ2_GR <- GRanges(seqnames = CENAthila$chr,
                           ranges = IRanges(start = CENAthila$genome_internalQ2_start,
                                            end = CENAthila$genome_internalQ2_end),
                           strand = CENAthila$direction,
                           name = CENAthila$gap_name,
                           class = CENAthila$phylo)
CENAthilaIQ2_GR <- CENAthilaIQ2_GR[seqnames(CENAthilaIQ2_GR) %in% chrName]
CENAthilaIQ2_bed <- data.frame(chr = as.character(seqnames(CENAthilaIQ2_GR)),
                               start = as.integer(start(CENAthilaIQ2_GR)-1),
                               end = as.integer(end(CENAthilaIQ2_GR)),
                               name = as.character(CENAthilaIQ2_GR$name),
                               score = as.character(CENAthilaIQ2_GR$class),
                               strand = as.character(strand(CENAthilaIQ2_GR)))
write.table(CENAthilaIQ2_bed,
            file = paste0(CENAthilaDir, "CENAthilaIQ2_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENAthilaIQ3 into GRanges
CENAthilaIQ3_GR <- GRanges(seqnames = CENAthila$chr,
                           ranges = IRanges(start = CENAthila$genome_internalQ3_start,
                                            end = CENAthila$genome_internalQ3_end),
                           strand = CENAthila$direction,
                           name = CENAthila$gap_name,
                           class = CENAthila$phylo)
CENAthilaIQ3_GR <- CENAthilaIQ3_GR[seqnames(CENAthilaIQ3_GR) %in% chrName]
CENAthilaIQ3_bed <- data.frame(chr = as.character(seqnames(CENAthilaIQ3_GR)),
                               start = as.integer(start(CENAthilaIQ3_GR)-1),
                               end = as.integer(end(CENAthilaIQ3_GR)),
                               name = as.character(CENAthilaIQ3_GR$name),
                               score = as.character(CENAthilaIQ3_GR$class),
                               strand = as.character(strand(CENAthilaIQ3_GR)))
write.table(CENAthilaIQ3_bed,
            file = paste0(CENAthilaDir, "CENAthilaIQ3_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENAthilaIQ4 into GRanges
CENAthilaIQ4_GR <- GRanges(seqnames = CENAthila$chr,
                           ranges = IRanges(start = CENAthila$genome_internalQ4_start,
                                            end = CENAthila$genome_internalQ4_end),
                           strand = CENAthila$direction,
                           name = CENAthila$gap_name,
                           class = CENAthila$phylo)
CENAthilaIQ4_GR <- CENAthilaIQ4_GR[seqnames(CENAthilaIQ4_GR) %in% chrName]
CENAthilaIQ4_bed <- data.frame(chr = as.character(seqnames(CENAthilaIQ4_GR)),
                               start = as.integer(start(CENAthilaIQ4_GR)-1),
                               end = as.integer(end(CENAthilaIQ4_GR)),
                               name = as.character(CENAthilaIQ4_GR$name),
                               score = as.character(CENAthilaIQ4_GR$class),
                               strand = as.character(strand(CENAthilaIQ4_GR)))
write.table(CENAthilaIQ4_bed,
            file = paste0(CENAthilaDir, "CENAthilaIQ4_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Convert nonCENAthilaIQ1 into GRanges
nonCENAthilaIQ1_GR <- GRanges(seqnames = nonCENAthila$chr,
                              ranges = IRanges(start = nonCENAthila$genome_internalQ1_start,
                                               end = nonCENAthila$genome_internalQ1_end),
                              strand = nonCENAthila$direction,
                              name = nonCENAthila$gap_name,
                              class = nonCENAthila$phylo)
nonCENAthilaIQ1_GR <- nonCENAthilaIQ1_GR[seqnames(nonCENAthilaIQ1_GR) %in% chrName]
nonCENAthilaIQ1_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaIQ1_GR)),
                                  start = as.integer(start(nonCENAthilaIQ1_GR)-1),
                                  end = as.integer(end(nonCENAthilaIQ1_GR)),
                                  name = as.character(nonCENAthilaIQ1_GR$name),
                                  score = as.character(nonCENAthilaIQ1_GR$class),
                                  strand = as.character(strand(nonCENAthilaIQ1_GR)))
write.table(nonCENAthilaIQ1_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthilaIQ1_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert nonCENAthilaIQ2 into GRanges
nonCENAthilaIQ2_GR <- GRanges(seqnames = nonCENAthila$chr,
                              ranges = IRanges(start = nonCENAthila$genome_internalQ2_start,
                                               end = nonCENAthila$genome_internalQ2_end),
                              strand = nonCENAthila$direction,
                              name = nonCENAthila$gap_name,
                              class = nonCENAthila$phylo)
nonCENAthilaIQ2_GR <- nonCENAthilaIQ2_GR[seqnames(nonCENAthilaIQ2_GR) %in% chrName]
nonCENAthilaIQ2_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaIQ2_GR)),
                                  start = as.integer(start(nonCENAthilaIQ2_GR)-1),
                                  end = as.integer(end(nonCENAthilaIQ2_GR)),
                                  name = as.character(nonCENAthilaIQ2_GR$name),
                                  score = as.character(nonCENAthilaIQ2_GR$class),
                                  strand = as.character(strand(nonCENAthilaIQ2_GR)))
write.table(nonCENAthilaIQ2_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthilaIQ2_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert nonCENAthilaIQ3 into GRanges
nonCENAthilaIQ3_GR <- GRanges(seqnames = nonCENAthila$chr,
                              ranges = IRanges(start = nonCENAthila$genome_internalQ3_start,
                                               end = nonCENAthila$genome_internalQ3_end),
                              strand = nonCENAthila$direction,
                              name = nonCENAthila$gap_name,
                              class = nonCENAthila$phylo)
nonCENAthilaIQ3_GR <- nonCENAthilaIQ3_GR[seqnames(nonCENAthilaIQ3_GR) %in% chrName]
nonCENAthilaIQ3_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaIQ3_GR)),
                                  start = as.integer(start(nonCENAthilaIQ3_GR)-1),
                                  end = as.integer(end(nonCENAthilaIQ3_GR)),
                                  name = as.character(nonCENAthilaIQ3_GR$name),
                                  score = as.character(nonCENAthilaIQ3_GR$class),
                                  strand = as.character(strand(nonCENAthilaIQ3_GR)))
write.table(nonCENAthilaIQ3_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthilaIQ3_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert nonCENAthilaIQ4 into GRanges
nonCENAthilaIQ4_GR <- GRanges(seqnames = nonCENAthila$chr,
                              ranges = IRanges(start = nonCENAthila$genome_internalQ4_start,
                                               end = nonCENAthila$genome_internalQ4_end),
                              strand = nonCENAthila$direction,
                              name = nonCENAthila$gap_name,
                              class = nonCENAthila$phylo)
nonCENAthilaIQ4_GR <- nonCENAthilaIQ4_GR[seqnames(nonCENAthilaIQ4_GR) %in% chrName]
nonCENAthilaIQ4_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaIQ4_GR)),
                                  start = as.integer(start(nonCENAthilaIQ4_GR)-1),
                                  end = as.integer(end(nonCENAthilaIQ4_GR)),
                                  name = as.character(nonCENAthilaIQ4_GR$name),
                                  score = as.character(nonCENAthilaIQ4_GR$class),
                                  strand = as.character(strand(nonCENAthilaIQ4_GR)))
write.table(nonCENAthilaIQ4_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthilaIQ4_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

