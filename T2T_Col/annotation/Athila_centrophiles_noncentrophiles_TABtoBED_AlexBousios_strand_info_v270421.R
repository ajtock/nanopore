#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_noncentrophiles_TABtoBED_AlexBousios_strand_info_v270421.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(RColorBrewer)
revSpectralScale11 <- rev(brewer.pal(11, "Spectral"))
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"
))}

CENAthilaDir <- "CENAthila/"
nonCENAthilaDir <- "nonCENAthila/"
CENgapDir <- "CENgap/"
system(paste0("[ -d ", CENAthilaDir, " ] || mkdir -p ", CENAthilaDir))
system(paste0("[ -d ", nonCENAthilaDir, " ] || mkdir -p ", nonCENAthilaDir))
system(paste0("[ -d ", CENgapDir, " ] || mkdir -p ", CENgapDir))

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
CENgap <- tab[tab$position == "in_centr",]

# Convert CENAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
CENAthilaGR <- GRanges(seqnames = CENAthila$chr,
                       ranges = IRanges(start = CENAthila$genome_left_coord_FL,
                                        end = CENAthila$genome_right_coord_FL),
                       strand = CENAthila$direction,
                       name = CENAthila$gap_name,
                       class = CENAthila$phylo)
CENAthilaGR <- CENAthilaGR[seqnames(CENAthilaGR) %in% chrName]
CENAthila_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                         start = as.integer(start(CENAthilaGR)-1),
                         end = as.integer(end(CENAthilaGR)),
                         name = as.character(CENAthilaGR$name),
                         score = as.character(CENAthilaGR$class),
                         strand = as.character(strand(CENAthilaGR)))
write.table(CENAthila_bed,
            file = paste0(CENAthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
CENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                                     start = as.integer(start(CENAthilaGR)-1),
                                     end = as.integer(end(CENAthilaGR)),
                                     name = as.character(CENAthilaGR$name),
                                     score = as.integer(0),
                                     strand = as.character(strand(CENAthilaGR)))
write.table(CENAthila_nofamily_bed,
            file = paste0(CENAthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), "_nofamily.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
CENfams <- sort(unique(CENAthila_bed$score))
print(CENfams)
CENAthila_colofamily_bed <- data.frame(CENAthila_bed,
                                       thickStart = as.integer(0),
                                       thichEnd = as.integer(0),
                                       itemRgb = ".")
CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")

CENAthila_colofamily_bed$score <- as.integer(0)
write.table(CENAthila_colofamily_bed,
            file = paste0(CENAthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), "_colofamily.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Convert nonCENAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
nonCENAthilaGR <- GRanges(seqnames = nonCENAthila$chr,
                          ranges = IRanges(start = nonCENAthila$genome_left_coord_FL,
                                           end = nonCENAthila$genome_right_coord_FL),
                          strand = nonCENAthila$direction,
                          name = nonCENAthila$TE_ID,
                          class = nonCENAthila$phylo)
nonCENAthilaGR <- nonCENAthilaGR[seqnames(nonCENAthilaGR) %in% chrName]
nonCENAthila_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                               start = as.integer(start(nonCENAthilaGR)-1),
                               end = as.integer(end(nonCENAthilaGR)),
                               name = as.character(nonCENAthilaGR$name),
                               score = as.character(nonCENAthilaGR$class),
                               strand = as.character(strand(nonCENAthilaGR)))
write.table(nonCENAthila_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
nonCENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                                        start = as.integer(start(nonCENAthilaGR)-1),
                                        end = as.integer(end(nonCENAthilaGR)),
                                        name = as.character(nonCENAthilaGR$name),
                                        score = as.integer(0),
                                        strand = as.character(strand(nonCENAthilaGR)))
write.table(nonCENAthila_nofamily_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), "_nofamily.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
nonCENfams <- sort(unique(nonCENAthila_bed$score))
print(nonCENfams)
nonCENAthila_colofamily_bed <- data.frame(nonCENAthila_bed,
                                          thickStart = as.integer(0),
                                          thichEnd = as.integer(0),
                                          itemRgb = ".")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "nonCENATHILA0_I",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[10])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[7])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[5])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")
nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA7A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[1])), collapse = ",")

nonCENAthila_colofamily_bed$score <- as.integer(0)
write.table(nonCENAthila_colofamily_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), "_colofamily.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENgap into GRanges and then BED format
CENgapGR <- GRanges(seqnames = CENgap$chr,
                    ranges = IRanges(start = CENgap$gap_start,
                                     end = CENgap$gap_stop),
                    strand = CENgap$direction,
                    name = CENgap$gap_name,
                    class = CENgap$phylo)
CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
CENgap_bed <- data.frame(chr = as.character(seqnames(CENgapGR)),
                      start = as.integer(start(CENgapGR)-1),
                      end = as.integer(end(CENgapGR)),
                      name = as.character(CENgapGR$name),
                      score = as.character(CENgapGR$class),
                      strand = as.character(strand(CENgapGR)))
write.table(CENgap_bed,
            file = paste0(CENgapDir, "CENgap_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
