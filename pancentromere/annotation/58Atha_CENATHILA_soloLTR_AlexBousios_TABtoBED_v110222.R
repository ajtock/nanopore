#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./58Atha_CENATHILA_soloLTR_AlexBousios_TABtoBED_v110222.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(scales)
# From https://github.com/dmarcelinobr/ggdecor/blob/master/R/color.R :
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
scale_colour_rich10 <- function(...) { discrete_scale("colour", "rich10", rich10(), ...) }
scale_fill_rich10 <- function(...) { discrete_scale("fill", "rich10", rich10(), ...) }

rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
scale_colour_rich12 <- function(...) { discrete_scale("colour", "rich12", rich12(), ...) }
scale_fill_rich12 <- function(...) { discrete_scale("fill", "rich12", rich12(), ...) }

CENATHILA_Dir <- "CENATHILA/"
CENsoloLTR_Dir <- "CENsoloLTR/"
system(paste0("[ -d ", CENATHILA_Dir, " ] || mkdir -p ", CENATHILA_Dir))
system(paste0("[ -d ", CENsoloLTR_Dir, " ] || mkdir -p ", CENsoloLTR_Dir))

# Load table of ATHILA sequence coordinates
tab <- read.table("58Atha_CENATHILA_soloLTR_AlexBousios_v110222.tsv",
                  header = T)
print(dim(tab))
tab <- tab[tab$chromosome %in% chrName,]
print(dim(tab))
tab$family_FL[which(tab$family_FL == "ATHILA2_check")] <- "ATHILA2"
CENATHILA <- tab[tab$quality == "intact",]
CENsoloLTR <- tab[tab$quality == "solo",]
stopifnot(nrow(CENATHILA) + nrow(CENsoloLTR) == nrow(tab))


tab$V1 <- gsub(pattern = "_RagTag_RagTag", replacement = ".RagTag.RagTag", x = tab$V1)
species <- gsub("_.+", "", tab$V1)
accession <- gsub("Atha_", "", tab$V1)
accession <- gsub("_.+", "", accession)
chromosome <- gsub("\\.\\d.+", "", tab$V1)
chromosome <- gsub(".+_", "", chromosome)
chromosome <- gsub(".RagTag.RagTag", "_RagTag_RagTag", chromosome)
coord_start <- gsub(".RagTag.RagTag", "_RagTag_RagTag", tab$V1)
coord_start <- gsub(".+\\.", "", coord_start)
coord_start <- as.integer(gsub("-.+", "", coord_start))
coord_end <- gsub(".+-", "", tab$V1)
coord_end <- as.integer(gsub("_.+", "", coord_end))
DP_strand <- gsub(".+_([DP])_.+", "\\1", tab$V1)
strand <- DP_strand
strand[which(strand == "D")] <- "+"
strand[which(strand == "P")] <- "-"
element <- gsub(".+_([DP])_", "", tab$V1)
element <- toupper(element)
element_phylo <- toupper(tab$V2)

CENATHILA <- data.frame(species = species,
                        accession = accession,
                        chr = chromosome,
                        start = coord_start,
                        end = coord_end,
                        DP_strand = DP_strand,
                        strand = strand,
                        element = element,
                        phylo = element_phylo)

CENATHILA <- CENATHILA[
               with( CENATHILA, order(species, accession, chr, start, end) ),
             ]

# Convert into GRanges
#CENATHILA_GR <- GRanges(seqnames = CENATHILA$chr,
#                        ranges = IRanges(start = CENATHILA$start, end = CENATHILA$end),
#                        strand = CENATHILA$strand,
#                        name = paste0(CENATHILA$accession, "_", CENATHILA$chr, "_",
#                                      CENATHILA$start, "_", CENATHILA$end, "_", CENATHILA$DP_strand),
#                        phylo = CENATHILA$phylo)
#CENATHILA_GR <- unique(CENATHILA_GR)
#CENATHILA_GR <- sortSeqlevels(CENATHILA_GR)
#CENATHILA_GR_sort <- sort(CENATHILA_GR, ignore.strand = T)

CENATHILA_BED <- data.frame(chr = CENATHILA$chr,
                            start = CENATHILA$start-1,
                            end = CENATHILA$end,
                            name = paste0(CENATHILA$accession, "_", CENATHILA$chr, "_",
                                          CENATHILA$start, "_", CENATHILA$end, "_", CENATHILA$DP_strand),
                            score = CENATHILA$phylo,
                            strand = CENATHILA$strand)

outDir <- "58Atha/"
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
accessionDir <- paste0(accession, "/")
for(i in accessionDir) {
  system(paste0("[ -d ", i, " ] || mkdir -p ", i))
}

write.table(CENATHILA_BED,
            file = paste0(outDir, "CENATHILA_in_58Atha_v110222.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
  nonCENfams <- sort(unique(nonCENAthila_bed$score))
  print(nonCENfams)
  CENATHILA_BED_colophylo <- data.frame(CENATHILA_BED,
                                            thickStart = as.integer(0),
                                            thickEnd = as.integer(0),
                                            itemRgb = ".")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA0",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[12])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[11])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[10])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[9])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[8])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA4A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[7])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA4C",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[6])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[5])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[4])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[3])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA7A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[2])), collapse = ",")
  CENATHILA_BED_colophylo[CENATHILA_BED_colophylo$score == "ATHILA8A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[1])), collapse = ",")
  
  nonCENAthila_colofamily_bed$score <- as.integer(0)
  write.table(nonCENAthila_colofamily_bed,
              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}


if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                                       start = as.integer(start(CENAthilaGR)-1),
                                       end = as.integer(end(CENAthilaGR)),
                                       name = as.character(CENAthilaGR$name),
                                       score = as.integer(0),
                                       strand = as.character(strand(CENAthilaGR)))
  write.table(CENAthila_nofamily_bed,
              file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENfams <- sort(unique(CENAthila_bed$score))
  print(CENfams)
  CENAthila_colofamily_bed <- data.frame(CENAthila_bed,
                                         thickStart = as.integer(0),
                                         thickEnd = as.integer(0),
                                         itemRgb = ".")
  CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
  CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
  CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
  CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
  CENAthila_colofamily_bed[CENAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")
  
  CENAthila_colofamily_bed$score <- as.integer(0)
  write.table(CENAthila_colofamily_bed,
              file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert nonCENAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
nonCENAthilaGR <- GRanges(seqnames = nonCENAthila$chr,
                          ranges = IRanges(start = nonCENAthila$genome_left_coord_FL,
                                           end = nonCENAthila$genome_right_coord_FL),
                          strand = nonCENAthila$direction,
                          name = nonCENAthila$TE_ID,
                          phylo = nonCENAthila$phylo)
nonCENAthilaGR <- unique(nonCENAthilaGR)
nonCENAthilaGR <- nonCENAthilaGR[seqnames(nonCENAthilaGR) %in% chrName]
nonCENAthilaGR <- sortSeqlevels(nonCENAthilaGR)
nonCENAthilaGR <- sort(nonCENAthilaGR, ignore.strand = TRUE)
nonCENAthila_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                               start = as.integer(start(nonCENAthilaGR)-1),
                               end = as.integer(end(nonCENAthilaGR)),
                               name = as.character(nonCENAthilaGR$name),
                               score = as.character(nonCENAthilaGR$phylo),
                               strand = as.character(strand(nonCENAthilaGR)))
write.table(nonCENAthila_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  nonCENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                                          start = as.integer(start(nonCENAthilaGR)-1),
                                          end = as.integer(end(nonCENAthilaGR)),
                                          name = as.character(nonCENAthilaGR$name),
                                          score = as.integer(0),
                                          strand = as.character(strand(nonCENAthilaGR)))
  write.table(nonCENAthila_nofamily_bed,
              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
  nonCENfams <- sort(unique(nonCENAthila_bed$score))
  print(nonCENfams)
  nonCENAthila_colofamily_bed <- data.frame(nonCENAthila_bed,
                                            thickStart = as.integer(0),
                                            thickEnd = as.integer(0),
                                            itemRgb = ".")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA0",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[10])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[7])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA7A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[1])), collapse = ",")
  
  nonCENAthila_colofamily_bed$score <- as.integer(0)
  write.table(nonCENAthila_colofamily_bed,
              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert CENsoloLTR into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries
CENsoloLTRGR <- GRanges(seqnames = CENsoloLTR$chr,
                        ranges = IRanges(start = CENsoloLTR$genome_left_coord_FL,
                                         end = CENsoloLTR$genome_right_coord_FL),
                        strand = CENsoloLTR$direction,
                        name = CENsoloLTR$TE_ID,
                        phylo = CENsoloLTR$phylo)
CENsoloLTRGR <- unique(CENsoloLTRGR)
CENsoloLTRGR <- CENsoloLTRGR[seqnames(CENsoloLTRGR) %in% chrName]
CENsoloLTRGR <- sortSeqlevels(CENsoloLTRGR)
CENsoloLTRGR <- sort(CENsoloLTRGR, ignore.strand = TRUE)
CENsoloLTR_bed <- data.frame(chr = as.character(seqnames(CENsoloLTRGR)),
                             start = as.integer(start(CENsoloLTRGR)-1),
                             end = as.integer(end(CENsoloLTRGR)),
                             name = as.character(CENsoloLTRGR$name),
                             score = as.character(CENsoloLTRGR$phylo),
                             strand = as.character(strand(CENsoloLTRGR)))
write.table(CENsoloLTR_bed,
            file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENsoloLTR_nofamily_bed <- data.frame(chr = as.character(seqnames(CENsoloLTRGR)),
                                        start = as.integer(start(CENsoloLTRGR)-1),
                                        end = as.integer(end(CENsoloLTRGR)),
                                        name = as.character(CENsoloLTRGR$name),
                                        score = as.integer(0),
                                        strand = as.character(strand(CENsoloLTRGR)))
  write.table(CENsoloLTR_nofamily_bed,
              file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENfams <- sort(unique(CENsoloLTR_bed$score))
  print(CENfams)
  CENsoloLTR_colofamily_bed <- data.frame(CENsoloLTR_bed,
                                          thickStart = as.integer(0),
                                          thickEnd = as.integer(0),
                                          itemRgb = ".")
  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA4-4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA6",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
  
  CENsoloLTR_colofamily_bed$score <- as.integer(0)
  write.table(CENsoloLTR_colofamily_bed,
              file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert CENfragmentAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries
CENfragmentAthilaGR <- GRanges(seqnames = CENfragmentAthila$chr,
                               ranges = IRanges(start = CENfragmentAthila$genome_left_coord_FL,
                                                end = CENfragmentAthila$genome_right_coord_FL),
                               strand = CENfragmentAthila$direction,
                               name = CENfragmentAthila$TE_ID,
                               phylo = CENfragmentAthila$phylo)
CENfragmentAthilaGR <- unique(CENfragmentAthilaGR)
CENfragmentAthilaGR <- CENfragmentAthilaGR[seqnames(CENfragmentAthilaGR) %in% chrName]
CENfragmentAthilaGR <- sortSeqlevels(CENfragmentAthilaGR)
CENfragmentAthilaGR <- sort(CENfragmentAthilaGR, ignore.strand = TRUE)
CENfragmentAthila_bed <- data.frame(chr = as.character(seqnames(CENfragmentAthilaGR)),
                                    start = as.integer(start(CENfragmentAthilaGR)-1),
                                    end = as.integer(end(CENfragmentAthilaGR)),
                                    name = as.character(CENfragmentAthilaGR$name),
                                    score = as.character(CENfragmentAthilaGR$phylo),
                                    strand = as.character(strand(CENfragmentAthilaGR)))
write.table(CENfragmentAthila_bed,
            file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENfragmentAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(CENfragmentAthilaGR)),
                                               start = as.integer(start(CENfragmentAthilaGR)-1),
                                               end = as.integer(end(CENfragmentAthilaGR)),
                                               name = as.character(CENfragmentAthilaGR$name),
                                               score = as.integer(0),
                                               strand = as.character(strand(CENfragmentAthilaGR)))
  write.table(CENfragmentAthila_nofamily_bed,
              file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENfams <- sort(unique(CENfragmentAthila_bed$score))
  print(CENfams)
  CENfragmentAthila_colofamily_bed <- data.frame(CENfragmentAthila_bed,
                                                 thickStart = as.integer(0),
                                                 thickEnd = as.integer(0),
                                                 itemRgb = ".")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[7])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA4-4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")
  
  CENfragmentAthila_colofamily_bed$score <- as.integer(0)
  write.table(CENfragmentAthila_colofamily_bed,
              file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert CENgap into GRanges and then BED format
CENgapGR <- GRanges(seqnames = CENgap$chr,
                    ranges = IRanges(start = as.integer(CENgap$gap_start),
                                     end = as.integer(CENgap$gap_stop)),
                    strand = CENgap$direction,
                    name = CENgap$gap_name,
                    phylo = CENgap$phylo)
CENgapGR <- unique(CENgapGR)
CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
CENgapGR <- sortSeqlevels(CENgapGR)
CENgapGR <- sort(CENgapGR, ignore.strand = TRUE)
CENgap_bed <- data.frame(chr = as.character(seqnames(CENgapGR)),
                         start = as.integer(start(CENgapGR)-1),
                         end = as.integer(end(CENgapGR)),
                         name = as.character(CENgapGR$name),
                         score = as.character(CENgapGR$phylo),
                         strand = as.character(strand(CENgapGR)))
write.table(CENgap_bed,
            file = paste0(CENgapDir, "CENgap_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

## Convert CENgapAll into GRanges and then BED format
#CENgapAllGR <- GRanges(seqnames = CENgapAll$chr,
#                       ranges = IRanges(start = as.integer(CENgapAll$gap_start),
#                                        end = as.integer(CENgapAll$gap_stop)),
#                       strand = CENgapAll$direction,
#                       name = CENgapAll$gap_name
#                       phylo = CENgapAll$phylo)
#CENgapAllGR <- unique(CENgapAllGR)
#CENgapAllGR <- CENgapAllGR[seqnames(CENgapAllGR) %in% chrName]
#CENgapAllGR <- sortSeqlevels(CENgapAllGR)
#CENgapAllGR <- sort(CENgapAllGR, ignore.strand = TRUE)
#CENgapAll_bed <- data.frame(chr = as.character(seqnames(CENgapAllGR)),
#                            start = as.integer(start(CENgapAllGR)-1),
#                            end = as.integer(end(CENgapAllGR)),
#                            name = as.character(CENgapAllGR$name),
#                            score = as.character(CENgapAllGR$phylo),
#                            strand = as.character(strand(CENgapAllGR)))
#write.table(CENgapAll_bed,
#            file = paste0(CENgapAllDir, "CENgapAll_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENgapAllAthila into GRanges and then BED format
CENgapAllAthilaGR <- GRanges(seqnames = CENgapAllAthila$chr,
                             ranges = IRanges(start = as.integer(CENgapAllAthila$gap_start),
                                              end = as.integer(CENgapAllAthila$gap_stop)),
                             strand = CENgapAllAthila$direction,
                             name = CENgapAllAthila$gap_name,
                             phylo = CENgapAllAthila$phylo)
CENgapAllAthilaGR <- unique(CENgapAllAthilaGR)
CENgapAllAthilaGR <- CENgapAllAthilaGR[seqnames(CENgapAllAthilaGR) %in% chrName]
CENgapAllAthilaGR <- sortSeqlevels(CENgapAllAthilaGR)
CENgapAllAthilaGR <- sort(CENgapAllAthilaGR, ignore.strand = TRUE)
CENgapAllAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllAthilaGR)),
                                  start = as.integer(start(CENgapAllAthilaGR)-1),
                                  end = as.integer(end(CENgapAllAthilaGR)),
                                  name = as.character(CENgapAllAthilaGR$name),
                                  score = as.character(CENgapAllAthilaGR$phylo),
                                  strand = as.character(strand(CENgapAllAthilaGR)))
write.table(CENgapAllAthila_bed,
            file = paste0(CENgapAllAthilaDir, "CENgapAllAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

## Convert CENgapAllNotAthila into GRanges and then BED format
#CENgapAllNotAthilaGR <- GRanges(seqnames = CENgapAllNotAthila$chr,
#                             ranges = IRanges(start = as.integer(CENgapAllNotAthila$gap_start),
#                                              end = as.integer(CENgapAllNotAthila$gap_stop)),
#                             strand = CENgapAllNotAthila$direction,
#                             name = CENgapAllNotAthila$gap_name,
#                             phylo = CENgapAllNotAthila$phylo)
#CENgapAllNotAthilaGR <- unique(CENgapAllNotAthilaGR)
#CENgapAllNotAthilaGR <- CENgapAllNotAthilaGR[seqnames(CENgapAllNotAthilaGR) %in% chrName]
#CENgapAllNotAthilaGR <- sortSeqlevels(CENgapAllNotAthilaGR)
#CENgapAllNotAthilaGR <- sort(CENgapAllNotAthilaGR, ignore.strand = TRUE)
#CENgapAllNotAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllNotAthilaGR)),
#                                  start = as.integer(start(CENgapAllNotAthilaGR)-1),
#                                  end = as.integer(end(CENgapAllNotAthilaGR)),
#                                  name = as.character(CENgapAllNotAthilaGR$name),
#                                  score = as.character(CENgapAllNotAthilaGR$phylo),
#                                  strand = as.character(strand(CENgapAllNotAthilaGR)))
#write.table(CENgapAllNotAthila_bed,
#            file = paste0(CENgapAllNotAthilaDir, "CENgapAllNotAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as CENgapAllAthila_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as CENgapAllAthilaGR
CENranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CENgapAllAthilaChrGR <- CENgapAllAthilaGR[seqnames(CENgapAllAthilaGR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  # Contract CENChrGR so that CENrandom loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(CENChrGR) <- end(CENChrGR)-max(width(CENgapAllAthilaChrGR))-2000
  start(CENChrGR) <- start(CENChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  CENranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {
                                                                start(CENChrGR[x]) : end(CENChrGR[x])
                                                              })),
                                         n = length(CENgapAllAthilaChrGR))
  CENranLocChrGR <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = CENranLocChrStart,
                                             width = width(CENgapAllAthilaChrGR)),
                            strand = strand(CENgapAllAthilaChrGR))
  CENranLocGR <- append(CENranLocGR, CENranLocChrGR)
}
stopifnot(identical(width(CENranLocGR), width(CENgapAllAthilaGR)))
stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CENgapAllAthilaGR))))
stopifnot(identical(strand(CENranLocGR), strand(CENgapAllAthilaGR)))
CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = start(CENranLocGR)-1,
                            end = end(CENranLocGR),
                            name = 1:length(CENranLocGR),
                            score = "NA",
                            strand = strand(CENranLocGR))
write.table(CENranLoc_bed,
            file = paste0(CENgapAllAthilaDir, "CENgapAllAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
