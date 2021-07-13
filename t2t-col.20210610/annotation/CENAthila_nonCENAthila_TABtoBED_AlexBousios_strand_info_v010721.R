#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./CENAthila_nonCENAthila_TABtoBED_AlexBousios_strand_info_v010721.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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

CENAthilaDir <- "CENAthila/"
nonCENAthilaDir <- "nonCENAthila/"
CENgapDir <- "CENgap/"
system(paste0("[ -d ", CENAthilaDir, " ] || mkdir -p ", CENAthilaDir))
system(paste0("[ -d ", nonCENAthilaDir, " ] || mkdir -p ", nonCENAthilaDir))
system(paste0("[ -d ", CENgapDir, " ] || mkdir -p ", CENgapDir))

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
tab <- read.table("t2t-col.20210610.fasta.blast.e180--91AthilaALL+info_0.75coverage_115Athila.bed",
                  header = F)
colnames(tab) <- c("chr", "start", "end", "catinfo", "score", "strand")

# Split "catinfo" column into multiple columns
tab$name <- sub(pattern = "_.+", replacement = "", x = tab$catinfo)
tab$position <- sub(pattern = "Chr\\d\\.\\d+-\\d+_", replacement = "", x = tab$catinfo)
tab$position <- sub(pattern = "_.+", replacement = "", x = tab$position)
tab$position <- paste0(tab$position, "_centr")
tab$phylo <- sub(pattern = "Chr\\d\\.\\d+-\\d+_[A-z]+_", replacement = "", x = tab$catinfo)
tab$phylo <- sub(pattern = "_.+", replacement = "", x = tab$phylo)
tab$phylo <- toupper(tab$phylo)
tab$compare <- sub(pattern = "Chr\\d\\.\\d+-\\d+_[A-z]+_Athila[0-9]", replacement = "", x = tab$catinfo)
tab$compare <- sub(pattern = "^[A-z]_", replacement = "", x = tab$compare)
tab$compare <- sub(pattern = "^_", replacement = "", x = tab$compare)
tab$compare <- sub(pattern = "_.+", replacement = "", x = tab$compare)
tab$share <- sub(pattern = ".+_", replacement = "", x = tab$catinfo)

CENAthila <- tab[tab$position == "in_centr",]
nonCENAthila <- tab[tab$position == "out_centr",]
# Placeholder CEN180 "gap_start" and "gap_end" coordinates (using ATHILA coordinates),
# while waiting for actual CEN180 gap coordinates from Ian and Alex
CENgap <- tab[tab$position == "in_centr",]
CENgap$gap_start <- CENgap$start
CENgap$gap_end <- CENgap$end
CENgap$gap_name <- CENgap$name

# Convert CENAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
CENAthilaGR <- GRanges(seqnames = CENAthila$chr,
                       ranges = IRanges(start = CENAthila$start,
                                        end = CENAthila$end),
                       strand = CENAthila$strand,
                       name = CENAthila$name,
                       phylo = CENAthila$phylo)
CENAthilaGR <- CENAthilaGR[seqnames(CENAthilaGR) %in% chrName]
CENAthila_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                         start = as.integer(start(CENAthilaGR)-1),
                         end = as.integer(end(CENAthilaGR)),
                         name = as.character(CENAthilaGR$name),
                         score = as.character(CENAthilaGR$phylo),
                         strand = as.character(strand(CENAthilaGR)))
write.table(CENAthila_bed,
            file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

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
                                         thichEnd = as.integer(0),
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
                          ranges = IRanges(start = nonCENAthila$start,
                                           end = nonCENAthila$end),
                          strand = nonCENAthila$strand,
                          name = nonCENAthila$name,
                          phylo = nonCENAthila$phylo)
nonCENAthilaGR <- nonCENAthilaGR[seqnames(nonCENAthilaGR) %in% chrName]
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
                                            thichEnd = as.integer(0),
                                            itemRgb = ".")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA0I",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[10])), collapse = ",")
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
              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Convert CENgap into GRanges and then BED format
CENgapGR <- GRanges(seqnames = CENgap$chr,
                    ranges = IRanges(start = CENgap$gap_start,
                                     end = CENgap$gap_end),
                    strand = CENgap$strand,
                    name = CENgap$gap_name,
                    phylo = CENgap$phylo)
CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
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
