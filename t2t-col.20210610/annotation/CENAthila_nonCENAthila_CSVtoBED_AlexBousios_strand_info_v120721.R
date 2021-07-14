#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./CENAthila_nonCENAthila_CSVtoBED_AlexBousios_strand_info_v120721.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
CENgapAllDir <- "CENgapAll/"
CENgapAllAthilaDir <- "CENgapAllAthila/"
CENgapAllNotAthilaDir <- "CENgapAllNotAthila/"
system(paste0("[ -d ", CENAthilaDir, " ] || mkdir -p ", CENAthilaDir))
system(paste0("[ -d ", nonCENAthilaDir, " ] || mkdir -p ", nonCENAthilaDir))
system(paste0("[ -d ", CENgapDir, " ] || mkdir -p ", CENgapDir))
system(paste0("[ -d ", CENgapAllDir, " ] || mkdir -p ", CENgapAllDir))
system(paste0("[ -d ", CENgapAllAthilaDir, " ] || mkdir -p ", CENgapAllAthilaDir))
system(paste0("[ -d ", CENgapAllNotAthilaDir, " ] || mkdir -p ", CENgapAllNotAthilaDir))

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

## Load table of centromeric gaps (Athila-containing discontinuities in CEN180 sequences)
# and Athila sequence coordinates
# "rows 117 to 121 are some fragments which can be ignored."
tab <- read.csv("t2t-col.2021061.ATHILA_v100721.csv",
                header = T)[1:110,]
tab$phylo <- toupper(tab$phylo)
tab$centr_pos <- paste0(tab$centr_pos, "_centr")

CENAthila <- tab[tab$centr_pos == "in_centr",]
nonCENAthila <- tab[tab$centr_pos == "out_centr",]
CENgap <- tab[tab$centr_pos == "in_centr",]

# Load table of centromeric gaps (both those that do and do not contain Athila)
CENgapAll <- read.csv("CENgapSNVboolean_in_t2t-col.20210610.csv",
                      header = T)
colnames(CENgapAll)[1] <- "gap_name"
colnames(CENgapAll)[2] <- "chr"
CENgapAll$chr <- gsub(pattern = "^", replacement = "Chr", x = CENgapAll$chr)
colnames(CENgapAll)[3] <- "gap_start"
colnames(CENgapAll)[4] <- "gap_stop"
colnames(CENgapAll)[5] <- "variant"
CENgapAll$variant <- gsub(pattern = "yes", replacement = as.logical(1), x = CENgapAll$variant)
CENgapAll$variant <- gsub(pattern = "no", replacement = as.logical(0), x = CENgapAll$variant)

CENgapAllAthila <- CENgapAll[which(CENgapAll$gap_start %in% as.integer(CENgap$gap_start) |
                                   CENgapAll$gap_stop %in% as.integer(CENgap$gap_stop)),]
CENgapAllNotAthila <- CENgapAll[which(!(CENgapAll$gap_start %in% as.integer(CENgap$gap_start)) &
                                      !(CENgapAll$gap_stop %in% as.integer(CENgap$gap_stop))),]

# Convert CENAthila into GRanges
# Use genome_left_coord_FL and genome_right_coord_FL as
# element boundaries so that start and end coordinates shown in meta-profiles
# correspond to 5' LTR start and 3' LTR end coordinates
CENAthilaGR <- GRanges(seqnames = CENAthila$chr,
                       ranges = IRanges(start = CENAthila$genome_left_coord_FL,
                                        end = CENAthila$genome_right_coord_FL),
                       strand = CENAthila$direction,
                       name = CENAthila$TE_ID,
                       phylo = CENAthila$phylo)
CENAthilaGR <- unique(CENAthilaGR)
CENAthilaGR <- CENAthilaGR[seqnames(CENAthilaGR) %in% chrName]
CENAthilaGR <- sortSeqlevels(CENAthilaGR)
CENAthilaGR <- sort(CENAthilaGR, ignore.strand = TRUE)
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
                                            thichEnd = as.integer(0),
                                            itemRgb = ".")
  nonCENAthila_colofamily_bed[nonCENAthila_colofamily_bed$score == "ATHILA0",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[10])), collapse = ",")
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

# Convert CENgapAll into GRanges and then BED format
CENgapAllGR <- GRanges(seqnames = CENgapAll$chr,
                       ranges = IRanges(start = as.integer(CENgapAll$gap_start),
                                        end = as.integer(CENgapAll$gap_stop)),
                       strand = "*",
                       name = CENgapAll$gap_name,
                       variant = CENgapAll$variant)
CENgapAllGR <- unique(CENgapAllGR)
CENgapAllGR <- CENgapAllGR[seqnames(CENgapAllGR) %in% chrName]
CENgapAllGR <- sortSeqlevels(CENgapAllGR)
CENgapAllGR <- sort(CENgapAllGR, ignore.strand = TRUE)
CENgapAll_bed <- data.frame(chr = as.character(seqnames(CENgapAllGR)),
                            start = as.integer(start(CENgapAllGR)-1),
                            end = as.integer(end(CENgapAllGR)),
                            name = as.character(CENgapAllGR$name),
                            score = as.character(CENgapAllGR$variant),
                            strand = as.character(strand(CENgapAllGR)))
write.table(CENgapAll_bed,
            file = paste0(CENgapAllDir, "CENgapAll_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENgapAllAthila into GRanges and then BED format
CENgapAllAthilaGR <- GRanges(seqnames = CENgapAllAthila$chr,
                             ranges = IRanges(start = as.integer(CENgapAllAthila$gap_start),
                                              end = as.integer(CENgapAllAthila$gap_stop)),
                             strand = "*",
                             name = CENgapAllAthila$gap_name,
                             variant = CENgapAllAthila$variant)
CENgapAllAthilaGR <- unique(CENgapAllAthilaGR)
CENgapAllAthilaGR <- CENgapAllAthilaGR[seqnames(CENgapAllAthilaGR) %in% chrName]
CENgapAllAthilaGR <- sortSeqlevels(CENgapAllAthilaGR)
CENgapAllAthilaGR <- sort(CENgapAllAthilaGR, ignore.strand = TRUE)
CENgapAllAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllAthilaGR)),
                                  start = as.integer(start(CENgapAllAthilaGR)-1),
                                  end = as.integer(end(CENgapAllAthilaGR)),
                                  name = as.character(CENgapAllAthilaGR$name),
                                  score = as.character(CENgapAllAthilaGR$variant),
                                  strand = as.character(strand(CENgapAllAthilaGR)))
write.table(CENgapAllAthila_bed,
            file = paste0(CENgapAllAthilaDir, "CENgapAllAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert CENgapAllNotAthila into GRanges and then BED format
CENgapAllNotAthilaGR <- GRanges(seqnames = CENgapAllNotAthila$chr,
                                ranges = IRanges(start = as.integer(CENgapAllNotAthila$gap_start),
                                                 end = as.integer(CENgapAllNotAthila$gap_stop)),
                                strand = "*",
                                name = CENgapAllNotAthila$gap_name,
                                variant = CENgapAllNotAthila$variant)
CENgapAllNotAthilaGR <- unique(CENgapAllNotAthilaGR)
CENgapAllNotAthilaGR <- CENgapAllNotAthilaGR[seqnames(CENgapAllNotAthilaGR) %in% chrName]
CENgapAllNotAthilaGR <- sortSeqlevels(CENgapAllNotAthilaGR)
CENgapAllNotAthilaGR <- sort(CENgapAllNotAthilaGR, ignore.strand = TRUE)
CENgapAllNotAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllNotAthilaGR)),
                                     start = as.integer(start(CENgapAllNotAthilaGR)-1),
                                     end = as.integer(end(CENgapAllNotAthilaGR)),
                                     name = as.character(CENgapAllNotAthilaGR$name),
                                     score = as.character(CENgapAllNotAthilaGR$variant),
                                     strand = as.character(strand(CENgapAllNotAthilaGR)))
write.table(CENgapAllNotAthila_bed,
            file = paste0(CENgapAllNotAthilaDir, "CENgapAllNotAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Define function to select randomly positioned loci of the same
# width distribution as CENgapAll_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as CENgapAllGR
CENranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CENgapAllChrGR <- CENgapAllGR[seqnames(CENgapAllGR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  # Contract CENChrGR so that CENrandom loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(CENChrGR) <- end(CENChrGR)-max(width(CENgapAllChrGR))-2000
  start(CENChrGR) <- start(CENChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  CENranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {
                                                                start(CENChrGR[x]) : end(CENChrGR[x])
                                                              })),
                                         n = length(CENgapAllChrGR))
  CENranLocChrGR <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = CENranLocChrStart,
                                             width = width(CENgapAllChrGR)),
                            strand = strand(CENgapAllChrGR))
  CENranLocGR <- append(CENranLocGR, CENranLocChrGR)
}
stopifnot(identical(width(CENranLocGR), width(CENgapAllGR)))
stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CENgapAllGR))))
stopifnot(identical(strand(CENranLocGR), strand(CENgapAllGR)))
CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = start(CENranLocGR)-1,
                            end = end(CENranLocGR),
                            name = 1:length(CENranLocGR),
                            score = "NA",
                            strand = strand(CENranLocGR))
write.table(CENranLoc_bed,
            file = paste0(CENgapAllDir, "CENgapAll_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
