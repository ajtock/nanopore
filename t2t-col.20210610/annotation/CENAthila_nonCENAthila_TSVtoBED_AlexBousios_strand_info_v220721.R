#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./CENAthila_nonCENAthila_TSVtoBED_AlexBousios_strand_info_v220721.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
CENsoloLTRDir <- "CENsoloLTR/"
CENfragmentAthilaDir <- "CENfragmentAthila/"
CENgapDir <- "CENgap/"
CENgapAllDir <- "CENgapAll/"
CENgapAllAthilaDir <- "CENgapAllAthila/"
CENgapAllNotAthilaDir <- "CENgapAllNotAthila/"
system(paste0("[ -d ", CENAthilaDir, " ] || mkdir -p ", CENAthilaDir))
system(paste0("[ -d ", nonCENAthilaDir, " ] || mkdir -p ", nonCENAthilaDir))
system(paste0("[ -d ", CENsoloLTRDir, " ] || mkdir -p ", CENsoloLTRDir))
system(paste0("[ -d ", CENfragmentAthilaDir, " ] || mkdir -p ", CENfragmentAthilaDir))
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
tab <- read.table("t2t-col.20210610_ATHILA.inout.FL.soloLTR.fragments.tsv",
                  header = T)
tab$phylo <- toupper(tab$phylo)
tab$centr_pos <- paste0(tab$centr_pos, "_centr")

CENAthila <- tab[tab$centr_pos == "in_centr" & tab$quality == "intact",]
nonCENAthila <- tab[tab$centr_pos == "out_centr" & tab$quality == "intact",]
CENsoloLTR <- tab[tab$centr_pos == "in_centr" & tab$quality == "solo",]
CENfragmentAthila <- tab[tab$centr_pos == "in_centr" & tab$quality == "fragment",]
CENgap <- tab[tab$centr_pos == "in_centr" & tab$quality == "intact",]

# Load table of centromeric gaps (both those that do contain Athila)
CENgapAllAthila <- read.csv("CEN180_adjusted_gaps.csv",
                      header = T)
colnames(CENgapAllAthila)[1] <- "chr"
colnames(CENgapAllAthila)[2] <- "gap_name"
colnames(CENgapAllAthila)[5] <- "old_gap_start"
colnames(CENgapAllAthila)[6] <- "old_gap_stop"
colnames(CENgapAllAthila)[3] <- "gap_start"
colnames(CENgapAllAthila)[4] <- "gap_stop"
colnames(CENgapAllAthila)[9] <- "old_gap_width"
colnames(CENgapAllAthila)[10] <- "gap_width"
colnames(CENgapAllAthila)[11] <- "direction"
colnames(CENgapAllAthila)[12] <- "phylo"
colnames(CENgapAllAthila)[13] <- "quality"
CENgapAllAthila$chr <- gsub(pattern = "^", replacement = "Chr", x = CENgapAllAthila$chr)
CENgapAllAthila$phylo <- toupper(CENgapAllAthila$phylo)

#CENgapAllAthila <- CENgapAll[which(CENgapAll$gap_start %in% as.integer(CENgap$gap_start) |
#                                   CENgapAll$gap_stop %in% as.integer(CENgap$gap_stop)),]
#CENgapAllNotAthila <- CENgapAll[which(!(CENgapAll$gap_start %in% as.integer(CENgap$gap_start)) &
#                                      !(CENgapAll$gap_stop %in% as.integer(CENgap$gap_stop))),]

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
#write.table(CENAthila_bed,
#            file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)

if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  CENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(CENAthilaGR)),
                                       start = as.integer(start(CENAthilaGR)-1),
                                       end = as.integer(end(CENAthilaGR)),
                                       name = as.character(CENAthilaGR$name),
                                       score = as.integer(0),
                                       strand = as.character(strand(CENAthilaGR)))
#  write.table(CENAthila_nofamily_bed,
#              file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
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
#  write.table(CENAthila_colofamily_bed,
#              file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
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
#write.table(nonCENAthila_bed,
#            file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)

if(length(chrName) > 1) {
  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
  nonCENAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(nonCENAthilaGR)),
                                          start = as.integer(start(nonCENAthilaGR)-1),
                                          end = as.integer(end(nonCENAthilaGR)),
                                          name = as.character(nonCENAthilaGR$name),
                                          score = as.integer(0),
                                          strand = as.character(strand(nonCENAthilaGR)))
#  write.table(nonCENAthila_nofamily_bed,
#              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
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
#  write.table(nonCENAthila_colofamily_bed,
#              file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
}


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
# CENranLocGR contains the same number of loci per chromosome as CENAthilaGR
CENranLocGR <- GRanges()
nonCENranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_CENGR <- CENGR[seqnames(CENGR) == chrs[j]]

  chr_CENAthilaGR <- CENAthilaGR[seqnames(CENAthilaGR) == chrs[j]]
  if(length(chr_CENAthilaGR) > 0) {
    ## Contract chr_CENGR so that CENranLocGR and 2-kb flanking regions
    ## do not extend beyond centromeric coordinates
    #end(chr_CENGR) <- end(chr_CENGR)-max(width(chr_CENAthilaGR))-2000
    #start(chr_CENGR) <- start(chr_CENGR)+2000

    # Define seed so that random selections are reproducible
    set.seed(76492749)
    chr_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_CENGR), function(x) {
                                                                    start(chr_CENGR[x]) : end(chr_CENGR[x])
                                                                  })),
                                             n = length(chr_CENAthilaGR))
    chr_CENranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_CENranLoc_Start,
                                                width = width(chr_CENAthilaGR)),
                               strand = strand(chr_CENAthilaGR))
    CENranLocGR <- append(CENranLocGR, chr_CENranLocGR)
  }

  chr_nonCENAthilaGR <- nonCENAthilaGR[seqnames(nonCENAthilaGR) == chrs[j]]
  if(length(chr_nonCENAthilaGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(76492749)
    chr_nonCENranLoc_Start <- ranLocStartSelect(coordinates = unique(unlist(lapply(seq_along(chr_CENGR), function(x) {
                                                                       c(1:chrLens[j])[
                                                                         which(
                                                                           !( 1:chrLens[j] %in%
                                                                              start(chr_CENGR[x]) : end(chr_CENGR[x])
                                                                            )
                                                                         )
                                                                       ]
                                                                     }))),
                                                n = length(chr_nonCENAthilaGR))
    chr_nonCENranLocGR <- GRanges(seqnames = chrs[j],
                                  ranges = IRanges(start = chr_nonCENranLoc_Start,
                                                   width = width(chr_nonCENAthilaGR)),
                                  strand = strand(chr_nonCENAthilaGR))
    nonCENranLocGR <- append(nonCENranLocGR, chr_nonCENranLocGR)
  }
}
stopifnot(identical(width(CENranLocGR), width(CENAthilaGR)))
stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CENAthilaGR))))
stopifnot(identical(strand(CENranLocGR), strand(CENAthilaGR)))
CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = start(CENranLocGR)-1,
                            end = end(CENranLocGR),
                            name = 1:length(CENranLocGR),
                            score = CENAthilaGR$phylo,
                            strand = strand(CENranLocGR))
write.table(CENranLoc_bed,
            file = paste0(CENAthilaDir, "CENAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

stopifnot(identical(width(nonCENranLocGR), width(nonCENAthilaGR)))
stopifnot(identical(as.character(seqnames(nonCENranLocGR)), as.character(seqnames(nonCENAthilaGR))))
stopifnot(identical(strand(nonCENranLocGR), strand(nonCENAthilaGR)))
nonCENranLoc_bed <- data.frame(chr = as.character(seqnames(nonCENranLocGR)),
                               start = start(nonCENranLocGR)-1,
                               end = end(nonCENranLocGR),
                               name = 1:length(nonCENranLocGR),
                               score = nonCENAthilaGR$phylo,
                               strand = strand(nonCENranLocGR))
write.table(nonCENranLoc_bed,
            file = paste0(nonCENAthilaDir, "nonCENAthila_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_nonCENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)



## Convert CENsoloLTR into GRanges
## Use genome_left_coord_FL and genome_right_coord_FL as
## element boundaries
#CENsoloLTRGR <- GRanges(seqnames = CENsoloLTR$chr,
#                        ranges = IRanges(start = CENsoloLTR$genome_left_coord_FL,
#                                         end = CENsoloLTR$genome_right_coord_FL),
#                        strand = CENsoloLTR$direction,
#                        name = CENsoloLTR$TE_ID,
#                        phylo = CENsoloLTR$phylo)
#CENsoloLTRGR <- unique(CENsoloLTRGR)
#CENsoloLTRGR <- CENsoloLTRGR[seqnames(CENsoloLTRGR) %in% chrName]
#CENsoloLTRGR <- sortSeqlevels(CENsoloLTRGR)
#CENsoloLTRGR <- sort(CENsoloLTRGR, ignore.strand = TRUE)
#CENsoloLTR_bed <- data.frame(chr = as.character(seqnames(CENsoloLTRGR)),
#                             start = as.integer(start(CENsoloLTRGR)-1),
#                             end = as.integer(end(CENsoloLTRGR)),
#                             name = as.character(CENsoloLTRGR$name),
#                             score = as.character(CENsoloLTRGR$phylo),
#                             strand = as.character(strand(CENsoloLTRGR)))
#write.table(CENsoloLTR_bed,
#            file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
#
#if(length(chrName) > 1) {
#  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
#  CENsoloLTR_nofamily_bed <- data.frame(chr = as.character(seqnames(CENsoloLTRGR)),
#                                        start = as.integer(start(CENsoloLTRGR)-1),
#                                        end = as.integer(end(CENsoloLTRGR)),
#                                        name = as.character(CENsoloLTRGR$name),
#                                        score = as.integer(0),
#                                        strand = as.character(strand(CENsoloLTRGR)))
#  write.table(CENsoloLTR_nofamily_bed,
#              file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
#  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
#  CENfams <- sort(unique(CENsoloLTR_bed$score))
#  print(CENfams)
#  CENsoloLTR_colofamily_bed <- data.frame(CENsoloLTR_bed,
#                                          thickStart = as.integer(0),
#                                          thickEnd = as.integer(0),
#                                          itemRgb = ".")
#  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
#  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
#  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA4-4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
#  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[4])), collapse = ",")
#  CENsoloLTR_colofamily_bed[CENsoloLTR_colofamily_bed$score == "ATHILA6",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
#  
#  CENsoloLTR_colofamily_bed$score <- as.integer(0)
#  write.table(CENsoloLTR_colofamily_bed,
#              file = paste0(CENsoloLTRDir, "CENsoloLTR_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
#}
#
## Convert CENfragmentAthila into GRanges
## Use genome_left_coord_FL and genome_right_coord_FL as
## element boundaries
#CENfragmentAthilaGR <- GRanges(seqnames = CENfragmentAthila$chr,
#                               ranges = IRanges(start = CENfragmentAthila$genome_left_coord_FL,
#                                                end = CENfragmentAthila$genome_right_coord_FL),
#                               strand = CENfragmentAthila$direction,
#                               name = CENfragmentAthila$TE_ID,
#                               phylo = CENfragmentAthila$phylo)
#CENfragmentAthilaGR <- unique(CENfragmentAthilaGR)
#CENfragmentAthilaGR <- CENfragmentAthilaGR[seqnames(CENfragmentAthilaGR) %in% chrName]
#CENfragmentAthilaGR <- sortSeqlevels(CENfragmentAthilaGR)
#CENfragmentAthilaGR <- sort(CENfragmentAthilaGR, ignore.strand = TRUE)
#CENfragmentAthila_bed <- data.frame(chr = as.character(seqnames(CENfragmentAthilaGR)),
#                                    start = as.integer(start(CENfragmentAthilaGR)-1),
#                                    end = as.integer(end(CENfragmentAthilaGR)),
#                                    name = as.character(CENfragmentAthilaGR$name),
#                                    score = as.character(CENfragmentAthilaGR$phylo),
#                                    strand = as.character(strand(CENfragmentAthilaGR)))
#write.table(CENfragmentAthila_bed,
#            file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
#
#if(length(chrName) > 1) {
#  # Write BED without family information, for use with pyGenomeTracks (BED specified in *.ini file)
#  CENfragmentAthila_nofamily_bed <- data.frame(chr = as.character(seqnames(CENfragmentAthilaGR)),
#                                               start = as.integer(start(CENfragmentAthilaGR)-1),
#                                               end = as.integer(end(CENfragmentAthilaGR)),
#                                               name = as.character(CENfragmentAthilaGR$name),
#                                               score = as.integer(0),
#                                               strand = as.character(strand(CENfragmentAthilaGR)))
#  write.table(CENfragmentAthila_nofamily_bed,
#              file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_nofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
#  # Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
#  CENfams <- sort(unique(CENfragmentAthila_bed$score))
#  print(CENfams)
#  CENfragmentAthila_colofamily_bed <- data.frame(CENfragmentAthila_bed,
#                                                 thickStart = as.integer(0),
#                                                 thickEnd = as.integer(0),
#                                                 itemRgb = ".")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[9])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[8])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[7])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA4-4C",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[6])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[3])), collapse = ",")
#  CENfragmentAthila_colofamily_bed[CENfragmentAthila_colofamily_bed$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich10()(10)[2])), collapse = ",")
#  
#  CENfragmentAthila_colofamily_bed$score <- as.integer(0)
#  write.table(CENfragmentAthila_colofamily_bed,
#              file = paste0(CENfragmentAthilaDir, "CENfragmentAthila_in_t2t-col.20210610_",
#                            paste0(chrName, collapse = "_"), "_colofamily.bed"),
#              quote = F, sep = "\t", row.names = F, col.names = F)
#}
#
## Convert CENgap into GRanges and then BED format
#CENgapGR <- GRanges(seqnames = CENgap$chr,
#                    ranges = IRanges(start = as.integer(CENgap$gap_start),
#                                     end = as.integer(CENgap$gap_stop)),
#                    strand = CENgap$direction,
#                    name = CENgap$gap_name,
#                    phylo = CENgap$phylo)
#CENgapGR <- unique(CENgapGR)
#CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
#CENgapGR <- sortSeqlevels(CENgapGR)
#CENgapGR <- sort(CENgapGR, ignore.strand = TRUE)
#CENgap_bed <- data.frame(chr = as.character(seqnames(CENgapGR)),
#                         start = as.integer(start(CENgapGR)-1),
#                         end = as.integer(end(CENgapGR)),
#                         name = as.character(CENgapGR$name),
#                         score = as.character(CENgapGR$phylo),
#                         strand = as.character(strand(CENgapGR)))
#write.table(CENgap_bed,
#            file = paste0(CENgapDir, "CENgap_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
#
### Convert CENgapAll into GRanges and then BED format
##CENgapAllGR <- GRanges(seqnames = CENgapAll$chr,
##                       ranges = IRanges(start = as.integer(CENgapAll$gap_start),
##                                        end = as.integer(CENgapAll$gap_stop)),
##                       strand = CENgapAll$direction,
##                       name = CENgapAll$gap_name
##                       phylo = CENgapAll$phylo)
##CENgapAllGR <- unique(CENgapAllGR)
##CENgapAllGR <- CENgapAllGR[seqnames(CENgapAllGR) %in% chrName]
##CENgapAllGR <- sortSeqlevels(CENgapAllGR)
##CENgapAllGR <- sort(CENgapAllGR, ignore.strand = TRUE)
##CENgapAll_bed <- data.frame(chr = as.character(seqnames(CENgapAllGR)),
##                            start = as.integer(start(CENgapAllGR)-1),
##                            end = as.integer(end(CENgapAllGR)),
##                            name = as.character(CENgapAllGR$name),
##                            score = as.character(CENgapAllGR$phylo),
##                            strand = as.character(strand(CENgapAllGR)))
##write.table(CENgapAll_bed,
##            file = paste0(CENgapAllDir, "CENgapAll_in_t2t-col.20210610_",
##                          paste0(chrName, collapse = "_"), ".bed"),
##            quote = F, sep = "\t", row.names = F, col.names = F)
#
## Convert CENgapAllAthila into GRanges and then BED format
#CENgapAllAthilaGR <- GRanges(seqnames = CENgapAllAthila$chr,
#                             ranges = IRanges(start = as.integer(CENgapAllAthila$gap_start),
#                                              end = as.integer(CENgapAllAthila$gap_stop)),
#                             strand = CENgapAllAthila$direction,
#                             name = CENgapAllAthila$gap_name,
#                             phylo = CENgapAllAthila$phylo)
#CENgapAllAthilaGR <- unique(CENgapAllAthilaGR)
#CENgapAllAthilaGR <- CENgapAllAthilaGR[seqnames(CENgapAllAthilaGR) %in% chrName]
#CENgapAllAthilaGR <- sortSeqlevels(CENgapAllAthilaGR)
#CENgapAllAthilaGR <- sort(CENgapAllAthilaGR, ignore.strand = TRUE)
#CENgapAllAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllAthilaGR)),
#                                  start = as.integer(start(CENgapAllAthilaGR)-1),
#                                  end = as.integer(end(CENgapAllAthilaGR)),
#                                  name = as.character(CENgapAllAthilaGR$name),
#                                  score = as.character(CENgapAllAthilaGR$phylo),
#                                  strand = as.character(strand(CENgapAllAthilaGR)))
#write.table(CENgapAllAthila_bed,
#            file = paste0(CENgapAllAthilaDir, "CENgapAllAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), ".bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
#
### Convert CENgapAllNotAthila into GRanges and then BED format
##CENgapAllNotAthilaGR <- GRanges(seqnames = CENgapAllNotAthila$chr,
##                             ranges = IRanges(start = as.integer(CENgapAllNotAthila$gap_start),
##                                              end = as.integer(CENgapAllNotAthila$gap_stop)),
##                             strand = CENgapAllNotAthila$direction,
##                             name = CENgapAllNotAthila$gap_name,
##                             phylo = CENgapAllNotAthila$phylo)
##CENgapAllNotAthilaGR <- unique(CENgapAllNotAthilaGR)
##CENgapAllNotAthilaGR <- CENgapAllNotAthilaGR[seqnames(CENgapAllNotAthilaGR) %in% chrName]
##CENgapAllNotAthilaGR <- sortSeqlevels(CENgapAllNotAthilaGR)
##CENgapAllNotAthilaGR <- sort(CENgapAllNotAthilaGR, ignore.strand = TRUE)
##CENgapAllNotAthila_bed <- data.frame(chr = as.character(seqnames(CENgapAllNotAthilaGR)),
##                                  start = as.integer(start(CENgapAllNotAthilaGR)-1),
##                                  end = as.integer(end(CENgapAllNotAthilaGR)),
##                                  name = as.character(CENgapAllNotAthilaGR$name),
##                                  score = as.character(CENgapAllNotAthilaGR$phylo),
##                                  strand = as.character(strand(CENgapAllNotAthilaGR)))
##write.table(CENgapAllNotAthila_bed,
##            file = paste0(CENgapAllNotAthilaDir, "CENgapAllNotAthila_in_t2t-col.20210610_",
##                          paste0(chrName, collapse = "_"), ".bed"),
##            quote = F, sep = "\t", row.names = F, col.names = F)
#
## Define function to select randomly positioned loci of the same
## width distribution as CENgapAllAthila_bed
#ranLocStartSelect <- function(coordinates, n) {
#  sample(x = coordinates,
#         size = n,
#         replace = FALSE)
#}
#
## Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
#options(scipen = 100)
#
## Apply ranLocStartSelect() on a per-chromosome basis so that
## CENranLocGR contains the same number of loci per chromosome as CENgapAllAthilaGR
#CENranLocGR <- GRanges()
#for(i in 1:length(chrs)) {
#  CENgapAllAthilaChrGR <- CENgapAllAthilaGR[seqnames(CENgapAllAthilaGR) == chrs[i]]
#  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
#  # Contract CENChrGR so that CENrandom loci and 2-kb flanking regions
#  # do not extend beyond chromosome ends
#  end(CENChrGR) <- end(CENChrGR)-max(width(CENgapAllAthilaChrGR))-2000
#  start(CENChrGR) <- start(CENChrGR)+2000
#  # Define seed so that random selections are reproducible
#  set.seed(76492749)
#  CENranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {
#                                                                start(CENChrGR[x]) : end(CENChrGR[x])
#                                                              })),
#                                         n = length(CENgapAllAthilaChrGR))
#  CENranLocChrGR <- GRanges(seqnames = chrs[i],
#                            ranges = IRanges(start = CENranLocChrStart,
#                                             width = width(CENgapAllAthilaChrGR)),
#                            strand = strand(CENgapAllAthilaChrGR))
#  CENranLocGR <- append(CENranLocGR, CENranLocChrGR)
#}
#stopifnot(identical(width(CENranLocGR), width(CENgapAllAthilaGR)))
#stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CENgapAllAthilaGR))))
#stopifnot(identical(strand(CENranLocGR), strand(CENgapAllAthilaGR)))
#CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
#                            start = start(CENranLocGR)-1,
#                            end = end(CENranLocGR),
#                            name = 1:length(CENranLocGR),
#                            score = "NA",
#                            strand = strand(CENranLocGR))
#write.table(CENranLoc_bed,
#            file = paste0(CENgapAllAthilaDir, "CENgapAllAthila_in_t2t-col.20210610_",
#                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
