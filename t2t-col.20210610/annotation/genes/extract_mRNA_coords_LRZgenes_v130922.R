#!/applications/R/R-4.0.0/bin/Rscript

# Create tables of LRZ genes and nonLRZ genes
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_mRNA_coords_LRZgenes_v130922.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(data.table)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CEN <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.centromeres", header = T)
CENstart <- CEN$start[which(fai$V1 %in% chrName)]
CENend <- CEN$end[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")
nonCENGR <- GRanges(seqnames = rep(seqnames(CENGR), 2),
                    ranges = IRanges(start = c(rep(1, length(start(CENGR))),
                                               end(CENGR)+1),
                                     end = c(start(CENGR)-1,
                                             chrLens)),
                    strand = "*")

# Get genes in low-recombining zones (LRZgenes)
mRNA_rep_LRZ <- read.csv("cold.genes.for.andy.csv", header = T)[,c(2:9, 11)]
colnames(mRNA_rep_LRZ) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "Parent")
mRNA_rep_LRZ <- mRNA_rep_LRZ[mRNA_rep_LRZ$seqid %in% chrName,]
#mRNA_rep_LRZ[mRNA_rep_LRZ == "."] <- "NA"
nrow(mRNA_rep_LRZ)
#[1] 712
head(mRNA_rep_LRZ)

# Get genes not in low-recombining zones (nonLRZgenes)
mRNA_rep_nonLRZ <- read.csv("hot.genes.for.andy.csv", header = T)[,c(2:10)]
colnames(mRNA_rep_nonLRZ) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "Parent")
mRNA_rep_nonLRZ$Parent <- gsub("ID=", "", mRNA_rep_nonLRZ$Parent)
mRNA_rep_nonLRZ$Parent <- gsub(";.+", "", mRNA_rep_nonLRZ$Parent)
mRNA_rep_nonLRZ <- mRNA_rep_nonLRZ[mRNA_rep_nonLRZ$seqid %in% chrName,]
#mRNA_rep_nonLRZ[mRNA_rep_nonLRZ == "."] <- "NA"
nrow(mRNA_rep_nonLRZ)
#[1] 27686
head(mRNA_rep_nonLRZ)


write.table(mRNA_rep_LRZ,
            file = paste0("t2t-col.20210610_representative_mRNA_LRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_LRZ_bed <- data.frame(chr = as.character(mRNA_rep_LRZ[,1]),
                               start = as.integer(mRNA_rep_LRZ[,4]-1),
                               end = as.integer(mRNA_rep_LRZ[,5]),
                               name = as.character(mRNA_rep_LRZ[,9]),
                               score = rep("NA", nrow(mRNA_rep_LRZ)),
                               strand = as.character(mRNA_rep_LRZ[,7]))
write.table(mRNA_rep_LRZ_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_LRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_LRZGR <- GRanges(seqnames = mRNA_rep_LRZ$seqid,
                          ranges = IRanges(start = mRNA_rep_LRZ$start,
                                           end = mRNA_rep_LRZ$end),
                          strand = mRNA_rep_LRZ$strand)

write.table(mRNA_rep_nonLRZ,
            file = paste0("t2t-col.20210610_representative_mRNA_nonLRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_nonLRZ_bed <- data.frame(chr = as.character(mRNA_rep_nonLRZ[,1]),
                               start = as.integer(mRNA_rep_nonLRZ[,4]-1),
                               end = as.integer(mRNA_rep_nonLRZ[,5]),
                               name = as.character(mRNA_rep_nonLRZ[,9]),
                               score = rep("NA", nrow(mRNA_rep_nonLRZ)),
                               strand = as.character(mRNA_rep_nonLRZ[,7]))
write.table(mRNA_rep_nonLRZ_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_nonLRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_nonLRZGR <- GRanges(seqnames = mRNA_rep_nonLRZ$seqid,
                          ranges = IRanges(start = mRNA_rep_nonLRZ$start,
                                           end = mRNA_rep_nonLRZ$end),
                          strand = mRNA_rep_nonLRZ$strand)



# Define function to select randomly positioned loci of the same
# width distribution as mRNA_rep_LRZGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as mRNA_rep_LRZGR
chrs <- chrs[chrs %in% chrName]
CENranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_CENGR <- CENGR[seqnames(CENGR) == chrs[j]]

  chr_mRNA_rep_LRZGR <- mRNA_rep_LRZGR[seqnames(mRNA_rep_LRZGR) == chrs[j]]

  if(length(chr_mRNA_rep_LRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_CENGR), function(x) {
                                                                    ( start(chr_CENGR[x]) + max(width(chr_mRNA_rep_LRZGR)) + 2000 ) :
                                                                    ( end(chr_CENGR[x]) - max(width(chr_mRNA_rep_LRZGR)) - 2000 )
                                                                  })),
                                             n = length(chr_mRNA_rep_LRZGR))
    chr_CENranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_CENranLoc_Start,
                                                width = width(chr_mRNA_rep_LRZGR)),
                               strand = strand(chr_mRNA_rep_LRZGR))
    CENranLocGR <- append(CENranLocGR, chr_CENranLocGR)
  }
}
stopifnot(identical(width(CENranLocGR), width(mRNA_rep_LRZGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# nonCENranLocGR contains the same number of loci per chromosome as mRNA_rep_nonLRZGR
chrs <- chrs[chrs %in% chrName]
nonCENranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_nonCENGR <- nonCENGR[seqnames(nonCENGR) == chrs[j]]

  chr_mRNA_rep_nonLRZGR <- mRNA_rep_nonLRZGR[seqnames(mRNA_rep_nonLRZGR) == chrs[j]]

  if(length(chr_mRNA_rep_nonLRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_nonCENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_nonCENGR), function(x) {
                                                                       ( start(chr_nonCENGR[x]) + max(width(chr_mRNA_rep_nonLRZGR)) + 2000 ) :
                                                                       ( end(chr_nonCENGR[x]) - max(width(chr_mRNA_rep_nonLRZGR)) - 2000 )
                                                                     })),
                                                n = length(chr_mRNA_rep_nonLRZGR))
    chr_nonCENranLocGR <- GRanges(seqnames = chrs[j],
                                  ranges = IRanges(start = chr_nonCENranLoc_Start,
                                                   width = width(chr_mRNA_rep_nonLRZGR)),
                                  strand = strand(chr_mRNA_rep_nonLRZGR))
    nonCENranLocGR <- append(nonCENranLocGR, chr_nonCENranLocGR)
  }
}
stopifnot(identical(width(nonCENranLocGR), width(mRNA_rep_nonLRZGR)))


CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = as.integer(start(CENranLocGR)-1),
                            end = as.integer(end(CENranLocGR)),
                            name = as.integer(1:length(CENranLocGR)),
                            score = rep("NA", length(CENranLocGR)),
                            strand = as.character(strand(CENranLocGR)))
write.table(CENranLoc_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_LRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

nonCENranLoc_bed <- data.frame(chr = as.character(seqnames(nonCENranLocGR)),
                               start = as.integer(start(nonCENranLocGR)-1),
                               end = as.integer(end(nonCENranLocGR)),
                               name = as.integer(1:length(nonCENranLocGR)),
                               score = rep("NA", length(nonCENranLocGR)),
                               strand = as.character(strand(nonCENranLocGR)))
write.table(nonCENranLoc_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_nonLRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
