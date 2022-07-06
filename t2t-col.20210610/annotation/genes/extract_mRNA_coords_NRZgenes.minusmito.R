#!/applications/R/R-4.0.0/bin/Rscript

# Create table of representative genes
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_mRNA_coords_NRZgenes.minusmito.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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

# Load table of gene coordinates in t2t-col.20210610 (GFF3 derived from Liftoff tool)
genes <- readGFF("t2t-col.20210610.genes.gff3")
genes <- genes[genes$seqid %in% chrName,]
mRNA <- genes[genes$type == "mRNA",]
print(dim(mRNA))
#[1] 35234    20

# Obtain frequency of occurrence of each gene parent ID
n_occur <- data.frame(table(unlist(mRNA$Parent)))

# Obtain mRNA records for which the gene parent ID occurs only once
mRNA_unique <-  as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq == 1],
  ]
)
# Obtain mRNA records for which the gene parent ID occurs more than once
mRNA_multi <- as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq > 1],
  ]
)

# For each gene parent ID in mRNA_multi, obtain the mRNA record with the
# longest transcript
# If multiple mRNA records have the longest transcript,
# keep the first reported one only
mRNA_multi_list <- mclapply(seq_along(mRNA_multi[,1]), function(h) {
  mRNA_multi_ID_all <- mRNA_multi[ unlist(mRNA_multi$Parent)
                         == unlist(mRNA_multi[h,]$Parent),
                       ]
  mRNA_multi_ID_all[ mRNA_multi_ID_all$end-mRNA_multi_ID_all$start
    == max(mRNA_multi_ID_all$end-mRNA_multi_ID_all$start),
  ][1,]
}, mc.cores = detectCores())

# Collapse mRNA_multi_list into single data.frame and remove duplicates
mRNA_multi_dup <- rbindlist(mRNA_multi_list)
mRNA_multi_rep <- unique(as.data.frame(mRNA_multi_dup))

# Combine into one representative set of mRNA entries, order,
# and output in GFF3 and BED formats
mRNA_rep <- rbind(mRNA_unique, mRNA_multi_rep)
mRNA_rep <- mRNA_rep[ order(mRNA_rep$seqid,
                            mRNA_rep$start,
                            mRNA_rep$end), ]

# Get mRNA_rep in non-recombining zones (NRZgenes)
NRZgenes <- read.table("t2t-col.20210610.NRZgenes.minusmito.txt", header = T)[,1]
length(NRZgenes)
#[1] 579
mRNA_rep_NRZ <- mRNA_rep[ which(unlist(mRNA_rep$Parent) %in% NRZgenes), ]
length(unique(unlist(mRNA_rep_NRZ$Parent)))
nrow(mRNA_rep_NRZ)
#[1] 551

# Get mRNA_rep not in non-recombining zones (NRZgenes)
mRNA_rep_notNRZ <- mRNA_rep[ -which(unlist(mRNA_rep$Parent) %in% NRZgenes), ]
length(unique(unlist(mRNA_rep_notNRZ$Parent)))
nrow(mRNA_rep_notNRZ)
#[1] 26707 
stopifnot( nrow(mRNA_rep) == ( nrow(mRNA_rep_NRZ) + nrow(mRNA_rep_notNRZ) ) )

write.table(mRNA_rep_NRZ[,1:9],
            file = paste0("t2t-col.20210610_representative_mRNA_NRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_NRZ_bed <- data.frame(chr = as.character(mRNA_rep_NRZ[,1]),
                               start = as.integer(mRNA_rep_NRZ[,4]-1),
                               end = as.integer(mRNA_rep_NRZ[,5]),
                               name = as.character(mRNA_rep_NRZ[,9]),
                               score = as.numeric(mRNA_rep_NRZ[,6]),
                               strand = as.character(mRNA_rep_NRZ[,7]))
write.table(mRNA_rep_NRZ_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_NRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_NRZGR <- GRanges(seqnames = mRNA_rep_NRZ$seqid,
                          ranges = IRanges(start = mRNA_rep_NRZ$start,
                                           end = mRNA_rep_NRZ$end),
                          strand = mRNA_rep_NRZ$strand)

write.table(mRNA_rep_notNRZ[,1:9],
            file = paste0("t2t-col.20210610_representative_mRNA_notNRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_notNRZ_bed <- data.frame(chr = as.character(mRNA_rep_notNRZ[,1]),
                               start = as.integer(mRNA_rep_notNRZ[,4]-1),
                               end = as.integer(mRNA_rep_notNRZ[,5]),
                               name = as.character(mRNA_rep_notNRZ[,9]),
                               score = as.numeric(mRNA_rep_notNRZ[,6]),
                               strand = as.character(mRNA_rep_notNRZ[,7]))
write.table(mRNA_rep_notNRZ_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_notNRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_notNRZGR <- GRanges(seqnames = mRNA_rep_notNRZ$seqid,
                          ranges = IRanges(start = mRNA_rep_notNRZ$start,
                                           end = mRNA_rep_notNRZ$end),
                          strand = mRNA_rep_notNRZ$strand)


# Define function to select randomly positioned loci of the same
# width distribution as mRNA_rep_NRZGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as mRNA_rep_NRZGR
chrs <- chrs[chrs %in% chrName]
CENranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_CENGR <- CENGR[seqnames(CENGR) == chrs[j]]

  chr_mRNA_rep_NRZGR <- mRNA_rep_NRZGR[seqnames(mRNA_rep_NRZGR) == chrs[j]]

  if(length(chr_mRNA_rep_NRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_CENGR), function(x) {
                                                                    ( start(chr_CENGR[x]) + max(width(chr_mRNA_rep_NRZGR)) + 2000 ) :
                                                                    ( end(chr_CENGR[x]) - max(width(chr_mRNA_rep_NRZGR)) - 2000 )
                                                                  })),
                                             n = length(chr_mRNA_rep_NRZGR))
    chr_CENranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_CENranLoc_Start,
                                                width = width(chr_mRNA_rep_NRZGR)),
                               strand = strand(chr_mRNA_rep_NRZGR))
    CENranLocGR <- append(CENranLocGR, chr_CENranLocGR)
  }
}
stopifnot(identical(width(CENranLocGR), width(mRNA_rep_NRZGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# nonCENranLocGR contains the same number of loci per chromosome as mRNA_rep_notNRZGR
chrs <- chrs[chrs %in% chrName]
nonCENranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_nonCENGR <- nonCENGR[seqnames(nonCENGR) == chrs[j]]

  chr_mRNA_rep_notNRZGR <- mRNA_rep_notNRZGR[seqnames(mRNA_rep_notNRZGR) == chrs[j]]

  if(length(chr_mRNA_rep_notNRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_nonCENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_nonCENGR), function(x) {
                                                                       ( start(chr_nonCENGR[x]) + max(width(chr_mRNA_rep_notNRZGR)) + 2000 ) :
                                                                       ( end(chr_nonCENGR[x]) - max(width(chr_mRNA_rep_notNRZGR)) - 2000 )
                                                                     })),
                                                n = length(chr_mRNA_rep_notNRZGR))
    chr_nonCENranLocGR <- GRanges(seqnames = chrs[j],
                                  ranges = IRanges(start = chr_nonCENranLoc_Start,
                                                   width = width(chr_mRNA_rep_notNRZGR)),
                                  strand = strand(chr_mRNA_rep_notNRZGR))
    nonCENranLocGR <- append(nonCENranLocGR, chr_nonCENranLocGR)
  }
}
stopifnot(identical(width(nonCENranLocGR), width(mRNA_rep_notNRZGR)))


CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = as.integer(start(CENranLocGR)-1),
                            end = as.integer(end(CENranLocGR)),
                            name = as.integer(1:length(CENranLocGR)),
                            score = rep("NA", length(CENranLocGR)),
                            strand = as.character(strand(CENranLocGR)))
write.table(CENranLoc_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_NRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

nonCENranLoc_bed <- data.frame(chr = as.character(seqnames(nonCENranLocGR)),
                               start = as.integer(start(nonCENranLocGR)-1),
                               end = as.integer(end(nonCENranLocGR)),
                               name = as.integer(1:length(nonCENranLocGR)),
                               score = rep("NA", length(nonCENranLocGR)),
                               strand = as.character(strand(nonCENranLocGR)))
write.table(nonCENranLoc_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_notNRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
