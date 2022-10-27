#!/applications/R/R-4.0.0/bin/Rscript

# Create tables of LRZ genes, NRZ genes, and arm genes,
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_gene_coords_LRZgenes_NRZgenes_v271022.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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

# Region definitions
regions <- read.csv("LRZ_NRZ/Table_S1_LRZ_NRZ.csv", header = T)[,-1]
regions <- regions[ with(regions, order(Chr, decreasing = F)), ]

LRZstart <- regions$LRZ.left
LRZend <- regions$LRZ.right
LRZGR <- GRanges(seqnames = regions$Chr,
                 ranges = IRanges(start = LRZstart,
                                  end = LRZend),
                 strand = "*")

NRZstart <- regions$NRZ.left
NRZend <- regions$NRZ.right
NRZGR <- GRanges(seqnames = regions$Chr,
                 ranges = IRanges(start = NRZstart,
                                  end = NRZend),
                 strand = "*")

armGR <- GRanges(seqnames = rep(seqnames(LRZGR), 2),
                 ranges = IRanges(start = c(rep(1, length(start(LRZGR))),
                                            end(LRZGR)+1),
                                  end = c(start(LRZGR)-1,
                                          chrLens)),
                 strand = "*")

# Load table of gene coordinates in t2t-col.20210610 (GFF3 derived from Liftoff tool)
genes <- readGFF("t2t-col.20210610.genes.gff3")
genes <- genes[genes$seqid %in% chrName,]
genes <- genes[genes$type == "gene",]
print(dim(genes))
#[1] 28504    20

# Obtain frequency of occurrence of each gene parent ID
n_occur <- data.frame(table(unlist(genes$Name)))

# Obtain gene records for which the gene parent ID occurs only once
genes_unique <-  as.data.frame(
  genes[ unlist(genes$Name)
    %in% n_occur$Var1[n_occur$Freq == 1],
  ]
)
# Obtain gene records for which the gene parent ID occurs more than once
genes_multi <- as.data.frame(
  genes[ unlist(genes$Name)
    %in% n_occur$Var1[n_occur$Freq > 1],
  ]
)

# For each gene parent ID in genes_multi, obtain the gene record with the
# longest transcript
# If multiple genes records have the longest transcript,
# keep the first reported one only
genes_multi_list <- mclapply(seq_along(genes_multi[,1]), function(h) {
  genes_multi_ID_all <- genes_multi[ unlist(genes_multi$Name)
                          == unlist(genes_multi[h,]$Name),
                        ]
  genes_multi_ID_all[ genes_multi_ID_all$end-genes_multi_ID_all$start
    == max(genes_multi_ID_all$end-genes_multi_ID_all$start),
  ][1,]
}, mc.cores = detectCores())

# Collapse genes_multi_list into single data.frame and remove duplicates
genes_multi_dup <- rbindlist(genes_multi_list)
genes_multi_rep <- unique(as.data.frame(genes_multi_dup))

# Combine into one representative set of gene entries, order,
# and output in GFF3 and BED formats
genes_rep <- rbind(genes_unique, genes_multi_rep)
genes_rep <- genes_rep[ order(genes_rep$seqid,
                              genes_rep$start,
                              genes_rep$end), ]
nrow(genes_rep)

# Get genes in low-recombining zones (LRZgenes)
LRZgenes <- read.table("LRZ_NRZ/all.lrz.genes.txt", header = T)[,1]
length(LRZgenes)
#[1] 542
genes_rep_LRZ <- genes_rep[ which(unlist(genes_rep$Name) %in% LRZgenes), ]
length(unique(unlist(genes_rep_LRZ$Name)))
#[1] 542
nrow(genes_rep_LRZ)
#[1] 542
write.table(genes_rep_LRZ[,c(1:8, 11)],
            file = paste0("t2t-col.20210610_representative_genes_LRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_LRZ_bed <- data.frame(chr = as.character(genes_rep_LRZ[,1]),
                                start = as.integer(genes_rep_LRZ[,4]-1),
                                end = as.integer(genes_rep_LRZ[,5]),
                                name = as.character(genes_rep_LRZ[,11]),
                                score = as.numeric(genes_rep_LRZ[,6]),
                                strand = as.character(genes_rep_LRZ[,7]))
write.table(genes_rep_LRZ_bed,
            file = paste0("t2t-col.20210610_representative_genes_LRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_LRZGR <- GRanges(seqnames = genes_rep_LRZ$seqid,
                           ranges = IRanges(start = genes_rep_LRZ$start,
                                            end = genes_rep_LRZ$end),
                           strand = genes_rep_LRZ$strand)

# Get genes in non-recombining zones (NRZgenes)
NRZgenes <- read.table("LRZ_NRZ/all.nrz.genes.txt", header = T)[,1]
length(NRZgenes)
#[1] 132
genes_rep_NRZ <- genes_rep[ which(unlist(genes_rep$Name) %in% NRZgenes), ]
length(unique(unlist(genes_rep_NRZ$Name)))
#[1] 132
nrow(genes_rep_NRZ)
#[1] 132
write.table(genes_rep_NRZ[,c(1:8, 11)],
            file = paste0("t2t-col.20210610_representative_genes_NRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_NRZ_bed <- data.frame(chr = as.character(genes_rep_NRZ[,1]),
                                start = as.integer(genes_rep_NRZ[,4]-1),
                                end = as.integer(genes_rep_NRZ[,5]),
                                name = as.character(genes_rep_NRZ[,11]),
                                score = as.numeric(genes_rep_NRZ[,6]),
                                strand = as.character(genes_rep_NRZ[,7]))
write.table(genes_rep_NRZ_bed,
            file = paste0("t2t-col.20210610_representative_genes_NRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_NRZGR <- GRanges(seqnames = genes_rep_NRZ$seqid,
                           ranges = IRanges(start = genes_rep_NRZ$start,
                                            end = genes_rep_NRZ$end),
                           strand = genes_rep_NRZ$strand)

# Get genes in arms (armgenes)
armgenes <- read.table("LRZ_NRZ/all.arm.genes.txt", header = T)[,1]
length(armgenes)
#[1] 27567
armgenes <- unique(armgenes)
length(armgenes)
#[1] 27493
genes_rep_arm <- genes_rep[ which(unlist(genes_rep$Name) %in% armgenes), ]
length(unique(unlist(genes_rep_arm$Name)))
#[1] 27493
nrow(genes_rep_arm)
#[1] 27493
write.table(genes_rep_arm[,c(1:8, 11)],
            file = paste0("t2t-col.20210610_representative_genes_arm_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_arm_bed <- data.frame(chr = as.character(genes_rep_arm[,1]),
                                start = as.integer(genes_rep_arm[,4]-1),
                                end = as.integer(genes_rep_arm[,5]),
                                name = as.character(genes_rep_arm[,11]),
                                score = as.numeric(genes_rep_arm[,6]),
                                strand = as.character(genes_rep_arm[,7]))
write.table(genes_rep_arm_bed,
            file = paste0("t2t-col.20210610_representative_genes_arm_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_armGR <- GRanges(seqnames = genes_rep_arm$seqid,
                           ranges = IRanges(start = genes_rep_arm$start,
                                            end = genes_rep_arm$end),
                           strand = genes_rep_arm$strand)


# Concatenate and write LRZgenes and NRZgenes
genes_rep_LRZ_NRZ <- rbind(genes_rep_LRZ, genes_rep_NRZ)
length(unique(unlist(genes_rep_LRZ_NRZ$Name)))
#[1] 674
nrow(genes_rep_LRZ_NRZ)
#[1] 674
write.table(genes_rep_LRZ_NRZ[,c(1:8, 11)],
            file = paste0("t2t-col.20210610_representative_genes_LRZ_NRZ_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_LRZ_NRZ_bed <- data.frame(chr = as.character(genes_rep_LRZ_NRZ[,1]),
                                    start = as.integer(genes_rep_LRZ_NRZ[,4]-1),
                                    end = as.integer(genes_rep_LRZ_NRZ[,5]),
                                    name = as.character(genes_rep_LRZ_NRZ[,11]),
                                    score = as.numeric(genes_rep_LRZ_NRZ[,6]),
                                    strand = as.character(genes_rep_LRZ_NRZ[,7]))
write.table(genes_rep_LRZ_NRZ_bed,
            file = paste0("t2t-col.20210610_representative_genes_LRZ_NRZ_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
genes_rep_LRZ_NRZGR <- GRanges(seqnames = genes_rep_LRZ_NRZ$seqid,
                               ranges = IRanges(start = genes_rep_LRZ_NRZ$start,
                                                end = genes_rep_LRZ_NRZ$end),
                               strand = genes_rep_LRZ_NRZ$strand)


# Define function to select randomly positioned loci of the same
# width distribution as genes_rep_LRZGR, genes_rep_NRZGR or genes_rep_armGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

chrs <- chrs[chrs %in% chrName]

# Apply ranLocStartSelect() on a per-chromosome basis so that
# LRZranLocGR contains the same number of loci per chromosome as genes_rep_LRZGR
LRZranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_LRZGR <- LRZGR[seqnames(LRZGR) == chrs[j]]

  chr_genes_rep_LRZGR <- genes_rep_LRZGR[seqnames(genes_rep_LRZGR) == chrs[j]]

  if(length(chr_genes_rep_LRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_LRZranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_LRZGR), function(x) {
                                                                    ( start(chr_LRZGR[x]) + max(width(chr_genes_rep_LRZGR)) + 2000 ) :
                                                                    ( end(chr_LRZGR[x]) - max(width(chr_genes_rep_LRZGR)) - 2000 )
                                                                  })),
                                             n = length(chr_genes_rep_LRZGR))
    chr_LRZranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_LRZranLoc_Start,
                                                width = width(chr_genes_rep_LRZGR)),
                               strand = strand(chr_genes_rep_LRZGR))
    LRZranLocGR <- append(LRZranLocGR, chr_LRZranLocGR)
  }
}
stopifnot(identical(width(LRZranLocGR), width(genes_rep_LRZGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# NRZranLocGR contains the same number of loci per chromosome as genes_rep_NRZGR
NRZranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_NRZGR <- NRZGR[seqnames(NRZGR) == chrs[j]]

  chr_genes_rep_NRZGR <- genes_rep_NRZGR[seqnames(genes_rep_NRZGR) == chrs[j]]

  if(length(chr_genes_rep_NRZGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_NRZranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_NRZGR), function(x) {
                                                                    ( start(chr_NRZGR[x]) + max(width(chr_genes_rep_NRZGR)) + 2000 ) :
                                                                    ( end(chr_NRZGR[x]) - max(width(chr_genes_rep_NRZGR)) - 2000 )
                                                                  })),
                                             n = length(chr_genes_rep_NRZGR))
    chr_NRZranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_NRZranLoc_Start,
                                                width = width(chr_genes_rep_NRZGR)),
                               strand = strand(chr_genes_rep_NRZGR))
    NRZranLocGR <- append(NRZranLocGR, chr_NRZranLocGR)
  }
}
stopifnot(identical(width(NRZranLocGR), width(genes_rep_NRZGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# armranLocGR contains the same number of loci per chromosome as genes_rep_armGR
armranLocGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_armGR <- armGR[seqnames(armGR) == chrs[j]]

  chr_genes_rep_armGR <- genes_rep_armGR[seqnames(genes_rep_armGR) == chrs[j]]

  if(length(chr_genes_rep_armGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_armranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_armGR), function(x) {
                                                                    ( start(chr_armGR[x]) + max(width(chr_genes_rep_armGR)) + 2000 ) :
                                                                    ( end(chr_armGR[x]) - max(width(chr_genes_rep_armGR)) - 2000 )
                                                                  })),
                                             n = length(chr_genes_rep_armGR))
    chr_armranLocGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_armranLoc_Start,
                                                width = width(chr_genes_rep_armGR)),
                               strand = strand(chr_genes_rep_armGR))
    armranLocGR <- append(armranLocGR, chr_armranLocGR)
  }
}
stopifnot(identical(width(armranLocGR), width(genes_rep_armGR)))

# Concatenate LRZranLoc and NRZranLoc
LRZ_NRZranLocGR <- c(LRZranLocGR, NRZranLocGR)
stopifnot(identical(width(LRZ_NRZranLocGR), width(genes_rep_LRZ_NRZGR)))

# Write ranLoc to BED
LRZranLoc_bed <- data.frame(chr = as.character(seqnames(LRZranLocGR)),
                            start = as.integer(start(LRZranLocGR)-1),
                            end = as.integer(end(LRZranLocGR)),
                            name = as.integer(1:length(LRZranLocGR)),
                            score = rep("NA", length(LRZranLocGR)),
                            strand = as.character(strand(LRZranLocGR)))
write.table(LRZranLoc_bed,
            file = paste0("t2t-col.20210610_representative_genes_LRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

NRZranLoc_bed <- data.frame(chr = as.character(seqnames(NRZranLocGR)),
                            start = as.integer(start(NRZranLocGR)-1),
                            end = as.integer(end(NRZranLocGR)),
                            name = as.integer(1:length(NRZranLocGR)),
                            score = rep("NA", length(NRZranLocGR)),
                            strand = as.character(strand(NRZranLocGR)))
write.table(NRZranLoc_bed,
            file = paste0("t2t-col.20210610_representative_genes_NRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

armranLoc_bed <- data.frame(chr = as.character(seqnames(armranLocGR)),
                            start = as.integer(start(armranLocGR)-1),
                            end = as.integer(end(armranLocGR)),
                            name = as.integer(1:length(armranLocGR)),
                            score = rep("NA", length(armranLocGR)),
                            strand = as.character(strand(armranLocGR)))
write.table(armranLoc_bed,
            file = paste0("t2t-col.20210610_representative_genes_arm_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

LRZ_NRZranLoc_bed <- data.frame(chr = as.character(seqnames(LRZ_NRZranLocGR)),
                                start = as.integer(start(LRZ_NRZranLocGR)-1),
                                end = as.integer(end(LRZ_NRZranLocGR)),
                                name = as.integer(1:length(LRZ_NRZranLocGR)),
                                score = rep("NA", length(LRZ_NRZranLocGR)),
                                strand = as.character(strand(LRZ_NRZranLocGR)))
write.table(LRZ_NRZranLoc_bed,
            file = paste0("t2t-col.20210610_representative_genes_LRZ_NRZ_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

