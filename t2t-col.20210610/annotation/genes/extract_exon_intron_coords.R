#!/applications/R/R-4.0.0/bin/Rscript

# Create tables of exons and introns within representative genes
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_exon_intron_coords.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
library(doParallel)
library(doFuture)
registerDoFuture()
plan(multicore)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14841110,3823792,13597188,4203902,11784131)[which(fai$V1 %in% chrName)]
CENend <- c(17559778,6045243,15733925,6977949,14551809)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of gene coordinates in t2t-col.20210610 (GFF3 derived from Liftoff tool)
mRNA_rep <- readGFF(paste0("t2t-col.20210610_representative_mRNA_", paste0(chrName, collapse = "_"), ".gff3"))
mRNA_rep <- mRNA_rep[mRNA_rep$seqid %in% chrName,]
mRNA_rep <- data.frame(mRNA_rep)
mRNA_rep$group <- as.character(mRNA_rep$group)
print(dim(mRNA_rep))
#[1] 27258     9

# Load genes including exons
genes <- readGFF("t2t-col.20210610.genes.gff3")
genes <- genes[genes$seqid %in% chrName,]

# Get exons for each representative gene model
exons <- genes[genes$type == "exon",]
exons_rep <- exons[which(as.character(exons$Parent) %in% mRNA_rep$group),]
exons_rep <- data.frame(exons_rep)
exons_rep$Parent <- unlist(exons_rep$Parent)
exons_rep$seqid <- as.character(exons_rep$seqid)
# Remove problem exons on Chr4
exons_rep <- exons_rep[-which(exons_rep$seqid == "Chr4" & exons_rep$Parent == "AT2G01021.1_1"),]

exons_rep <- exons_rep[ order(exons_rep$seqid,
                              exons_rep$start,
                              exons_rep$end), ]

exons_rep_gff <- data.frame(exons_rep[,c(1:8, 16)],
                            stringsAsFactors = F)
write.table(exons_rep_gff,
            file = paste0("t2t-col.20210610_representative_exons_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)

exons_rep_bed <- data.frame(chr = as.character(exons_rep_gff[,1]),
                            start = as.integer(exons_rep_gff[,4]-1),
                            end = as.integer(exons_rep_gff[,5]),
                            name = as.character(exons_rep_gff[,9]),
                            score = as.numeric(exons_rep_gff[,6]),
                            strand = as.character(exons_rep_gff[,7]))
write.table(exons_rep_bed,
            file = paste0("t2t-col.20210610_representative_exons_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Get introns for each representative gene model
parentIDs <- as.character(mRNA_rep$group)
introns_rep_gff <- foreach(x = iter(parentIDs),
                           .combine = "rbind",
                           .multicombine = T,
                           .maxcombine = length(parentIDs)+1e1,
                           .inorder = T,
                           .errorhandling = "pass") %dopar% {

  exonParent <- exons_rep_gff[which(as.character(exons_rep_gff$Parent) == x),]
  if(nrow(exonParent) > 1) {
    intronParent <- data.frame()
    for(i in 1:(nrow(exonParent)-1)) {
      intron_i_start <- exonParent[i,]$end+1
      intron_i_end <- exonParent[i+1,]$start-1
      intron_i <- exonParent[i,]
      intron_i$type <- "intron"
      intron_i$start <- intron_i_start
      intron_i$end <- intron_i_end
      intronParent <- rbind(intronParent, intron_i)
    }
    intronParent
  }

}
write.table(introns_rep_gff,
            file = paste0("t2t-col.20210610_representative_introns_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)

introns_rep_bed <- data.frame(chr = as.character(introns_rep_gff[,1]),
                              start = as.integer(introns_rep_gff[,4]-1),
                              end = as.integer(introns_rep_gff[,5]),
                              name = as.character(introns_rep_gff[,9]),
                              score = as.numeric(introns_rep_gff[,6]),
                              strand = as.character(introns_rep_gff[,7]))
write.table(introns_rep_bed,
            file = paste0("t2t-col.20210610_representative_introns_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
