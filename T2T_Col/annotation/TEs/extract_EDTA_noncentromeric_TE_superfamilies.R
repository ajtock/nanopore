#!/applications/R/R-4.0.0/bin/Rscript

# Create tables of elements within each transposable element superfamily
# not located within T2T_Col centromeres
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_EDTA_noncentromeric_TE_superfamilies.R 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#superfamName <- "Gypsy_LTR"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
superfamName <- args[2]

options(stringsAsFactors = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(doParallel)
library(data.table)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load TEs
TEs <- read.table(file = paste0("T2T_Col_TEs_All_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                  header = T)

superfamNames <- sort(unique(TEs$superfamName))
print(superfamNames)
# [1] "Copia_LTR"         "EnSpm"             "Gypsy_LTR"
# [4] "Harbinger"         "hAT"               "Helitron"
# [7] "LINE1"             "MuDR"              "Pogo_Tc1_Mariner"
#[10] "SINE"              "Unclassified"      "Unclassified_DNA"
#[13] "Unclassified_LINE" "Unclassified_LTR"

print(dim(TEs))
#[1] 42770    22

# Create a list of TE superfamily DataFrames
superfamList <- lapply(seq_along(superfamNames), function(x) {
  TEs[TEs$superfamName == superfamNames[x],]
})
tmp <- 0
for(x in seq_along(superfamList)) {
  tmp <- tmp + dim(superfamList[[x]])[1]
}
print(tmp)
#[1] 42770

stopifnot(identical(dim(TEs)[1], as.integer(tmp)))

# Remove features shorter than 100 bp
TEs <- TEs[which((TEs$end - TEs$start) + 1 >= 100),]

TEs <- TEs[TEs$seqid %in% chrName,]
print(dim(TEs))
#[1] 39034    22

# Order by chromosome, start coordinates and end coordinates
TEs <- TEs[ order(TEs$seqid,
                  TEs$start,
                  TEs$end), ]

TEsGR <- GRanges(seqnames = TEs$seqid,
                 ranges = IRanges(start = TEs$start,
                                  end = TEs$end),
                 strand = TEs$strand,
                 name = paste0(1:dim(TEs)[1], "_", TEs[,10]),
                 score = as.character(TEs[,6]))

TEs_CEN_overlaps <- findOverlaps(query = CENGR,
                                 subject = TEsGR,
                                 type = "any",
                                 select = "all",
                                 ignore.strand = TRUE)
if(length(TEs_CEN_overlaps) > 0) {
  nonCEN_TEsGR <- TEsGR[-subjectHits(TEs_CEN_overlaps)]
} else {
  nonCEN_TEsGR <- TEsGR
}

# Convert into BED format
nonCEN_TEs_bed <- data.frame(chr = as.character(seqnames(nonCEN_TEsGR)),
                             start = as.integer(start(nonCEN_TEsGR)-1),
                             end = as.integer(end(nonCEN_TEsGR)),
                             name = as.character(nonCEN_TEsGR$name),
                             score = as.character(nonCEN_TEsGR$score),
                             strand = as.character(strand(nonCEN_TEsGR)))
write.table(nonCEN_TEs_bed,
            file = paste0("T2T_Col_nonCEN_TEs_All_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Extract noncentromeric TEs for each superfamily or a given superfamily

if(superfamName != "All") {
  superfamNames <- superfamName
}

# Extract elements in each superfamily and define matched random loci
registerDoParallel(cores = length(superfamNames))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(h = seq_along(superfamNames)) %dopar% {
#for(h in seq_along(superfamNames)) {
  print(superfamNames[h])
  TEssf <- TEs[TEs$superfamName == superfamNames[h],]

  # Remove features shorter than 100 bp
  TEssf <- TEssf[which((TEssf$end - TEssf$start) + 1 >= 100),]

  # Order by chromosome, start coordinates and end coordinates
  TEssf <- TEssf[ order(TEssf$seqid,
                        TEssf$start,
                        TEssf$end), ]
  
  TEssfGR <- GRanges(seqnames = TEssf$seqid,
                     ranges = IRanges(start = TEssf$start,
                                      end = TEssf$end),
                     strand = TEssf$strand,
                     name = paste0(1:dim(TEssf)[1], "_", TEssf[,10]),
                     score = as.character(TEssf[,6]))
  
  TEssf_CEN_overlaps <- findOverlaps(query = CENGR,
                                     subject = TEssfGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
  if(length(TEssf_CEN_overlaps) > 0) {
    nonCEN_TEssfGR <- TEssfGR[-subjectHits(TEssf_CEN_overlaps)]
  } else {
    nonCEN_TEssfGR <- TEssfGR
  }
  
  # Convert into BED format
  nonCEN_TEssf_bed <- data.frame(chr = as.character(seqnames(nonCEN_TEssfGR)),
                                 start = as.integer(start(nonCEN_TEssfGR)-1),
                                 end = as.integer(end(nonCEN_TEssfGR)),
                                 name = as.character(nonCEN_TEssfGR$name),
                                 score = as.character(nonCEN_TEssfGR$score),
                                 strand = as.character(strand(nonCEN_TEssfGR)))
  write.table(nonCEN_TEssf_bed,
              file = paste0("T2T_Col_nonCEN_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
