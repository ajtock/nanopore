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

TEs <- TEs[TEs$seqid %in% chrName,]
print(dim(TEs))
#[1] 42770    22

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

# Order by chromosome, start coordinates and end coordinates
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


# Define function to select randomly positioned loci of the same
# width distribution as nonCEN_TEsGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as TEsGR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  nonCEN_TEsChrGR <- nonCEN_TEsGR[seqnames(nonCEN_TEsGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(nonCEN_TEsChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(nonCEN_TEsChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(nonCEN_TEsChrGR)),
                         strand = strand(nonCEN_TEsChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0("T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


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

  # Order by chromosome, start coordinates and end coordinates
  TEssf <- TEssf[ order(TEssf$seqid,
                        TEssf$start,
                        TEssf$end), ]
  write.table(TEssf,
              file = paste0("T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  TEssf_bed <- data.frame(chr = as.character(TEssf[,1]),
                          start = as.integer(TEssf[,4]-1),
                          end = as.integer(TEssf[,5]),
                          name = paste0(1:dim(TEssf)[1], "_", TEssf[,10]),
                          score = as.character(TEssf[,6]),
                          strand = as.character(TEssf[,7]))
  write.table(TEssf_bed,
              file = paste0("T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  TEssfGR <- GRanges(seqnames = TEssf$seqid,
                     ranges = IRanges(start = TEssf$start,
                                      end = TEssf$end),
                     strand = TEssf$strand)
  
  # Define function to select randomly positioned loci of the same
  # width distribution as TEssfGR
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as TEssfGR
  ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    TEssfChrGR <- TEssfGR[seqnames(TEssfGR) == chrs[i]]
    regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
    # Contract regionChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(regionChrGR) <- end(regionChrGR)-max(width(TEssfChrGR))-2000
    start(regionChrGR) <- start(regionChrGR)+2000
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                               start(regionChrGR[x]) : end(regionChrGR[x])          
                                                             })),
                                        n = length(TEssfChrGR))
    ranLocChrGR <- GRanges(seqnames = chrs[i],
                           ranges = IRanges(start = ranLocChrStart,
                                            width = width(TEssfChrGR)),
                           strand = strand(TEssfChrGR))
    ranLocGR <- append(ranLocGR, ranLocChrGR)
  }
  
  ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                           start = as.integer(start(ranLocGR)-1),
                           end = as.integer(end(ranLocGR)),
                           name = as.integer(1:length(ranLocGR)),
                           score = rep("NA", length(ranLocGR)),
                           strand = as.character(strand(ranLocGR)))
  write.table(ranLoc_bed,
              file = paste0("T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), "_randomLoci.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

superfamDNARNAUnclass <- sort(unique(TEs$DNA_RNA))
foreach(h = seq_along(superfamDNARNAUnclass)) %dopar% {
#for(h in seq_along(superfamDNARNAUnclass)) {
  print(superfamDNARNAUnclass[h])
  TEssf <- TEs[TEs$DNA_RNA == superfamDNARNAUnclass[h],]

  # Order by chromosome, start coordinates and end coordinates
  TEssf <- TEssf[ order(TEssf$seqid,
                        TEssf$start,
                        TEssf$end), ]
  write.table(TEssf,
              file = paste0("T2T_Col_TEs_", superfamDNARNAUnclass[h], "_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  TEssf_bed <- data.frame(chr = as.character(TEssf[,1]),
                          start = as.integer(TEssf[,4]-1),
                          end = as.integer(TEssf[,5]),
                          name = paste0(1:dim(TEssf)[1], "_", TEssf[,10]),
                          score = as.character(TEssf[,6]),
                          strand = as.character(TEssf[,7]))
  write.table(TEssf_bed,
              file = paste0("T2T_Col_TEs_", superfamDNARNAUnclass[h], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  TEssfGR <- GRanges(seqnames = TEssf$seqid,
                     ranges = IRanges(start = TEssf$start,
                                      end = TEssf$end),
                     strand = TEssf$strand)
  
  # Define function to select randomly positioned loci of the same
  # width distribution as TEssfGR
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as TEssfGR
  ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    TEssfChrGR <- TEssfGR[seqnames(TEssfGR) == chrs[i]]
    regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
    # Contract regionChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(regionChrGR) <- end(regionChrGR)-max(width(TEssfChrGR))-2000
    start(regionChrGR) <- start(regionChrGR)+2000
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                               start(regionChrGR[x]) : end(regionChrGR[x])          
                                                             })),
                                        n = length(TEssfChrGR))
    ranLocChrGR <- GRanges(seqnames = chrs[i],
                           ranges = IRanges(start = ranLocChrStart,
                                            width = width(TEssfChrGR)),
                           strand = strand(TEssfChrGR))
    ranLocGR <- append(ranLocGR, ranLocChrGR)
  }
  
  ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                           start = as.integer(start(ranLocGR)-1),
                           end = as.integer(end(ranLocGR)),
                           name = as.integer(1:length(ranLocGR)),
                           score = rep("NA", length(ranLocGR)),
                           strand = as.character(strand(ranLocGR)))
  write.table(ranLoc_bed,
              file = paste0("T2T_Col_TEs_", superfamDNARNAUnclass[h], "_",
                            paste0(chrName, collapse = "_"), "_randomLoci.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
