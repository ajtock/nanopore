#!/applications/R/R-4.0.0/bin/Rscript

# Create tables of elements within each transposable element superfamily
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_EDTA_TE_superfamilies.R 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR

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

# Load table of TE coordinates in T2T_Col (derived from EDTA)
TEs <- data.frame(readGFF("t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3"))
superfamNamesFull <- sort(unique(TEs$Classification))
print(superfamNamesFull)
# [1] "DNA"           "DNA/DTA"       "DNA/DTC"       "DNA/DTH"      
# [5] "DNA/DTM"       "DNA/DTT"       "DNA/En-Spm"    "DNA/Harbinger"
# [9] "DNA/HAT"       "DNA/Helitron"  "DNA/Mariner"   "DNA/MuDR"     
#[13] "DNA/Pogo"      "DNA/Tc1"       "LINE?"         "LINE/L1"      
#[17] "LINE/unknown"  "LTR/Copia"     "LTR/Gypsy"     "LTR/unknown"  
#[21] "MITE/DTA"      "MITE/DTC"      "MITE/DTH"      "MITE/DTM"     
#[25] "MITE/DTT"      "RathE1_cons"   "RathE2_cons"   "RathE3_cons"  
#[29] "RC/Helitron"   "SINE"          "Unassigned"    "Unknown"   

# Get DNA elements
DNA <- TEs[grep("DNA", TEs$Classification),]
DNA <- rbind(DNA,
             TEs[grep("RC/Helitron", TEs$Classification),],
             TEs[grep("MITE", TEs$Classification),]
            )
DNA <- data.frame(DNA,
                  DNA_RNA = "DNA")

# Get unclassified elements
Unclassified <- TEs[TEs$Classification %in% c("Unassigned", "Unknown"),]
Unclassified <- data.frame(Unclassified,
                           DNA_RNA = "Unclassified")

# Get RNA elements
RNA <- TEs[!(TEs$Classification %in% DNA$Classification)
         & !(TEs$Classification %in% Unclassified$Classification),]
RNA <- data.frame(RNA,
                  DNA_RNA = "RNA")

# Combine
TEs <- rbind(DNA, RNA, Unclassified)

# Make superfamily names filename-friendly
TEs <- rbind(
  data.frame(TEs[TEs$Classification == "DNA/En-Spm"
               | TEs$Classification == "DNA/DTC"
               | TEs$Classification == "MITE/DTC",],
             superfamName = "EnSpm"),
  data.frame(TEs[TEs$Classification == "DNA/Harbinger"
               | TEs$Classification == "DNA/DTH"
               | TEs$Classification == "MITE/DTH",],
             superfamName = "Harbinger"),
  data.frame(TEs[TEs$Classification == "DNA/HAT"
               | TEs$Classification == "DNA/DTA"
               | TEs$Classification == "MITE/DTA",],
             superfamName = "hAT"),
  data.frame(TEs[TEs$Classification == "DNA/Helitron"
               | TEs$Classification == "RC/Helitron",],
             superfamName = "Helitron"),
  data.frame(TEs[TEs$Classification == "DNA/MuDR"
               | TEs$Classification == "DNA/DTM"
               | TEs$Classification == "MITE/DTM",],
             superfamName = "MuDR"),
  data.frame(TEs[TEs$Classification == "DNA/Pogo"
               | TEs$Classification == "DNA/Tc1"
               | TEs$Classification == "DNA/Mariner"
               | TEs$Classification == "DNA/DTT"
               | TEs$Classification == "MITE/DTT",],
             superfamName = "Pogo_Tc1_Mariner"),
  data.frame(TEs[TEs$Classification == "DNA",],
             superfamName = "Unclassified_DNA"),
  data.frame(TEs[TEs$Classification == "LTR/Copia",],
             superfamName = "Copia_LTR"),
  data.frame(TEs[TEs$Classification == "LTR/Gypsy",],
             superfamName = "Gypsy_LTR"),
  data.frame(TEs[TEs$Classification == "LTR/unknown",],
             superfamName = "Unclassified_LTR"),
  data.frame(TEs[TEs$Classification == "LINE/L1",],
             superfamName = "LINE1"),
  data.frame(TEs[TEs$Classification == "LINE?"
               | TEs$Classification == "LINE/unknown",],
             superfamName = "Unclassified_LINE"),
  data.frame(TEs[TEs$Classification == "SINE"
               | TEs$Classification == "RathE1_cons"
               | TEs$Classification == "RathE2_cons"
               | TEs$Classification == "RathE3_cons",],
             superfamName = "SINE"),
  data.frame(TEs[TEs$Classification == "Unassigned"
               | TEs$Classification == "Unknown",],
             superfamName = "Unclassified")
)

superfamNames <- sort(unique(TEs$superfamName))
print(superfamNames)
# [1] "Copia_LTR"         "EnSpm"             "Gypsy_LTR"
# [4] "Harbinger"         "hAT"               "Helitron"
# [7] "LINE1"             "MuDR"              "Pogo_Tc1_Mariner"
#[10] "SINE"              "Unclassified"      "Unclassified_DNA"
#[13] "Unclassified_LINE" "Unclassified_LTR"

print(dim(TEs))
#[1] 42927    22

# Create a list of TE superfamily DataFrames
superfamList <- lapply(seq_along(superfamNames), function(x) {
  TEs[TEs$superfamName == superfamNames[x],]
})
tmp <- 0
for(x in seq_along(superfamList)) {
  tmp <- tmp + dim(superfamList[[x]])[1]
}
print(tmp)
#[1] 42927

stopifnot(identical(dim(TEs)[1], as.integer(tmp)))

TEs <- TEs[TEs$seqid %in% chrName,]
print(dim(TEs))
#[1] 42770    22

TEs$Parent <- "."

# Order by chromosome, start coordinates and end coordinates
TEs <- TEs[ order(TEs$seqid,
                  TEs$start,
                  TEs$end), ]
#TEsDNA <- TEs[TEs$DNA_RNA == "DNA",]
#TEsRNA <- TEs[TEs$DNA_RNA == "RNA",]
write.table(TEs,
            file = paste0("T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
TEs_bed <- data.frame(chr = as.character(TEs[,1]),
                      start = as.integer(TEs[,4]-1),
                      end = as.integer(TEs[,5]),
                      name = paste0(1:dim(TEs)[1], "_", TEs[,10]),
                      score = as.character(TEs[,6]),
                      strand = as.character(TEs[,7]))
write.table(TEs_bed,
            file = paste0("T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

TEsGR <- GRanges(seqnames = TEs$seqid,
                 ranges = IRanges(start = TEs$start,
                                  end = TEs$end),
                 strand = TEs$strand)

# Define function to select randomly positioned loci of the same
# width distribution as TEsGR
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
  TEsChrGR <- TEsGR[seqnames(TEsGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(TEsChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(TEsChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(TEsChrGR)),
                         strand = strand(TEsChrGR))
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
