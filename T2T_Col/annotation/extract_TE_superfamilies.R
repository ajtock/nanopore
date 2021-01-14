#!/applications/R/R-3.5.0/bin/Rscript

# Create tables of elements within each transposable element superfamily
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_TE_superfamilies.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(doParallel)
library(data.table)

outDir <- "TEs/"
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of TE coordinates in TAIR10 (including Buisine superfamily annotations)
TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
superfamLevels <- sort(unique(TEs$Transposon_Super_Family))
print(superfamLevels)
# [1] "DNA"           "DNA/En-Spm"    "DNA/Harbinger" "DNA/HAT"
# [5] "DNA/Mariner"   "DNA/MuDR"      "DNA/Pogo"      "DNA/Tc1"
# [9] "LINE?"         "LINE/L1"       "LTR/Copia"     "LTR/Gypsy"
#[13] "RathE1_cons"   "RathE2_cons"   "RathE3_cons"   "RC/Helitron"
#[17] "SINE"          "Unassigned"

DNA <- TEs[grep("DNA", TEs$Transposon_Super_Family),]
DNA <- rbind(DNA,
             TEs[grep("Helitron", TEs$Transposon_Super_Family),])
DNA <- data.frame(DNA,
                  DNA_RNA = "DNA")

Unclassified <- TEs[TEs$Transposon_Super_Family == "Unassigned",]
Unclassified <- data.frame(Unclassified,
                           DNA_RNA = "Unclassified")

RNA <- TEs[!(TEs$Transposon_Name %in% DNA$Transposon_Name)
         & !(TEs$Transposon_Name %in% Unclassified$Transposon_Name),]
RNA <- data.frame(RNA,
                  DNA_RNA = "RNA")

TEs <- rbind(DNA, RNA, Unclassified)

TEs <- rbind(
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA/En-Spm",],
             superfamName = "EnSpm"),
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA/Harbinger",],
             superfamName = "Harbinger"),
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA/HAT",],
             superfamName = "hAT"),
  data.frame(TEs[TEs$Transposon_Super_Family == "RC/Helitron",],
             superfamName = "Helitron"),
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA/MuDR",],
             superfamName = "MuDR"),
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA/Pogo"
               | TEs$Transposon_Super_Family == "DNA/Tc1"
               | TEs$Transposon_Super_Family == "DNA/Mariner",],
             superfamName = "Pogo_Tc1_Mariner"),
  data.frame(TEs[TEs$Transposon_Super_Family == "DNA",],
             superfamName = "Unclassified_DNA"),
  data.frame(TEs[TEs$Transposon_Super_Family == "LTR/Copia",],
             superfamName = "Copia_LTR"),
  data.frame(TEs[TEs$Transposon_Super_Family == "LTR/Gypsy",],
             superfamName = "Gypsy_LTR"),
  data.frame(TEs[TEs$Transposon_Super_Family == "LINE/L1",],
             superfamName = "LINE1"),
  data.frame(TEs[TEs$Transposon_Super_Family == "LINE?",],
             superfamName = "Putative_LINE"),
  data.frame(TEs[TEs$Transposon_Super_Family == "SINE"
               | TEs$Transposon_Super_Family == "RathE1_cons"
               | TEs$Transposon_Super_Family == "RathE2_cons"
               | TEs$Transposon_Super_Family == "RathE3_cons",],
             superfamName = "SINE"),
  data.frame(TEs[TEs$Transposon_Super_Family == "Unassigned",],
             superfamName = "Unclassified")
)

superfamNames <- sort(unique(TEs$superfamName))

# Create a list of TE superfamily DataFrames
superfamList <- lapply(seq_along(superfamNames), function(x) {
  TEs[TEs$superfamName == superfamNames[x],]
})
TEsDNA <- TEs[TEs$DNA_RNA == "DNA",]
TEsRNA <- TEs[TEs$DNA_RNA == "RNA",]

# Load table of TE coordinates in T2T_Col (derived from Liftoff tool)
TEsT2T <- read.csv("TE_liftover_TAIR10_T2T.csv",
                   header = T)
TEsT2T <- data.frame(chr = TEsT2T$chr,
                     start = TEsT2T$start,
                     end = TEsT2T$stop,
                     strand = TEsT2T$strand,
                     name = TEsT2T$atg.gene)
TEsT2T <- TEsT2T[TEsT2T$chr %in% chrName,]
print(dim(TEsT2T))
#[1] 27370     5

TEsT2T_list <- lapply(seq_along(superfamNames), function(x) {
  data.frame(TEsT2T[TEsT2T$name %in%
                    TEs[TEs$superfamName == superfamNames[x],]$Transposon_Name,],
             DNA_RNA = TEs[TEs$superfamName == superfamNames[x],]$DNA_RNA[1],
             superfamName = superfamNames[x])
})
TEsT2T <- do.call(rbind, TEsT2T_list)

# Order by chromosome, start coordinates and end coordinates
TEsT2T <- TEsT2T[ order(TEsT2T$chr,
                        TEsT2T$start,
                        TEsT2T$end), ]

write.table(TEsT2T,
            file = paste0(outDir, "T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
TEsT2T_bed <- data.frame(chr = as.character(TEsT2T[,1]),
                         start = as.integer(TEsT2T[,2]-1),
                         end = as.integer(TEsT2T[,3]),
                         name = as.integer(1:dim(TEsT2T)[1]),
                         score = as.character(TEsT2T[,5]),
                         strand = as.character(TEsT2T[,4]))
write.table(TEsT2T_bed,
            file = paste0(outDir, "T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

TEsT2TGR <- GRanges(seqnames = TEsT2T$chr,
                    ranges = IRanges(start = TEsT2T$start,
                                     end = TEsT2T$end),
                    strand = TEsT2T$strand)

# Define function to select randomly positioned loci of the same
# width distribution as TEsT2TGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as TEsT2TGR
chrs <- chrs[chrs %in% chrName]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  TEsT2TChrGR <- TEsT2TGR[seqnames(TEsT2TGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(TEsT2TChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(TEsT2TChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(TEsT2TChrGR)),
                         strand = strand(TEsT2TChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0(outDir, "T2T_Col_TEs_All_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Extract elements in each superfamily and define matched random loci
registerDoParallel(cores = length(superfamNames))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(h = seq_along(superfamNames)) %dopar% {
#for(h in seq_along(superfamNames)) {
  print(superfamNames[h])
  TEsT2Tsf <- TEsT2T[TEsT2T$superfamName == superfamNames[h],]

  # Order by chromosome, start coordinates and end coordinates
  TEsT2Tsf <- TEsT2Tsf[ order(TEsT2Tsf$chr,
                              TEsT2Tsf$start,
                              TEsT2Tsf$end), ]
  write.table(TEsT2Tsf,
              file = paste0(outDir, "T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  TEsT2Tsf_bed <- data.frame(chr = as.character(TEsT2Tsf[,1]),
                             start = as.integer(TEsT2Tsf[,2]-1),
                             end = as.integer(TEsT2Tsf[,3]),
                             name = as.integer(1:dim(TEsT2Tsf)[1]),
                             score = as.character(TEsT2Tsf[,5]),
                             strand = as.character(TEsT2Tsf[,4]))
  write.table(TEsT2Tsf_bed,
              file = paste0(outDir, "T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  TEsT2TsfGR <- GRanges(seqnames = TEsT2Tsf$chr,
                        ranges = IRanges(start = TEsT2Tsf$start,
                                         end = TEsT2Tsf$end),
                        strand = TEsT2Tsf$strand)
  
  # Define function to select randomly positioned loci of the same
  # width distribution as TEsT2TsfGR
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as TEsT2TsfGR
  chrs <- chrs[chrs %in% chrName]
  ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    TEsT2TsfChrGR <- TEsT2TsfGR[seqnames(TEsT2TsfGR) == chrs[i]]
    regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
    # Contract regionChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(regionChrGR) <- end(regionChrGR)-max(width(TEsT2TsfChrGR))-2000
    start(regionChrGR) <- start(regionChrGR)+2000
    ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                               start(regionChrGR[x]) : end(regionChrGR[x])          
                                                             })),
                                        n = length(TEsT2TsfChrGR))
    ranLocChrGR <- GRanges(seqnames = chrs[i],
                           ranges = IRanges(start = ranLocChrStart,
                                            width = width(TEsT2TsfChrGR)),
                           strand = strand(TEsT2TsfChrGR))
    ranLocGR <- append(ranLocGR, ranLocChrGR)
  }
  
  ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                           start = as.integer(start(ranLocGR)-1),
                           end = as.integer(end(ranLocGR)),
                           name = as.integer(1:length(ranLocGR)),
                           score = rep("NA", length(ranLocGR)),
                           strand = as.character(strand(ranLocGR)))
  write.table(ranLoc_bed,
              file = paste0(outDir, "T2T_Col_TEs_", superfamNames[h], "_",
                            paste0(chrName, collapse = "_"), "_randomLoci.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
