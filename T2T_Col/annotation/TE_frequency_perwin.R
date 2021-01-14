#!/applications/R/R-3.5.0/bin/Rscript

###################################################################################
# Generate genome profiles of feature frequency in adjacent windows of given size # 
###################################################################################

# Usage:
# ./TE_frequency_perwin.R 10000 10kb

args <- commandArgs(trailingOnly = T)
winSize <- as.numeric(args[1])
winName <- as.character(args[2])

library(GenomicRanges)
library(parallel)
library(doParallel)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)

# Count number of features in each adjacent genomic windows using countOverlaps()
TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
superfamLevels <- levels(TEs$Transposon_Super_Family)
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

superfamNames <- levels(TEs$superfamName)

# Create a list of TE superfamily DataFrames
superfamList <- lapply(seq_along(superfamNames), function(x) {
  TEs[TEs$superfamName == superfamNames[x],]
})
TEsDNA <- TEs[TEs$DNA_RNA == "DNA",]
TEsRNA <- TEs[TEs$DNA_RNA == "RNA",]

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = winSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = winSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

registerDoParallel(cores = length(superfamList))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(h = 1:length(superfamList)) %dopar% {
  superfamProfile <- NULL
  for(i in 1:length(chrs)) {
    # Count superfamily TEs within windows
    chrTEs <- superfamList[[h]][superfamList[[h]]$Chr == chrs[i],]
    chrTEsGR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chrTEs$start,
                                         end = chrTEs$end),
                        strand = chrTEs$strand)
    chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
    winTEs <- countOverlaps(chrWindowsGR,
                            chrTEsGR,
                            ignore.strand = T)
    chrProfile <- data.frame(chr = as.character(chrs[i]),
                             window = as.integer(start(chrWindowsGR)),
                             winTEs = as.integer(winTEs))
    superfamProfile <- rbind(superfamProfile, chrProfile)
  }
  write.table(superfamProfile,
              file = paste0("TE_frequency_per_", winName,
                            "_superfamily_", superfamNames[h],
                            ".txt"))
}

# DNA TEs
TEsDNAprofile <- NULL
for(i in 1:length(chrs)) {
  # Count superfamily TEs within windows
  chrTEs <- TEsDNA[TEsDNA$Chr == chrs[i],]
  chrTEsGR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = chrTEs$start,
                                       end = chrTEs$end),
                      strand = chrTEs$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  winTEs <- countOverlaps(chrWindowsGR,
                          chrTEsGR,
                          ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           winTEs = as.integer(winTEs))
  TEsDNAprofile <- rbind(TEsDNAprofile, chrProfile)
}
write.table(TEsDNAprofile,
            file = paste0("TE_frequency_per_", winName,
                          "_DNA.txt"))

# RNA TEs
TEsRNAprofile <- NULL
for(i in 1:length(chrs)) {
  # Count superfamily TEs within windows
  chrTEs <- TEsRNA[TEsRNA$Chr == chrs[i],]
  chrTEsGR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = chrTEs$start,
                                       end = chrTEs$end),
                      strand = chrTEs$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  winTEs <- countOverlaps(chrWindowsGR,
                          chrTEsGR,
                          ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           winTEs = as.integer(winTEs))
  TEsRNAprofile <- rbind(TEsRNAprofile, chrProfile)
}
write.table(TEsRNAprofile,
            file = paste0("TE_frequency_per_", winName,
                          "_RNA.txt"))

# TEs
TEsProfile <- NULL
for(i in 1:length(chrs)) {
  # Count superfamily TEs within windows
  chrTEs <- TEs[TEs$Chr == chrs[i],]
  chrTEsGR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = chrTEs$start,
                                       end = chrTEs$end),
                      strand = chrTEs$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  winTEs <- countOverlaps(chrWindowsGR,
                          chrTEsGR,
                          ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           winTEs = as.integer(winTEs))
  TEsProfile <- rbind(TEsProfile, chrProfile)
}
write.table(TEsProfile,
            file = paste0("TE_frequency_per_", winName,
                          ".txt"))
