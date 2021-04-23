#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./Athila_centrophiles_TABtoBED_AlexBousios_strand_info_v280321.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

AthilaDir <- "CENAthila/"
soloLTRDir <- "CENsoloLTR/" 
system(paste0("[ -d ", AthilaDir, " ] || mkdir -p ", AthilaDir))
system(paste0("[ -d ", soloLTRDir, " ] || mkdir -p ", soloLTRDir))

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of centromeric gap and Athila sequence coordinates
tab <- read.csv("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/updated_Athila_Andy_v280321.csv",
                header = T)
colnames(tab) <- c("chr", "class", "start", "end", "length_bp", "strand",
                   "left_TSD", "right_TSD", "5prime_LTR", "3prime_LTR",
                   "percent_identity", "5prime_LTR_PBS", "PPT_3prime_LTR", "type")
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
Athila <- tab[tab$type == "full length",]
soloLTR <- tab[tab$type == "soloLTR",]

# Convert Athila into GRanges and then BED format
AthilaGR <- GRanges(seqnames = Athila$chr,
                    ranges = IRanges(start = Athila$start,
                                     end = Athila$end),
                    strand = Athila$strand,
                    class = Athila$class)
AthilaGR <- AthilaGR[seqnames(AthilaGR) %in% chrName]
Athila_bed <- data.frame(chr = as.character(seqnames(AthilaGR)),
                         start = as.integer(start(AthilaGR)-1),
                         end = as.integer(end(AthilaGR)),
                         name = as.integer(1:length(AthilaGR)),
                         score = as.character(AthilaGR$class),
                         strand = as.character(strand(AthilaGR)))
write.table(Athila_bed,
            file = paste0(AthilaDir, "CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Convert soloLTR into GRanges and then BED format
soloLTRGR <- GRanges(seqnames = soloLTR$chr,
                     ranges = IRanges(start = soloLTR$start,
                                      end = soloLTR$end),
                     strand = soloLTR$strand,
                     class = soloLTR$class)
soloLTRGR <- soloLTRGR[seqnames(soloLTRGR) %in% chrName]
if(length(soloLTRGR) > 0) {
  soloLTR_bed <- data.frame(chr = as.character(seqnames(soloLTRGR)),
                            start = as.integer(start(soloLTRGR)-1),
                            end = as.integer(end(soloLTRGR)),
                            name = as.integer(1:length(soloLTRGR)),
                            score = as.character(soloLTRGR$class),
                            strand = as.character(strand(soloLTRGR)))
  write.table(soloLTR_bed,
              file = paste0(soloLTRDir, "CENsoloLTR_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
