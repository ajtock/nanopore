#!/applications/R/R-4.0.0/bin/Rscript

# Create table of representative genes
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_mRNA_coords.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
CENstart <- c(14841110,3823792,13597188,4203902,11784131)[which(fai$V1 %in% chrName)]
CENend <- c(17559778,6045243,15733925,6977949,14551809)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
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
write.table(mRNA_rep[,1:9],
            file = paste0("t2t-col.20210610_representative_mRNA_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_bed <- data.frame(chr = as.character(mRNA_rep[,1]),
                           start = as.integer(mRNA_rep[,4]-1),
                           end = as.integer(mRNA_rep[,5]),
                           name = as.character(mRNA_rep[,9]),
                           score = as.numeric(mRNA_rep[,6]),
                           strand = as.character(mRNA_rep[,7]))
write.table(mRNA_rep_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

mRNA_repGR <- GRanges(seqnames = mRNA_rep$seqid,
                      ranges = IRanges(start = mRNA_rep$start,
                                       end = mRNA_rep$end),
                      strand = mRNA_rep$strand)

# Define function to select randomly positioned loci of the same
# width distribution as mRNA_repGR
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
# ranLocGR contains the same number of loci per chromosome as mRNA_repGR
chrs <- chrs[chrs %in% chrName]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  mRNA_repChrGR <- mRNA_repGR[seqnames(mRNA_repGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(mRNA_repChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(mRNA_repChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(mRNA_repChrGR)),
                         strand = strand(mRNA_repChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0("t2t-col.20210610_representative_mRNA_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
