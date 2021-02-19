#!/applications/R/R-4.0.0/bin/Rscript

# Convert gene Liftoff (TAIR10 to T2T) from GFF3 to GTF

# Usage:
# ./convert_gene_Liftoff_GFF3_to_GTF.R

options(stringsAsFactors = F)
library(rtracklayer)

outDir <- "genes/"
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

# Load table of gene coordinates in T2T_Col (GFF3 derived from Liftoff tool)
gff3 <- import.gff3("T2T_Col.genes.gff3")
export.gff2(object = gff3,
            con = "genes/T2T_Col.genes.gtf")
