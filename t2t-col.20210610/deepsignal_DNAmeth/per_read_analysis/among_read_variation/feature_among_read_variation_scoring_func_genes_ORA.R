#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Gene set enrichment analysis of genes grouped by both among-read agreement and methylation proportion

# Usage on hydrogen node7:
# csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_genes_ORA.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene'"
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]

#if(context == "CpG") {
#  min_Cs <- 2
#  max_Cs <- Inf
#  min_reads <- 2
#  max_reads <- Inf 
#} else if(context == "CHG") {
#  min_Cs <- 2
#  max_Cs <- Inf
#  min_reads <- 2
#  max_reads <- Inf 
#} else if(context == "CHH") {
#  min_Cs <- 2
#  max_Cs <- Inf
#  min_reads <- 2
#  max_reads <- Inf 
#}

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.At.tair.db", character.only = TRUE)
library("org.At.tair.db", character.only = TRUE)
#BiocManager::install("pathview") # installation of package ‘Rgraphviz’ had non-zero exit status; installation of package ‘pathview’ had non-zero exit status
#library(pathview)
#BiocManager::install("enrichplot")
library(enrichplot)
#BiocManager::install("DOSE")
library(DOSE)
library(ggplot2)

#library(parallel)
#library(GenomicRanges)
#library(irr)
#library(dplyr)
#library(tidyr)
#library(cluster)
#library(fpc)
##library(data.table)
##library(segmentSeq)
#library(ComplexHeatmap)
##library(RColorBrewer)
#library(scales)
##library(circlize)
# 
#library(ggplot2)
#library(cowplot)
##library(ggcorrplot)
#library(viridis)
#library(ggthemes)
#library(tidyquant)
##library(grid)

outDir <- paste0(featName, "/")
plotDir <- paste0(outDir, "plots/GO_BP_ORA/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Read in feature annotation
if(featName == "CEN180") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/CEN180/CEN180_in_", refbase,
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand",
                      "HORlengthsSum", "HORcount", "percentageIdentity")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$HORlengthsSum)
} else if(featName == "gene") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/genes/", refbase, "_representative_mRNA",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$score)
} else if(featName == "GYPSY") {
  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
  featGR <- GRanges(seqnames = feat$chr,
                    ranges = IRanges(start = feat$start0based+1,
                                     end = feat$end),
                    strand = feat$strand,
                    name = feat$name,
                    score = feat$score)
} else {
  stop(print("featName not one of CEN180, gene or GYPSY"))
}

# Load coordinates for mitochondrial insertion on Chr2, in BED format
mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
                       header = F)
colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
mito_ins_GR <- GRanges(seqnames = "Chr2",
                       ranges = IRanges(start = min(mito_ins$start)+1,
                                        end = max(mito_ins$end)),
                       strand = "*")

# Mask out featGR within mitochondrial insertion on Chr2
fOverlaps_feat_mito_ins <- findOverlaps(query = featGR,
                                        subject = mito_ins_GR,
                                        type = "any",
                                        select = "all",
                                        ignore.strand = T)
if(length(fOverlaps_feat_mito_ins) > 0) {
  featGR <- featGR[-unique(queryHits(fOverlaps_feat_mito_ins))]
}


if(context == "CpG") {
  fk_kappa_all_high <- 0.55
  fk_kappa_all_mid  <- 0.35
  fk_kappa_all_low  <- 0.04
  mean_stocha_all_high <- 0.28
  mean_stocha_all_mid  <- 0.17
  mean_stocha_all_low  <- 0.08
  mean_min_acf_all_high <- -0.05
  mean_min_acf_all_mid  <- -0.10
  mean_min_acf_all_low  <- -0.15
  mean_mC_all_high  <- 0.75
  mean_mC_all_mid   <- 0.25
  mean_mC_all_low   <- 0.10
} else if(context == "CHG") {
  fk_kappa_all_high <- 0.05623413
  fk_kappa_all_mid  <- 0.01778279 
  fk_kappa_all_low  <- 0.01000000 
  mean_stocha_all_high <- 0.28
  mean_stocha_all_mid  <- 0.17
  mean_stocha_all_low  <- 0.08
  mean_min_acf_all_high <- -0.05
  mean_min_acf_all_mid  <- -0.10
  mean_min_acf_all_low  <- -0.15
  mean_mC_all_high  <- 0.05623413
  mean_mC_all_mid   <- 0.02511886
  mean_mC_all_low   <- 0.01778279
} else if(context == "CHH") {
  fk_kappa_all_high <- 0.05623413
  fk_kappa_all_mid  <- 0.01778279 
  fk_kappa_all_low  <- 0.01000000 
  mean_stocha_all_high <- 0.28
  mean_stocha_all_mid  <- 0.17
  mean_stocha_all_low  <- 0.08
  mean_min_acf_all_high <- -0.05
  mean_min_acf_all_mid  <- -0.10
  mean_min_acf_all_low  <- -0.15
  mean_mC_all_high  <- 0.02500000
  mean_mC_all_mid   <- 0.01000000
  mean_mC_all_low   <- 0.00500000
}


# Load feature groups (based on trend plots) to enable enrichment analysis

#write.table(con_fk_df_all_filt_kappa_low_mC_low_group1,

filt_kappa_mC_groups <- lapply(seq_along(1:8), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_mC_all_group", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T)
  tmp$name <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
  tmp$name <- sub(pattern = "_\\d+", replacement = "", x = tmp$name)
  tmp[ with(tmp, order(fk_kappa_all, decreasing = T)), ]
})


filt_kappa_mC_groups_featID <- lapply(seq_along(filt_kappa_mC_groups), function(x) {
  tmp <- filt_kappa_mC_groups[[x]]$fk_kappa_all
  names(tmp) <- filt_kappa_mC_groups[[x]]$name
  na.omit(tmp)
})

filt_kappa_mC_groups_featID_universe <- unlist(filt_kappa_mC_groups_featID)

keytypes(org.At.tair.db)
head(select(org.At.tair.db, keys = keys(org.At.tair.db), columns = c("TAIR")))

filt_kappa_mC_groups_enrichGO <- mclapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
  enrichGO(gene = names(filt_kappa_mC_groups_featID[[x]]),
           universe = names(filt_kappa_mC_groups_featID_universe),
           OrgDb = org.At.tair.db,
           keyType = "TAIR",
           readable = T,
           ont = "BP",
#           minGSSize = 10,
#           maxGSSize = 500,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.10,
           pAdjustMethod = "BH")
}, mc.cores = length(filt_kappa_mC_groups_featID), mc.preschedule = F)


for(x in 1:length(filt_kappa_mC_groups_enrichGO)) {
  dp_enrichGO <- dotplot(filt_kappa_mC_groups_enrichGO[[x]],
                         showCategory = 50,
                         title = paste0("Fleiss' kappa and m", context, " Gene Group ", x),
                         font.size = 12)
  ggsave(paste0(plotDir,
                featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                "_NAmax", NAmax,
                "_filt_df_fk_kappa_all_mean_mC_all_group", x , "_",
                paste0(chrName, collapse = "_"), "_enrichGO_dotplot.pdf"),
         plot = dp_enrichGO, height = 20, width = 50, limitsize = F)
}


write.table(con_fk_df_all_filt_kappa_mid_mC_low_group2,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mid_mC_mid_group3,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_high_mC_mid_group4,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_high_mC_high_group5,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_vhigh_mC_high_group6,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_low_mC_vhigh_group7,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_kappa_mid_mC_vhigh_group8,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


# Filter by mean_stocha_all and mean_mC_all
con_fk_df_all_filt_stocha_low_mC_low_group1 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_low)

con_fk_df_all_filt_stocha_mid_mC_low_group2 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_low)

con_fk_df_all_filt_stocha_mid_mC_mid_group3 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

con_fk_df_all_filt_stocha_high_mC_mid_group4 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_low) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_mid)

con_fk_df_all_filt_stocha_high_mC_high_group5 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <=  mean_stocha_all_high) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_high)

con_fk_df_all_filt_stocha_vhigh_mC_high_group6 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all >  mean_stocha_all_high) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_mid) %>%
  dplyr::filter(mean_mC_all     <= mean_mC_all_high)

con_fk_df_all_filt_stocha_mid_mC_vhigh_group7 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all <= mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_high)

con_fk_df_all_filt_stocha_high_mC_vhigh_group8 <- con_fk_df_all_filt %>%
  dplyr::filter(mean_stocha_all > mean_stocha_all_mid) %>%
  dplyr::filter(mean_mC_all     >  mean_mC_all_high)

write.table(con_fk_df_all_filt_stocha_low_mC_low_group1,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group1_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_low_group2,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group2_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_mid_group3,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group3_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_mid_group4,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group4_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_high_group5,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group5_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_vhigh_mC_high_group6,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group6_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_mid_mC_vhigh_group7,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group7_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt_stocha_high_mC_vhigh_group8,
            paste0(outDir,
                   featName, "_", sampleName, "_MappedOn_", refbase, "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_mean_stocha_all_mean_mC_all_group8_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)


