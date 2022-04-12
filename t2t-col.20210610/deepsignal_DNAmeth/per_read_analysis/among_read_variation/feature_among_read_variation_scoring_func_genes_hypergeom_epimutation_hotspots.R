#!/usr/bin/env Rscript

# Analysis:
# 1. Epimutation hotspot over- and under-representation analysis of genes grouped by both among-read agreement and mean methylation proportion
# 1. Epimutation hotspot over- and under-representation analysis of genes grouped by both mean within-read stochasticity and mean methylation proportion

# Usage:
# conda activate R-4.0.0
# ./feature_among_read_variation_scoring_func_genes_hypergeom_epimutation_hotspots.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions' 10000 1000
# conda deactivate
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "regions"
#genomeBinSize <- 10000
#genomeStepSize <- 1000

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
chrName <- unlist(strsplit(args[5], split = ","))
featName <- args[6]
featRegion <- args[7]
genomeBinSize <- as.numeric(args[8])
genomeStepSize <- as.numeric(args[9])

options(stringsAsFactors = F)
library(parallel)
library(dplyr)
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

library(methods)
library(plotrix)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

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


if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
  genomeBinNamePlot <- paste0(genomeBinSize, "-bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e3, "-kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e6, "-Mb")
}

if(floor(log10(genomeStepSize)) + 1 < 4) {
  genomeStepName <- paste0(genomeStepSize, "bp")
  genomeStepNamePlot <- paste0(genomeStepSize, "-bp")
} else if(floor(log10(genomeStepSize)) + 1 >= 4 &
          floor(log10(genomeStepSize)) + 1 <= 6) {
  genomeStepName <- paste0(genomeStepSize/1e3, "kb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e3, "-kb")
} else if(floor(log10(genomeStepSize)) + 1 >= 7) {
  genomeStepName <- paste0(genomeStepSize/1e6, "Mb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e6, "-Mb")
}


outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir_kappa_mC <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_stocha_mC/")
plotDir_kappa_stocha <- paste0(outDir, "plots/hypergeom_epimutation_hotspots_", context, "_kappa_stocha/")
system(paste0("[ -d ", plotDir_kappa_mC, " ] || mkdir -p ", plotDir_kappa_mC))
system(paste0("[ -d ", plotDir_stocha_mC, " ] || mkdir -p ", plotDir_stocha_mC))
system(paste0("[ -d ", plotDir_kappa_stocha, " ] || mkdir -p ", plotDir_kappa_stocha))

## Genomic definitions
#fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
#if(!grepl("Chr", fai[,1][1])) {
#  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
#} else {
#  chrs <- fai[,1][which(fai[,1] %in% chrName)]
#}
#chrLens <- fai[,2][which(fai[,1] %in% chrName)]
#
## Load coordinates for mitochondrial insertion on Chr2, in BED format
#mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
#                       header = F)
#colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
#mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
#mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
#mito_ins_GR <- GRanges(seqnames = "Chr2",
#                       ranges = IRanges(start = min(mito_ins$start)+1,
#                                        end = max(mito_ins$end)),
#                       strand = "*")
#
## Mask out featGR within mitochondrial insertion on Chr2
#fOverlaps_feat_mito_ins <- findOverlaps(query = featGR,
#                                        subject = mito_ins_GR,
#                                        type = "any",
#                                        select = "all",
#                                        ignore.strand = T)
#if(length(fOverlaps_feat_mito_ins) > 0) {
#  featGR <- featGR[-unique(queryHits(fOverlaps_feat_mito_ins))]
#}

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]



# Load coordinates for mitochondrial insertion on Chr2, in BED format
mito_ins <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/", refbase , ".mitochondrial_insertion.bed"),
                       header = F)
colnames(mito_ins) <- c("chr", "start", "end", "name", "score", "strand")
mito_ins <- mito_ins[ mito_ins$chr %in% "Chr2",]
mito_ins <- mito_ins[ with(mito_ins, order(chr, start, end)) , ]
mito_ins_chr <- unique(mito_ins$chr)
mito_ins_start <- min(mito_ins$start)+1
mito_ins_end <- max(mito_ins$end)
#mito_ins_GR <- GRanges(seqnames = "Chr2",
#                       ranges = IRanges(start = min(mito_ins$start)+1,
#                                        end = max(mito_ins$end)),
#                       strand = "*")

# Read in methylation divergence (mD) files
mDdir <- "/home/ajt200/analysis/BSseq_leaf_Hazarika_Johannes_2022_NatPlants/snakemake_BSseq_t2t-col.20210610/MA1_2/coverage/report/alphabeta/"
MA1_2_mD_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0(mDdir,
                    "mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                    "_MA1_2_MappedOn_", refbase, "_", chrName[x], "_", context, ".tsv"),
             header = T)
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  MA1_2_mD <- dplyr::bind_rows(MA1_2_mD_list)
} else {
  MA1_2_mD <- MA1_2_mD_list[[1]]
}
rm(MA1_2_mD_list); gc()

tab_mD <- MA1_2_mD
rm(MA1_2_mD); gc()

tab_mD <- tab_mD[
  with(tab_mD, order(chr, start, end)),
]

# Remove genomic bins in tab_mD overlapping the mitochondrial insertion on Chr2,
# as we cannot be sure that these BS-seq reads come from the nuclear genome
# Note that long reads wholly contained within the boundaries of the Chr2 mito insertion
# have already been removed for the among-read and within-read mC variation analysis,
# by among_read_variation_scoring.R
tab_mD_mito <- tab_mD %>%
  dplyr::filter(chr == mito_ins_chr) %>%
  dplyr::filter( (start >= mito_ins_start &
                  start <= mito_ins_end) |
                 (end >= mito_ins_start &
                  end <= mito_ins_end) )

stopifnot(unique(tab_mD_mito$chr) == mito_ins_chr)
stopifnot(tab_mD_mito$end[1] >= mito_ins_start)
stopifnot(tab_mD_mito$start[nrow(tab_mD_mito)] <= mito_ins_end)
print(head(tab_mD_mito[,1:6]))
print(tail(tab_mD_mito[,1:6]))
print(paste0("Mitochondrial insetion on ", mito_ins_chr, ":", mito_ins_start, "-", mito_ins_end))

tab_mD <- tab_mD %>%
  dplyr::anti_join(tab_mD_mito)

tab_mD <- tab_mD[!is.na(tab_mD$MA1_2_mean.D),]

tab_mD$rank <- rank(tab_mD$MA1_2_mean.D)/nrow(tab_mD)

mD_hs <- tab_mD[tab_mD$rank >= 0.9,]



# Load feature groups (defined based on Fleiss' kappa vs mean mC trend plots) to enable enrichment analysis

filt_kappa_mC_groups <- lapply(seq_along(1:8), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
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
ids <- select(org.At.tair.db, keys = keys(org.At.tair.db), columns = c("TAIR"))[,1]

filt_kappa_mC_groups_enrichGO <- lapply(seq_along(filt_kappa_mC_groups_featID), function(x) {
  print(x)
  enrichGO(gene = names(filt_kappa_mC_groups_featID[[x]]),
           universe = names(filt_kappa_mC_groups_featID_universe),
           OrgDb = org.At.tair.db,
           keyType = "TAIR",
           readable = T,
           ont = ontology,
#           minGSSize = 10,
#           maxGSSize = 500,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.10,
           pAdjustMethod = "BH")
})

for(x in 1:length(filt_kappa_mC_groups_enrichGO)) {
  if( !is.null(filt_kappa_mC_groups_enrichGO[[x]]) ) {
    if( sum(filt_kappa_mC_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 0 ) {
      print(x)
      dp_enrichGO <- dotplot(filt_kappa_mC_groups_enrichGO[[x]],
                             showCategory = 50,
                             title = paste0("Fleiss' kappa and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                             font.size = 12)
      ggsave(paste0(plotDir_kappa_mC,
                    featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                    "_", context,
                    "_NAmax", NAmax,
                    "_filt_df_fk_kappa_all_mean_mC_all_group", x , "_",
                    paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_dotplot.pdf"),
             plot = dp_enrichGO,
             height = 10, width = 12,
             limitsize = F)

      if(sum(filt_kappa_mC_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 1) {
        emp_enrichGO <- emapplot(filt_kappa_mC_groups_enrichGO[[x]],
                                 showCategory = 50,
                                 title = paste0("Fleiss' kappa and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                                 font.size = 12)
        ggsave(paste0(plotDir_kappa_mC,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_fk_kappa_all_mean_mC_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_emapplot.pdf"),
               plot = emp_enrichGO,
               height = 10, width = 12,
               limitsize = F)
  
        gp_enrichGO <- goplot(filt_kappa_mC_groups_enrichGO[[x]],
                              showCategory = 50,
                              title = paste0("Fleiss' kappa and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                              font.size = 12)
        ggsave(paste0(plotDir_kappa_mC,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_fk_kappa_all_mean_mC_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_goplot.pdf"),
               plot = gp_enrichGO,
               height = 10, width = 12,
               limitsize = F)
      }
    }
  }
}


# Load feature groups (defined based on mean stocha vs mean mC trend plots) to enable enrichment analysis

filt_stocha_mC_groups <- lapply(seq_along(1:8), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_mean_stocha_all_mean_mC_all_group", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T)
  tmp$name <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
  tmp$name <- sub(pattern = "_\\d+", replacement = "", x = tmp$name)
  tmp[ with(tmp, order(mean_stocha_all, decreasing = T)), ]
})

filt_stocha_mC_groups_featID <- lapply(seq_along(filt_stocha_mC_groups), function(x) {
  tmp <- filt_stocha_mC_groups[[x]]$mean_stocha_all
  names(tmp) <- filt_stocha_mC_groups[[x]]$name
  na.omit(tmp)
})

filt_stocha_mC_groups_featID_universe <- unlist(filt_stocha_mC_groups_featID)

filt_stocha_mC_groups_enrichGO <- lapply(seq_along(filt_stocha_mC_groups_featID), function(x) {
  print(x)
  enrichGO(gene = names(filt_stocha_mC_groups_featID[[x]]),
           universe = names(filt_stocha_mC_groups_featID_universe),
           OrgDb = org.At.tair.db,
           keyType = "TAIR",
           readable = T,
           ont = ontology,
#           minGSSize = 10,
#           maxGSSize = 500,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.10,
           pAdjustMethod = "BH")
})

for(x in 1:length(filt_stocha_mC_groups_enrichGO)) {
  if( !is.null(filt_stocha_mC_groups_enrichGO[[x]]) ) {
    if( sum(filt_stocha_mC_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 0 ) {
      print(x)
      dp_enrichGO <- dotplot(filt_stocha_mC_groups_enrichGO[[x]],
                             showCategory = 50,
                             title = paste0("Stochasticity and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                             font.size = 12)
      ggsave(paste0(plotDir_stocha_mC,
                    featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                    "_", context,
                    "_NAmax", NAmax,
                    "_filt_df_mean_stocha_all_mean_mC_all_group", x , "_",
                    paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_dotplot.pdf"),
             plot = dp_enrichGO,
             height = 10, width = 12,
             limitsize = F)

      if(sum(filt_stocha_mC_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 1) {
        emp_enrichGO <- emapplot(filt_stocha_mC_groups_enrichGO[[x]],
                                 showCategory = 50,
                                 title = paste0("Stochasticity and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                                 font.size = 12)
        ggsave(paste0(plotDir_stocha_mC,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_mean_stocha_all_mean_mC_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_emapplot.pdf"),
               plot = emp_enrichGO,
               height = 10, width = 12,
               limitsize = F)

        gp_enrichGO <- goplot(filt_stocha_mC_groups_enrichGO[[x]],
                              showCategory = 50,
                              title = paste0("Stochasticity and mean m", context, " in ", featName, " ", featRegion, " Group ", x),
                              font.size = 12)
        ggsave(paste0(plotDir_stocha_mC,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_mean_stocha_all_mean_mC_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_goplot.pdf"),
               plot = gp_enrichGO,
               height = 10, width = 12,
               limitsize = F)
      }
    }
  }
}


# Load feature groups (defined based on Fleiss' kappa vs mean stocha trend plots) to enable enrichment analysis

filt_kappa_stocha_groups <- lapply(seq_along(1:8), function(x) {
  tmp <- read.table(paste0(outDir,
                           featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                           "_", context,
                           "_NAmax", NAmax,
                           "_filt_df_fk_kappa_all_mean_stocha_all_group", x , "_",
                           paste0(chrName, collapse = "_"), ".tsv"),
                    header = T)
  tmp$name <- sub(pattern = "\\.\\d+", replacement = "", x = tmp$name)
  tmp$name <- sub(pattern = "_\\d+", replacement = "", x = tmp$name)
  tmp[ with(tmp, order(fk_kappa_all, decreasing = T)), ]
})

filt_kappa_stocha_groups_featID <- lapply(seq_along(filt_kappa_stocha_groups), function(x) {
  tmp <- filt_kappa_stocha_groups[[x]]$fk_kappa_all
  names(tmp) <- filt_kappa_stocha_groups[[x]]$name
  na.omit(tmp)
})

filt_kappa_stocha_groups_featID_universe <- unlist(filt_kappa_stocha_groups_featID)

keytypes(org.At.tair.db)
ids <- select(org.At.tair.db, keys = keys(org.At.tair.db), columns = c("TAIR"))[,1]

filt_kappa_stocha_groups_enrichGO <- lapply(seq_along(filt_kappa_stocha_groups_featID), function(x) {
  print(x)
  enrichGO(gene = names(filt_kappa_stocha_groups_featID[[x]]),
           universe = names(filt_kappa_stocha_groups_featID_universe),
           OrgDb = org.At.tair.db,
           keyType = "TAIR",
           readable = T,
           ont = ontology,
#           minGSSize = 10,
#           maxGSSize = 500,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.10,
           pAdjustMethod = "BH")
})

for(x in 1:length(filt_kappa_stocha_groups_enrichGO)) {
  if( !is.null(filt_kappa_stocha_groups_enrichGO[[x]]) ) {
    if( sum(filt_kappa_stocha_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 0 ) {
      print(x)
      dp_enrichGO <- dotplot(filt_kappa_stocha_groups_enrichGO[[x]],
                             showCategory = 50,
                             title = paste0("Fleiss' kappa and mean stochasticity (m", context, ") in ", featName, " ", featRegion, " Group ", x),
                             font.size = 12)
      ggsave(paste0(plotDir_kappa_stocha,
                    featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                    "_", context,
                    "_NAmax", NAmax,
                    "_filt_df_fk_kappa_all_mean_stocha_all_group", x , "_",
                    paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_dotplot.pdf"),
             plot = dp_enrichGO,
             height = 10, width = 12,
             limitsize = F)

      if(sum(filt_kappa_stocha_groups_enrichGO[[x]]@result$p.adjust <= 0.05) > 1) {
        emp_enrichGO <- emapplot(filt_kappa_stocha_groups_enrichGO[[x]],
                                 showCategory = 50,
                                 title = paste0("Fleiss' kappa and mean stochasticity (m", context, ") in ", featName, " ", featRegion, " Group ", x),
                                 font.size = 12)
        ggsave(paste0(plotDir_kappa_stocha,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_fk_kappa_all_mean_stocha_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_emapplot.pdf"),
               plot = emp_enrichGO,
               height = 10, width = 12,
               limitsize = F)
  
        gp_enrichGO <- goplot(filt_kappa_stocha_groups_enrichGO[[x]],
                              showCategory = 50,
                              title = paste0("Fleiss' kappa and mean stochasticity (m", context, ") in ", featName, " ", featRegion, " Group ", x),
                              font.size = 12)
        ggsave(paste0(plotDir_kappa_stocha,
                      featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                      "_", context,
                      "_NAmax", NAmax,
                      "_filt_df_fk_kappa_all_mean_stocha_all_group", x , "_",
                      paste0(chrName, collapse = "_"), "_enrichGO_", ontology, "_goplot.pdf"),
               plot = gp_enrichGO,
               height = 10, width = 12,
               limitsize = F)
      }
    }
  }
}
