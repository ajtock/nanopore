#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. GO term over-representation analysis of genes grouped by both among-read agreement and mean methylation proportion
# 1. GO term over-representation analysis of genes grouped by both mean within-read stochasticity and mean methylation proportion

# Usage:
# /applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_genes_ORA.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'bodies' 'BP'
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "bodies"
#ontology = "BP"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]
featRegion <- args[8]
ontology <- args[9]

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

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir_kappa_mC <- paste0(outDir, "plots/ORA_GO_", ontology, "_", context, "_kappa_mC/")
plotDir_stocha_mC <- paste0(outDir, "plots/ORA_GO_", ontology, "_", context, "_stocha_mC/")
plotDir_kappa_stocha <- paste0(outDir, "plots/ORA_GO_", ontology, "_", context, "_kappa_stocha/")
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
## Read in feature annotation
#if(featName == "CEN180") {
#  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
#                            "/annotation/CEN180/CEN180_in_", refbase,
#                            "_", paste0(chrName, collapse = "_"), ".bed"),
#                     header = F)
#  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand",
#                      "HORlengthsSum", "HORcount", "percentageIdentity")
#  featGR <- GRanges(seqnames = feat$chr,
#                    ranges = IRanges(start = feat$start0based+1,
#                                     end = feat$end),
#                    strand = feat$strand,
#                    name = feat$name,
#                    score = feat$HORlengthsSum)
#} else if(featName == "gene") {
#  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
#                            "/annotation/genes/", refbase, "_representative_mRNA",
#                            "_", paste0(chrName, collapse = "_"), ".bed"),
#                     header = F)
#  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
#  featGR <- GRanges(seqnames = feat$chr,
#                    ranges = IRanges(start = feat$start0based+1,
#                                     end = feat$end),
#                    strand = feat$strand,
#                    name = feat$name,
#                    score = feat$score)
#} else if(featName == "GYPSY") {
#  feat <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
#                            "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
#                            "_", paste0(chrName, collapse = "_"), ".bed"),
#                     header = F)
#  colnames(feat) <- c("chr", "start0based", "end", "name", "score", "strand")
#  featGR <- GRanges(seqnames = feat$chr,
#                    ranges = IRanges(start = feat$start0based+1,
#                                     end = feat$end),
#                    strand = feat$strand,
#                    name = feat$name,
#                    score = feat$score)
#} else {
#  stop(print("featName not one of CEN180, gene or GYPSY"))
#}
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
