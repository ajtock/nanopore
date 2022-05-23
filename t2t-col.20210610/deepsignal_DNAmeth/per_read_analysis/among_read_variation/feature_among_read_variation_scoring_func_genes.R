#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Score among-read variation/agreement (e.g., Fleiss' kappa) for each feature
# 2. Examine relationships between feature among-read agreement and other metrics

# Usage on hydrogen node7:
# csmit -m 200G -c 48 "/applications/R/R-4.0.0/bin/Rscript feature_among_read_variation_scoring_func_genes.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'gene' 'regions'"
 
# Divide each read into adjacent segments each consisting of a given number of consecutive cytosines,
# and calculate the methylation proportion for each segment of each read

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#featName <- "gene"
#featRegion <- "regions"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
context <- args[3]
NAmax <- as.numeric(args[4])
CPUpc <- as.numeric(args[5])
chrName <- unlist(strsplit(args[6], split = ","))
featName <- args[7]
featRegion <- args[8]

if(context == "CpG") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
} else if(context == "CHG") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
} else if(context == "CHH") {
  min_Cs <- 2
  max_Cs <- Inf
  min_reads <- 2
  max_reads <- Inf 
}

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
source("feature_among_read_variation_scoring_func_genes_function.R")
source("/projects/meiosis/ajt200/Rfunctions/TTSplus.R")
library(parallel)
library(GenomicRanges)
library(irr) #
library(dplyr)
library(tidyr)
library(cluster)
library(fpc) #
#library(data.table)
#library(segmentSeq)
library(ComplexHeatmap)
#library(RColorBrewer)
library(scales)
#library(circlize)
library(ggplot2)
library(cowplot)
#library(ggcorrplot)
library(viridis)
#library(ggthemes)
library(tidyquant) #
#library(grid)
library(doParallel)
library(doFuture)
registerDoFuture()
plan(multicore, workers = round(detectCores()*CPUpc))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())
options(future.globals.maxSize = 1000*1024^2)

outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
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

featextGR <- GRanges(seqnames = seqnames(featGR),
                     ranges = IRanges(start = start(featGR)-1000,
                                      end = end(featGR)+1000),
                     strand = strand(featGR),
                     name = featGR$name,
                     score = featGR$score)

# Mask out featGR within mitochondrial insertion on Chr2
fOverlaps_feat_mito_ins <- findOverlaps(query = featextGR,
                                        subject = mito_ins_GR,
                                        type = "any",
                                        select = "all",
                                        ignore.strand = T)
if(length(fOverlaps_feat_mito_ins) > 0) {
  featGR <- featGR[-unique(queryHits(fOverlaps_feat_mito_ins))]
}

# Define separate GRanges object for use with calculating exon and intron stats
genesGR <- featGR

# Get ranges corresponding to featRegion
if(featRegion == "bodies") {
  featGR <- featGR
} else if(featRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  featGR <- promoters(featGR, upstream = 1000, downstream = 0)
} else if(featRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/meiosis/ajt200/Rfunctions/TTSplus.R")
  featGR <- TTSplus(featGR, upstream = -1, downstream = 1000)
} else if(featRegion == "regions") {
  featGR <- GRanges(seqnames = seqnames(featGR),
                    ranges = IRanges(start = start(featGR)-1000,
                                     end = end(featGR)+1000),
                    strand = strand(featGR),
                    name = featGR$name,
                    score = featGR$score)
} else {
  stop("featRegion is none of bodies, promoters, terminators or regions")
}

# Load exons and introns
exons <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/genes/", refbase, "_representative_exons",
                           "_", paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
colnames(exons) <- c("chr", "start0based", "end", "name", "score", "strand")
exons <- exons[which(exons$name %in% genesGR$name),]
exonsGR <- GRanges(seqnames = exons$chr,
                   ranges = IRanges(start = exons$start0based+1,
                                    end = exons$end),
                   strand = exons$strand,
                   name = exons$name,
                   score = exons$score)

introns <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                             "/annotation/genes/", refbase, "_representative_introns",
                             "_", paste0(chrName, collapse = "_"), ".bed"),
                      header = F)
colnames(introns) <- c("chr", "start0based", "end", "name", "score", "strand")
introns <- introns[which(introns$name %in% genesGR$name),]
intronsGR <- GRanges(seqnames = introns$chr,
                     ranges = IRanges(start = introns$start0based+1,
                                      end = introns$end),
                     strand = introns$strand,
                     name = introns$name,
                     score = introns$score)

# Calculate number and total width of all exons or introns within each gene
genesIDs <- as.character(genesGR$name)
featGR_c <- foreach(x = iter(genesIDs),
                    .combine = "c",
                    .multicombine = T,
                    .maxcombine = length(genesIDs)+1e1,
                    .inorder = T,
                    .errorhandling = "pass") %dopar% {

  feat_genesID_x <- featGR[featGR$name == x]

  genes_genesID_x <- genesGR[genesGR$name == x]
  genes_genesID_x_width <- width(genes_genesID_x)

  exons_genesID_x <- exonsGR[exonsGR$name == x]
  exons_genesID_x_width <- sum(width(exons_genesID_x), na.rm = T)
  exons_genesID_x_width_prop <- exons_genesID_x_width / genes_genesID_x_width
  exons_genesID_x_count <- length(exons_genesID_x)
  exons_genesID_x_count_per_kb <- exons_genesID_x_count / (genes_genesID_x_width / 1e3)

  introns_genesID_x <- intronsGR[intronsGR$name == x]
  introns_genesID_x_width <- sum(width(introns_genesID_x), na.rm = T)
  introns_genesID_x_width_prop <- introns_genesID_x_width / genes_genesID_x_width
  introns_genesID_x_count <- length(introns_genesID_x)
  introns_genesID_x_count_per_kb <- introns_genesID_x_count / (genes_genesID_x_width / 1e3)

  GRanges(feat_genesID_x,
          name = x,
          score = NA,
          feature_width = genes_genesID_x_width,
          exons_width = exons_genesID_x_width,
          exons_width_prop = exons_genesID_x_width_prop,
          exons_count = exons_genesID_x_count,
          exons_count_per_kb = exons_genesID_x_count_per_kb,
          introns_width = introns_genesID_x_width,
          introns_width_prop = introns_genesID_x_width_prop,
          introns_count = introns_genesID_x_count,
          introns_count_per_kb = introns_genesID_x_count_per_kb)

} 
featGR <- featGR_c
stopifnot(identical(featGR$name, genesIDs))
rm(featGR_c); gc()

# Read in the raw output .tsv file from Deepsignal methylation model
tab_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                    sampleName, "_MappedOn_", refbase, "_", context, "_raw_", chrName[x], ".tsv"),
             header = F)
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  tab <- dplyr::bind_rows(tab_list)
} else {
  tab <- tab_list[[1]]
}
rm(tab_list); gc()

source("feature_among_read_variation_scoring_func_genes_function.R")

con_fk_df_all <- data.frame()
for(chrIndex in 1:length(chrName)) {

  print(chrName[chrIndex])

  chr_featGR <- featGR[seqnames(featGR) == chrName[chrIndex]]
  chr_tab <- tab[tab[,1] == chrName[chrIndex],]
  chr_tabGR <- GRanges(seqnames = chrName[chrIndex],
                       ranges = IRanges(start = chr_tab[,2]+1,
                                        width = 1),
                       strand = chr_tab[,3],
                       read = chr_tab[,5],
                       call = chr_tab[,9])
  chr_tabGR_fwd <- chr_tabGR[strand(chr_tabGR) == "+"]
  chr_tabGR_rev <- chr_tabGR[strand(chr_tabGR) == "-"]

  fOverlaps_fwd <- fOverlapsStrand(chr_tabGR_str = chr_tabGR_fwd, chr_featGR = chr_featGR)
  fOverlaps_rev <- fOverlapsStrand(chr_tabGR_str = chr_tabGR_rev, chr_featGR = chr_featGR)

  # Analyse each strand separately
  # fwd
#  makeDFx_list_fwd <- mclapply(1:length(chr_featGR), function(x) {
  makeDFx_list_fwd <- foreach(x = 1:length(chr_featGR),
                              .combine = "list",
                              .multicombine = T,
                              .maxcombine = length(chr_featGR)+1e1,
                              .inorder = T) %dopar% {
    makeDFx_strand(fOverlaps_str = fOverlaps_fwd,
                   chr_tabGR_str = chr_tabGR_fwd,
                   chr_featGR = chr_featGR,
                   featNum = x)
  }
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_fwd <- dplyr::bind_rows(makeDFx_list_fwd, .id = "column_label")
  
  chr_fk_df_fwd <- data.frame(chr_fk_df_fwd,
                              fk_adj_pval_str = p.adjust(chr_fk_df_fwd$fk_pval_str, method = "BH"))

  # rev  
#  makeDFx_list_rev <- mclapply(1:length(chr_featGR), function(x) {
  makeDFx_list_rev <- foreach(x = 1:length(chr_featGR),
                              .combine = "list",
                              .multicombine = T,
                              .maxcombine = length(chr_featGR)+1e1,
                              .inorder = T) %dopar% {
    makeDFx_strand(fOverlaps_str = fOverlaps_rev,
                   chr_tabGR_str = chr_tabGR_rev,
                   chr_featGR = chr_featGR,
                   featNum = x)
  }
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
   
  chr_fk_df_rev <- dplyr::bind_rows(makeDFx_list_rev, .id = "column_label")
  
  chr_fk_df_rev <- data.frame(chr_fk_df_rev,
                              fk_adj_pval_str = p.adjust(chr_fk_df_rev$fk_pval_str, method = "BH"))
  
  stopifnot(identical(chr_fk_df_fwd[,1:16], chr_fk_df_rev[,1:16]))
 
  # Take mean of fwd and rev equivalent columns 
  chr_fk_df_all_mean_list <- lapply(17:ncol(chr_fk_df_fwd), function(x) {
    sapply(1:nrow(chr_fk_df_fwd), function(y) {
      mean(c(chr_fk_df_fwd[y, x], chr_fk_df_rev[y, x]), na.rm = T)
    })
  })
  
  chr_fk_df_all <- data.frame(chr_fk_df_fwd[,1:16],
                              dplyr::bind_cols(chr_fk_df_all_mean_list))
  colnames(chr_fk_df_all) <- sub("_str", "_all", colnames(chr_fk_df_fwd))

  con_fk_df_all <- rbind(con_fk_df_all, chr_fk_df_all)

}

if(context == "CpG") {
  min_Cs <- 10
  max_Cs <- Inf
  min_reads <- 10
  max_reads <- Inf
} else if(context == "CHG") {
  min_Cs <- 10
  max_Cs <- Inf
  min_reads <- 10
  max_reads <- Inf
} else if(context == "CHH") {
  min_Cs <- 20
  max_Cs <- Inf
  min_reads <- 10
  max_reads <- Inf
}


con_fk_df_all_filt <- con_fk_df_all %>%
  dplyr::filter(fk_reads_all >= min_reads) %>%
  dplyr::filter(fk_Cs_all >= min_Cs)


write.table(con_fk_df_all,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_unfilt_df_fk_kappa_all_mean_mC_all_complete_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(con_fk_df_all_filt,
            paste0(outDir,
                   featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                   "_", context,
                   "_NAmax", NAmax,
                   "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                   paste0(chrName, collapse = "_"), ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
