#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Genome-wide autocorrelation of context-specific cytosine methylation status
# at increasing physical distances (e.g., 1 to 10,000 nucleotides)

# Usage on hydrogen node7:
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript autocorrelation_genes.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CpG 1e3 1.00 Chr1,Chr2,Chr3,Chr4,Chr5 1000"

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#nperm <- 1e3
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#maxDist <- 1000
#min_pval <- 1 - ( (nperm - 1) / nperm )

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
genomeBinSize <- as.integer(args[3])
genomeStepSize <- as.integer(args[4])
context <- args[5]
nperm <- as.numeric(args[6])
CPUpc <- as.numeric(args[7])
chrName <- unlist(strsplit(args[8], split = ","))
maxDist <- as.numeric(args[9])
min_pval <- 1 - ( (nperm - 1) / nperm )

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
#library(irr)
library(dplyr)
#library(tidyr)
#library(fpc)
##library(data.table)
##library(segmentSeq)
#library(ComplexHeatmap)
##library(RColorBrewer)
##library(viridis)
##library(scales)
##library(circlize)
 
library(ggplot2)
library(cowplot)
#library(ggcorrplot)
library(viridis)
library(ggthemes)
library(tidyquant)
#library(grid)

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

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

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
mito_ins_GR <- GRanges(seqnames = "Chr2",
                       ranges = IRanges(start = min(mito_ins$start)+1,
                                        end = max(mito_ins$end)),
                       strand = "*")

# Read in CEN180 annotation
CEN180 <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                            "/annotation/CEN180/CEN180_in_", refbase,
                            "_", paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
colnames(CEN180) <- c("chr", "start0based", "end", "name", "score", "strand", "HORlengthsSum", "HORcount")
CEN180GR <- GRanges(seqnames = CEN180$chr,
                    ranges = IRanges(start = CEN180$start0based+1,
                                     end = CEN180$end),
                    strand = CEN180$strand)

# Read in gene annotation
genes <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/genes/", refbase, "_representative_mRNA",
                           "_", paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
colnames(genes) <- c("chr", "start0based", "end", "name", "score", "strand")
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start0based+1,
                                    end = genes$end),
                   strand = genes$strand)

# Read in gypsy annotation
gypsy <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                           "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
                           "_", paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
colnames(gypsy) <- c("chr", "start0based", "end", "name", "score", "strand")
gypsyGR <- GRanges(seqnames = gypsy$chr,
                   ranges = IRanges(start = gypsy$start0based+1,
                                    end = gypsy$end),
                   strand = gypsy$strand)


# Read in the "methylation frequency" output .tsv file from Deepsignal methylation model
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                         sampleName, "_MappedOn_", refbase, "_", context, ".tsv"),
                  header = F)
tab <- tab[ order(tab[,1], tab[,2], tab[,3], decreasing = F), ]
tab <- tab[tab[,1] %in% chrName,]
tabGR <- GRanges(seqnames = tab[,1],
                 ranges = IRanges(start = tab[,2],
                                  end = tab[,2]),
                 strand = tab[,3],
                 prop = tab[,10])

# Mask out methylation info within mitochondrial insertion on Chr2
fOverlaps_tab_mito_ins <- findOverlaps(query = tabGR,
                                       subject = mito_ins_GR,
                                       type = "any",
                                       select = "all",
                                       ignore.strand = T)
if(length(fOverlaps_tab_mito_ins) > 0) {
  tabGR <- tabGR[-unique(queryHits(fOverlaps_tab_mito_ins))]
}

# Find context-specific cytosine sites that overlap genes
fOverlaps_tab_genes <- findOverlaps(query = tabGR,
                                    subject = genesGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = T)
tabGR_genes <- tabGR[unique(queryHits(fOverlaps_tab_genes))]

# Analyse each strand separately
tabGR_genes_fwd <- tabGR_genes[strand(tabGR_genes) == "+"]
tabGR_genes_fwd <- sortSeqlevels(tabGR_genes_fwd)
tabGR_genes_fwd <- sort(tabGR_genes_fwd, by = ~ seqnames + start + end)

# Randomly shuffle methylation proportion values over the cytosine coordinates
tabGR_genes_fwd_random <- mclapply(1:nperm, function(y) {
  tabGR_genes_fwd_chr_list_y <- lapply(seq_along(chrName), function(x) {
    GRanges(seqnames = seqnames(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]]),
            ranges = ranges(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]]),
            strand = strand(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]]),
            prop = sample(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]]$prop))
  })
  do.call(c, as(tabGR_genes_fwd_chr_list_y, "GRangesList")) 
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

#  sort( unique ( start( tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == x][2:(length(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == x]))] ) -
#                 start( tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == x][1:(length(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == x])-1)] ) ) )

# Get context-specific inter-cytosine bp distances that occur more than once
tabGR_genes_fwd_dists_list <- mclapply(seq_along(chrName), function(x) {
  as.integer( names (
    which( table(
      sort( diff( start( tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]] ) ) )
   ) > 2 )
  ) )
}, mc.cores = length(chrName), mc.preschedule = F)


# Analyse each strand separately
tabGR_genes_rev <- tabGR_genes[strand(tabGR_genes) == "-"]
tabGR_genes_rev <- sortSeqlevels(tabGR_genes_rev)
tabGR_genes_rev <- sort(tabGR_genes_rev, by = ~ seqnames + start + end)

# Randomly shuffle methylation proportion values over the cytosine coordinates
tabGR_genes_rev_random <- mclapply(1:nperm, function(y) {
  tabGR_genes_rev_chr_list_y <- lapply(seq_along(chrName), function(x) {
    GRanges(seqnames = seqnames(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
            ranges = ranges(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
            strand = strand(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
            prop = sample(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]$prop))
  })
  do.call(c, as(tabGR_genes_rev_chr_list_y, "GRangesList")) 
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

#  sort( unique ( start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == x][2:(length(tabGR_genes_rev[seqnames(tabGR_genes_rev) == x]))] ) -
#                 start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == x][1:(length(tabGR_genes_rev[seqnames(tabGR_genes_rev) == x])-1)] ) ) )

# Get context-specific inter-cytosine bp distances that occur more than once
tabGR_genes_rev_dists_list <- mclapply(seq_along(chrName), function(x) {
  as.integer( names (
    which( table(
      sort( diff( start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]] ) ) )
   ) > 2 )
  ) )
}, mc.cores = length(chrName), mc.preschedule = F)

# Get the distance that is the lowest common denominator among the chromsomes
tabGR_genes_all_dists_min_max <- min(sapply(c(tabGR_genes_fwd_dists_list, tabGR_genes_rev_dists_list), function(x) {
  max(x, na.rm = T)
}))

if(tabGR_genes_all_dists_min_max > maxDist) {
  tabGR_genes_all_dists_min_max <- maxDist
}

tabGR_genes_all_dists_list <- lapply(seq_along(chrName), function(x) {
  unique(tabGR_genes_fwd_dists_list[[x]], tabGR_genes_rev_dists_list[[x]])
})

tabGR_genes_all_dists_list <- lapply(tabGR_genes_all_dists_list, function(x) {
  x[which(x <= tabGR_genes_all_dists_min_max)]
})

acfDistance <- function(DSfreqGR, bpDistance) {
  cor(x = DSfreqGR[ ( which( diff( start(DSfreqGR)) == bpDistance) )]$prop,
      y = DSfreqGR[ ( which( diff( start(DSfreqGR)) == bpDistance) + 1)]$prop,
      method = "pearson")
}

tabGR_genes_all_acf <- lapply(seq_along(chrName), function(x) {
  unlist(mclapply(tabGR_genes_all_dists_list[[x]], function(y) {
    acfDistance(DSfreqGR = c(tabGR_genes_fwd[seqnames(tabGR_genes_fwd) == chrName[x]],
                             tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
                bpDistance = y)
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
})

tabGR_genes_all_random_acf <- mclapply(1:nperm, function(w) {
  lapply(seq_along(chrName), function(x) {
    unlist(lapply(tabGR_genes_all_dists_list[[x]], function(y) {
      acfDistance(DSfreqGR = c(tabGR_genes_fwd_random[[w]][seqnames(tabGR_genes_fwd_random[[w]]) == chrName[x]],
                               tabGR_genes_rev_random[[w]][seqnames(tabGR_genes_rev_random[[w]]) == chrName[x]]),
                  bpDistance = y)
    }))
  })
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)

tabGR_genes_all_acf_permTest_exp <- lapply(seq_along(chrName), function(x) {
  unlist(lapply(seq_along(tabGR_genes_all_dists_list[[x]]), function(y) {
    mean(
      sapply(seq_along(tabGR_genes_all_random_acf), function(w) {
        tabGR_genes_all_random_acf[[w]][[x]][y]
      })
    , na.rm = T)
  }))
})

tabGR_genes_all_acf_permTest_pval <- lapply(seq_along(chrName), function(x) {
  unlist(lapply(seq_along(tabGR_genes_all_dists_list[[x]]), function(y) {
    1 - ( sum(
      sapply(seq_along(tabGR_genes_all_random_acf), function(w) {
        tabGR_genes_all_acf[[x]][y] > tabGR_genes_all_random_acf[[w]][[x]][y]
      })
    , na.rm = T) / nperm )
  }))
})

# Set minimum p-value
for(x in seq_along(chrName)) {
  tabGR_genes_all_acf_permTest_pval[[x]][which(tabGR_genes_all_acf_permTest_pval[[x]] == 0)] <- min_pval
}

tabGR_genes_all_acf_df <- dplyr::bind_rows(lapply(seq_along(tabGR_genes_all_acf), function(x) {
  data.frame(chr = chrName[x],
             distance = tabGR_genes_all_dists_list[[x]],
             acf = tabGR_genes_all_acf[[x]],
             fft = fft(tabGR_genes_all_acf[[x]]),
             pval = -log10(tabGR_genes_all_acf_permTest_pval[[x]]),
             exp = tabGR_genes_all_acf_permTest_pval[[x]])
}))


## Analyse each strand separately
#tabGR_genes_rev <- tabGR_genes[strand(tabGR_genes) == "-"]
#tabGR_genes_rev <- sortSeqlevels(tabGR_genes_rev)
#tabGR_genes_rev <- sort(tabGR_genes_rev, by = ~ seqnames + start + end)
#
## Randomly shuffle methylation proportion values over the cytosine coordinates
#tabGR_genes_rev_random <- mclapply(1:nperm, function(y) {
#  tabGR_genes_rev_chr_list_y <- lapply(seq_along(chrName), function(x) {
#    GRanges(seqnames = seqnames(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
#            ranges = ranges(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
#            strand = strand(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]),
#            prop = sample(tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]]$prop))
#  })
#  do.call(c, as(tabGR_genes_rev_chr_list_y, "GRangesList")) 
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
#
##  sort( unique ( start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == x][2:(length(tabGR_genes_rev[seqnames(tabGR_genes_rev) == x]))] ) -
##                 start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == x][1:(length(tabGR_genes_rev[seqnames(tabGR_genes_rev) == x])-1)] ) ) )
#
## Get context-specific inter-cytosine bp distances that occur more than once
#tabGR_genes_rev_dists_list <- mclapply(seq_along(chrName), function(x) {
#  as.integer( names (
#    which( table(
#      sort( diff( start( tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]] ) ) )
#   ) > 2 )
#  ) )
#}, mc.cores = length(chrName), mc.preschedule = F)
#
## Get the distance that is the lowest common denominator among the chromsomes
#tabGR_genes_rev_dists_min_max <- min(sapply(tabGR_genes_rev_dists_list, function(x) {
#  max(x, na.rm = T)
#}))
#
#tabGR_genes_rev_dists_list <- lapply(tabGR_genes_rev_dists_list, function(x) {
#  x[which(x <= tabGR_genes_rev_dists_min_max)]
#})
#
#acfDistance <- function(DSfreqGR, bpDistance) {
#  cor.test(x = DSfreqGR[ ( which( diff( start(DSfreqGR)) == bpDistance) )]$prop,
#           y = DSfreqGR[ ( which( diff( start(DSfreqGR)) == bpDistance) + 1)]$prop,
#           alternative = "two.sided",
#           method = "pearson")$estimate
#}
#
#tabGR_genes_rev_acf <- lapply(seq_along(chrName), function(x) {
#  unlist(mclapply(tabGR_genes_rev_dists_list[[x]], function(y) {
#    acfDistance(DSfreqGR = tabGR_genes_rev[seqnames(tabGR_genes_rev) == chrName[x]],
#                bpDistance = y)
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
#})
#
#tabGR_genes_rev_random_acf <- mclapply(seq_along(tabGR_genes_rev_random), function(w) {
#  lapply(seq_along(chrName), function(x) {
#    unlist(lapply(tabGR_genes_rev_dists_list[[x]], function(y) {
#      acfDistance(DSfreqGR = tabGR_genes_rev_random[[w]][seqnames(tabGR_genes_rev_random[[w]]) == chrName[x]],
#                  bpDistance = y)
#    }))
#  })
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)
#
#tabGR_genes_rev_acf_permTest_exp <- lapply(seq_along(chrName), function(x) {
#  unlist(lapply(seq_along(tabGR_genes_rev_dists_list[[x]]), function(y) {
#    mean(
#      sapply(seq_along(tabGR_genes_rev_random_acf), function(w) {
#        tabGR_genes_rev_random_acf[[w]][[x]][y]
#      })
#    , na.rm = T)
#  }))
#})
#
#tabGR_genes_rev_acf_permTest_pval <- lapply(seq_along(chrName), function(x) {
#  unlist(lapply(seq_along(tabGR_genes_rev_dists_list[[x]]), function(y) {
#    1 - ( sum(
#      sapply(seq_along(tabGR_genes_rev_random_acf), function(w) {
#        tabGR_genes_rev_acf[[x]][y] > tabGR_genes_rev_random_acf[[w]][[x]][y]
#      })
#    , na.rm = T) / nperm )
#  }))
#})
#
## Set minimum p-value
#for(x in seq_along(chrName)) {
#  tabGR_genes_rev_acf_permTest_pval[[x]][which(tabGR_genes_rev_acf_permTest_pval[[x]] == 0)] <- min_pval
#}
#
#tabGR_genes_rev_acf_df <- dplyr::bind_rows(lapply(seq_along(tabGR_genes_rev_acf), function(x) {
#  data.frame(chr = chrName[x],
#             distance = tabGR_genes_rev_dists_list[[x]],
#             acf = tabGR_genes_rev_acf[[x]],
#             pval = -log10(tabGR_genes_rev_acf_permTest_pval[[x]]),
#             exp = tabGR_genes_rev_acf_permTest_pval[[x]])
#}))


#tabGR_genes_all_acf_df <- rbind(tabGR_genes_fwd_acf_df, tabGR_genes_rev_acf_df)
#tabGR_genes_all_acf_df <- tabGR_genes_all_acf_df[ order(tabGR_genes_all_acf_df[,1], tabGR_genes_all_acf_df[,2], decreasing = F) , ]


chrPlot <- function(dataFrame, xvar, yvar, xlab, ylab, colour) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) +
#  geom_line(colour = colour, size = 1) +
  geom_ma(ma_fun = SMA, n = 6, colour = colour, linetype = 1, size = 2.5) +
  scale_x_continuous(
                     labels = function(x) x) +
  labs(x = xlab,
       y = ylab) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 1.5, colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        panel.grid = element_blank(),
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}


gg_tabGR_genes_all_acf <- chrPlot(dataFrame = tabGR_genes_all_acf_df,
                                   xvar = distance,
                                   yvar = acf,
                                   xlab = bquote("Distance between genic cytosines (bp)"),
                                   ylab = bquote("Observed correlation (m"*.(context)*")"),
                                   colour = "dodgerblue")
gg_tabGR_genes_all_acf <- gg_tabGR_genes_all_acf +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_genes_all_pval <- chrPlot(dataFrame = tabGR_genes_all_acf_df,
                                    xvar = distance,
                                    yvar = pval,
                                    xlab = bquote("Distance between genic cytosines (bp)"),
                                    ylab = bquote("-"*Log[10]*"("*italic(P)*"-value) (m"*.(context)*")"),
                                    colour = "red")
gg_tabGR_genes_all_pval <- gg_tabGR_genes_all_pval +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_genes_all_exp <- chrPlot(dataFrame = tabGR_genes_all_acf_df,
                                   xvar = distance,
                                   yvar = exp,
                                   xlab = bquote("Distance between genic cytosines (bp)"),
                                   ylab = bquote("Mean permuted correlation (m"*.(context)*")"),
                                   colour = "lightseagreen")
gg_tabGR_genes_all_exp <- gg_tabGR_genes_all_exp +
  facet_grid(cols = vars(chr), scales = "free_x")

#gg_tabGR_genes_all_fft <- chrPlot(dataFrame = tabGR_genes_all_acf_df,
#                                   xvar = distance,
#                                   yvar = fft,
#                                   xlab = bquote("Distance between genic cytosines (bp)"),
#                                   ylab = bquote("FFT (m"*.(context)*")"),
#                                   colour = "lightseagreen")
#gg_tabGR_genes_all_fft <- gg_tabGR_genes_all_fft +
#  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_all_list <- list(
                        gg_tabGR_genes_all_acf,
                        gg_tabGR_genes_all_pval,
#                        gg_tabGR_genes_all_fft,
                        gg_tabGR_genes_all_exp
                       )
gg_cow_all <- plot_grid(plotlist = gg_cow_all_list,
                        labels = c("AUTO"), label_size = 30,
                        align = "hv",
                        axis = "l",
                        nrow = length(gg_cow_all_list), ncol = 1)

ggsave(paste0(plotDir,
              sampleName, "_MappedOn_", refbase, "_", context,
              "_autocorrelation_all_genes_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow_all,
       height = 5*length(gg_cow_all_list), width = 10*length(chrName), limitsize = F)


#gg_tabGR_genes_rev_acf <- chrPlot(dataFrame = tabGR_genes_rev_acf_df,
#                                   xvar = distance,
#                                   yvar = acf,
#                                   xlab = bquote("Distance between genic cytosines (bp)"),
#                                   ylab = bquote("Observed correlation (m"*.(context)*")"),
#                                   colour = "dodgerblue")
#gg_tabGR_genes_rev_acf <- gg_tabGR_genes_rev_acf +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_tabGR_genes_rev_pval <- chrPlot(dataFrame = tabGR_genes_rev_acf_df,
#                                    xvar = distance,
#                                    yvar = pval,
#                                    xlab = bquote("Distance between genic cytosines (bp)"),
#                                    ylab = bquote("-"*Log[10]*"("*italic(P)*"-value) (m"*.(context)*")"),
#                                    colour = "red")
#gg_tabGR_genes_rev_pval <- gg_tabGR_genes_rev_pval +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_tabGR_genes_rev_exp <- chrPlot(dataFrame = tabGR_genes_rev_acf_df,
#                                   xvar = distance,
#                                   yvar = exp,
#                                   xlab = bquote("Distance between genic cytosines (bp)"),
#                                   ylab = bquote("Mean permuted correlation (m"*.(context)*")"),
#                                   colour = "lightseagreen")
#gg_tabGR_genes_rev_exp <- gg_tabGR_genes_rev_exp +
#  facet_grid(cols = vars(chr), scales = "free_x")
#
#gg_cow_rev_list <- list(
#                        gg_tabGR_genes_rev_acf,
#                        gg_tabGR_genes_rev_pval,
#                        gg_tabGR_genes_rev_exp
#                       )
#gg_cow_rev <- plot_grid(plotlist = gg_cow_rev_list,
#                        labels = c("AUTO"), label_size = 30,
#                        align = "hv",
#                        axis = "l",
#                        nrow = length(gg_cow_rev_list), ncol = 1)
#
#ggsave(paste0(plotDir,
#              sampleName, "_MappedOn_", refbase, "_", context,
#              "_autocorrelation_rev_genes_", paste0(chrName, collapse = "_"),
#              ".pdf"),
#       plot = gg_cow_rev,
#       height = 5*length(gg_cow_rev_list), width = 10*length(chrName), limitsize = F)

