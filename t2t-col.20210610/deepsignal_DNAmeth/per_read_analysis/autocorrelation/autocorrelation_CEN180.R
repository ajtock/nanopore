#!/applications/R/R-4.0.0/bin/Rscript

# Analysis:
# 1. Genome-wide autocorrelation of context-specific cytosine methylation status
# at increasing physical distances (e.g., 1 to 10,000 nucleotides)

# Usage on hydrogen node7:
# csmit -m 200G -c 47 "/applications/R/R-4.0.0/bin/Rscript autocorrelation_CEN180.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 10000 10000 CpG 1e4 1.00 Chr1,Chr2,Chr3,Chr4,Chr5 200 CEN180"

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#nperm <- 1e3
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#maxDist <- 200
#featName <- "CEN180"
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
featName <- args[10]
min_pval <- 1 - ( (nperm - 1) / nperm )

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggthemes)
library(tidyquant)
library(fastmatch)

# Define more efficient matching function %fin%
`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
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

# Read in feature annotation
if(featName == "CEN180") {
  CEN180 <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                              "/annotation/CEN180/CEN180_in_", refbase,
                              "_", paste0(chrName, collapse = "_"), ".bed"),
                       header = F)
  colnames(CEN180) <- c("chr", "start0based", "end", "name", "score", "strand", "HORlengthsSum", "HORcount")
  CEN180GR <- GRanges(seqnames = CEN180$chr,
                      ranges = IRanges(start = CEN180$start0based+1,
                                       end = CEN180$end),
                      strand = CEN180$strand)
} else if(featName == "gene") {
  genes <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                             "/annotation/genes/", refbase, "_representative_mRNA",
                             "_", paste0(chrName, collapse = "_"), ".bed"),
                      header = F)
  colnames(genes) <- c("chr", "start0based", "end", "name", "score", "strand")
  genesGR <- GRanges(seqnames = genes$chr,
                     ranges = IRanges(start = genes$start0based+1,
                                      end = genes$end),
                     strand = genes$strand)
} else if(featName == "GYPSY") {
  gypsy <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                             "/annotation/TEs_EDTA/", refbase, "_TEs_Gypsy_LTR",
                             "_", paste0(chrName, collapse = "_"), ".bed"),
                      header = F)
  colnames(gypsy) <- c("chr", "start0based", "end", "name", "score", "strand")
  gypsyGR <- GRanges(seqnames = gypsy$chr,
                     ranges = IRanges(start = gypsy$start0based+1,
                                      end = gypsy$end),
                     strand = gypsy$strand)
} else {
  stop(print("featName not one of CEN180, gene or GYPSY"))
}


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

# Find context-specific cytosine sites that overlap CEN180
fOverlaps_tab_CEN180 <- findOverlaps(query = tabGR,
                                     subject = CEN180GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = T)
tabGR_CEN180 <- tabGR[unique(queryHits(fOverlaps_tab_CEN180))]


# Analyse each strand separately
tabGR_CEN180_fwd <- tabGR_CEN180[strand(tabGR_CEN180) == "+"]
tabGR_CEN180_fwd <- sortSeqlevels(tabGR_CEN180_fwd)
tabGR_CEN180_fwd <- sort(tabGR_CEN180_fwd, by = ~ seqnames + start + end)

# Randomly shuffle methylation proportion values over the cytosine coordinates
tabGR_CEN180_fwd_random <- mclapply(1:nperm, function(y) {
  tabGR_CEN180_fwd_chr_list_y <- lapply(seq_along(chrName), function(x) {
    GRanges(seqnames = seqnames(tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]]),
            ranges = ranges(tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]]),
            strand = strand(tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]]),
            prop = sample(tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]]$prop))
  })
  do.call(c, as(tabGR_CEN180_fwd_chr_list_y, "GRangesList")) 
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

# For each chromosome, get context-specific inter-cytosine distances
tabGR_CEN180_fwd_dists_list <- lapply(seq_along(chrName), function(x) {
  tabGR_CEN180_fwd_chr <- tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]]
  mclapply(seq_along(tabGR_CEN180_fwd_chr), function(y) {
    start(tabGR_CEN180_fwd_chr) - start(tabGR_CEN180_fwd_chr[y])
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
})

# NOT SURE DOUBLE-mclapplying IS A GOOD IDEA HERE - NEED TO REDUCE CPUpc (E.G. TO 0.50 ?)
# TO AVOID EXCEEDING RAM LIMIT
# For each chromosome, get row (GRanges) indices for each inter-cytosine distance
tabGR_CEN180_fwd_dists_bool_list <- mclapply(seq_along(chrName), function(x) {
  lapply(1:maxDist, function(z) {
    unlist(mclapply(seq_along(tabGR_CEN180_fwd_dists_list[[x]]), function(y) {
      z %fin% tabGR_CEN180_fwd_dists_list[[x]][[y]]
    }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
  })
}, mc.cores = length(chrName), mc.preschedule = F)

# For each chromosome, get inter-cytosine distances for which there are > 2
# data points to allow correlation coefficient calculation
tabGR_CEN180_fwd_dists_bool_list_gt2 <- unlist( lapply(seq_along(chrName), function(x) {
  which(sapply(tabGR_CEN180_fwd_dists_bool_list[[x]], function(z) {
    sum(z) > 2
  }))
}) )


# Analyse each strand separately
tabGR_CEN180_rev <- tabGR_CEN180[strand(tabGR_CEN180) == "-"]
tabGR_CEN180_rev <- sortSeqlevels(tabGR_CEN180_rev)
tabGR_CEN180_rev <- sort(tabGR_CEN180_rev, by = ~ seqnames + start + end)

# Randomly shuffle methylation proportion values over the cytosine coordinates
tabGR_CEN180_rev_random <- mclapply(1:nperm, function(y) {
  tabGR_CEN180_rev_chr_list_y <- lapply(seq_along(chrName), function(x) {
    GRanges(seqnames = seqnames(tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]]),
            ranges = ranges(tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]]),
            strand = strand(tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]]),
            prop = sample(tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]]$prop))
  })
  do.call(c, as(tabGR_CEN180_rev_chr_list_y, "GRangesList")) 
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)

# For each chromosome, get context-specific inter-cytosine distances
tabGR_CEN180_rev_dists_list <- lapply(seq_along(chrName), function(x) {
  tabGR_CEN180_rev_chr <- tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]]
  mclapply(seq_along(tabGR_CEN180_rev_chr), function(y) {
    start(tabGR_CEN180_rev_chr) - start(tabGR_CEN180_rev_chr[y])
  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
})

# NOT SURE DOUBLE-mclapplying IS A GOOD IDEA HERE - NEED TO REDUCE CPUpc (E.G. TO 0.50 ?)
# TO AVOID EXCEEDING RAM LIMIT
# For each chromosome, get row (GRanges) indices for each inter-cytosine distance
tabGR_CEN180_rev_dists_bool_list <- mclapply(seq_along(chrName), function(x) {
  lapply(1:maxDist, function(z) {
    unlist(mclapply(seq_along(tabGR_CEN180_rev_dists_list[[x]]), function(y) {
      z %fin% tabGR_CEN180_rev_dists_list[[x]][[y]]
    }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
  })
}, mc.cores = length(chrName), mc.preschedule = F)

# For each chromosome, get inter-cytosine distances for which there are > 2
# data points to allow correlation coefficient calculation
tabGR_CEN180_rev_dists_bool_list_gt2 <- unlist( lapply(seq_along(chrName), function(x) {
  which(sapply(tabGR_CEN180_rev_dists_bool_list[[x]], function(z) {
    sum(z) > 2
  }))
}) )


# Get inter-cytosine distances that are shared across all chromosomes and
# across fwd and rev analyses
tabGR_CEN180_all_dists_bool_list_gt2 <- c(tabGR_CEN180_fwd_dists_bool_list_gt2,
                                          tabGR_CEN180_rev_dists_bool_list_gt2)
tabGR_CEN180_all_dists_bool_list_gt2 <- as.integer( names(
    which( table( tabGR_CEN180_all_dists_bool_list_gt2 ) == length(chrName) * 2 )
) )


# Define "autocorrelation" function 
acfDistance <- function(fwd_DSfreqGR, rev_DSfreqGR, rev_dists_bool_list, fwd_dists_bool_list, bpDistance) { 
  cor(x = c( fwd_DSfreqGR[ ( which( fwd_dists_bool_list[[bpDistance]] ) ) ]$prop,
             rev_DSfreqGR[ ( which( rev_dists_bool_list[[bpDistance]] ) ) ]$prop ),
      y = c( fwd_DSfreqGR[ start(fwd_DSfreqGR) %in%
                             ( start(fwd_DSfreqGR[ ( which( fwd_dists_bool_list[[bpDistance]] ) ) ]) + bpDistance ) &
                           strand(fwd_DSfreqGR) ==
                             ( unique( as.character( strand(fwd_DSfreqGR[ ( which( fwd_dists_bool_list[[bpDistance]] ) ) ]) ) ) ) ]$prop,
             rev_DSfreqGR[ start(rev_DSfreqGR) %in%
                             ( start(rev_DSfreqGR[ ( which( rev_dists_bool_list[[bpDistance]] ) ) ]) + bpDistance ) & 
                           strand(rev_DSfreqGR) ==
                             ( unique( as.character( strand(rev_DSfreqGR[ ( which( rev_dists_bool_list[[bpDistance]] ) ) ]) ) ) ) ]$prop ),
      method = "pearson")
}

# Apply "autocorrelation" function
tabGR_CEN180_all_acf <- mclapply(seq_along(chrName), function(x) {
 sapply(tabGR_CEN180_all_dists_bool_list_gt2, function(y) {
    acfDistance(fwd_DSfreqGR = tabGR_CEN180_fwd[seqnames(tabGR_CEN180_fwd) == chrName[x]],
                rev_DSfreqGR = tabGR_CEN180_rev[seqnames(tabGR_CEN180_rev) == chrName[x]],
                fwd_dists_bool_list = tabGR_CEN180_fwd_dists_bool_list[[x]],
                rev_dists_bool_list = tabGR_CEN180_rev_dists_bool_list[[x]],
                bpDistance = y)
  })
}, mc.cores = length(chrName), mc.preschedule = F)

tabGR_CEN180_all_random_acf <- mclapply(1:nperm, function(w) {
  lapply(seq_along(chrName), function(x) {
    sapply(tabGR_CEN180_all_dists_bool_list_gt2, function(y) {
      acfDistance(fwd_DSfreqGR = tabGR_CEN180_fwd_random[[w]][seqnames(tabGR_CEN180_fwd_random[[w]]) == chrName[x]],
                  rev_DSfreqGR = tabGR_CEN180_rev_random[[w]][seqnames(tabGR_CEN180_rev_random[[w]]) == chrName[x]],
                  fwd_dists_bool_list = tabGR_CEN180_fwd_dists_bool_list[[x]],
                  rev_dists_bool_list = tabGR_CEN180_rev_dists_bool_list[[x]],
                  bpDistance = y)
    })
  })
}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)

tabGR_CEN180_all_acf_permTest_exp <- lapply(seq_along(chrName), function(x) {
  sapply(seq_along(tabGR_CEN180_all_dists_bool_list_gt2), function(y) {
    mean(
      sapply(seq_along(tabGR_CEN180_all_random_acf), function(w) {
        tabGR_CEN180_all_random_acf[[w]][[x]][y]
      })
    , na.rm = T)
  })
})

tabGR_CEN180_all_acf_permTest_pval <- lapply(seq_along(chrName), function(x) {
  sapply(seq_along(tabGR_CEN180_all_dists_bool_list_gt2), function(y) {
    1 - ( sum(
      sapply(seq_along(tabGR_CEN180_all_random_acf), function(w) {
        tabGR_CEN180_all_acf[[x]][y] > tabGR_CEN180_all_random_acf[[w]][[x]][y]
      })
    , na.rm = T) / nperm )
  })
})

# Set minimum p-value
for(x in seq_along(chrName)) {
  tabGR_CEN180_all_acf_permTest_pval[[x]][which(tabGR_CEN180_all_acf_permTest_pval[[x]] == 0)] <- min_pval
}

tabGR_CEN180_all_acf_df <- dplyr::bind_rows(lapply(seq_along(tabGR_CEN180_all_acf), function(x) {
  data.frame(chr = chrName[x],
             distance = tabGR_CEN180_all_dists_bool_list_gt2,
             acf = tabGR_CEN180_all_acf[[x]],
             pval = -log10(tabGR_CEN180_all_acf_permTest_pval[[x]]),
             exp = tabGR_CEN180_all_acf_permTest_pval[[x]])
}))

write.table(tabGR_CEN180_all_acf_df,
            paste0(sampleName, "_MappedOn_", refbase, "_", context,
                   "_all_autocorrelation_", featName, "_", paste0(chrName, collapse = "_"),
                   ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# Chromosome-scale plotting function
chrPlot <- function(dataFrame, xvar, yvar, xlab, ylab, colour) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = aes(x = !!xvar,
                       y = !!yvar)) +
#  geom_line(colour = colour, size = 1) +
  geom_ma(ma_fun = SMA, n = 6, colour = colour, linetype = 1, size = 2) +
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
        strip.text.x = element_text(size = 30, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"))
}


gg_tabGR_CEN180_all_acf <- chrPlot(dataFrame = tabGR_CEN180_all_acf_df,
                                   xvar = distance,
                                   yvar = acf,
                                   xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                   ylab = bquote("Observed correlation (m"*.(context)*")"),
                                   colour = "dodgerblue")
gg_tabGR_CEN180_all_acf <- gg_tabGR_CEN180_all_acf +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_CEN180_all_pval <- chrPlot(dataFrame = tabGR_CEN180_all_acf_df,
                                    xvar = distance,
                                    yvar = pval,
                                    xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                    ylab = bquote("-"*Log[10]*"("*italic(P)*"-value) (m"*.(context)*")"),
                                    colour = "red")
gg_tabGR_CEN180_all_pval <- gg_tabGR_CEN180_all_pval +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_tabGR_CEN180_all_exp <- chrPlot(dataFrame = tabGR_CEN180_all_acf_df,
                                   xvar = distance,
                                   yvar = exp,
                                   xlab = bquote("Distance between "*italic(.(featName))*" cytosines (bp)"),
                                   ylab = bquote("Mean permuted correlation (m"*.(context)*")"),
                                   colour = "lightseagreen")
gg_tabGR_CEN180_all_exp <- gg_tabGR_CEN180_all_exp +
  facet_grid(cols = vars(chr), scales = "free_x")

gg_cow_all_list <- list(
                        gg_tabGR_CEN180_all_acf,
                        gg_tabGR_CEN180_all_pval,
                        gg_tabGR_CEN180_all_exp
                       )
gg_cow_all <- plot_grid(plotlist = gg_cow_all_list,
                        labels = c("AUTO"), label_size = 30,
                        align = "hv",
                        axis = "l",
                        nrow = length(gg_cow_all_list), ncol = 1)

ggsave(paste0(plotDir,
              sampleName, "_MappedOn_", refbase, "_", context,
              "_all_autocorrelation_", featName, "_", paste0(chrName, collapse = "_"),
              ".pdf"),
       plot = gg_cow_all,
       height = 5*length(gg_cow_all_list), width = 10*length(chrName), limitsize = F)




#tabGR_CEN180_fwd_chr_dist_bool_list <- mclapply(seq_along(tabGR_CEN180_fwd_chr_dist_list), function(y) {
#  sapply(1:maxDist, function(z) {
#    z %fin% tabGR_CEN180_fwd_chr_dist_list[[y]]
#  })
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T)
#
#tabGR_CEN180_fwd_chr_dist_bool_list <- lapply(seq_alonglapply(1:maxDist, function(z) {
#  unlist(mclapply(seq_along(tabGR_CEN180_fwd_chr_dist_list), function(y) {
#    z %fin% tabGR_CEN180_fwd_chr_dist_list[[y]]
#  }, mc.cores = round(detectCores()*CPUpc), mc.preschedule = T))
#})
#
#tabGR_CEN180_fwd_chr_dist_bool_list <- mclapply(1:maxDist, function(z) {
#  sapply(seq_along(tabGR_CEN180_fwd_chr_dist_list), function(y) {
#    z %fin% tabGR_CEN180_fwd_chr_dist_list[[y]]
#  })
#}, mc.cores = round(detectCores()*CPUpc), mc.preschedule = F)
