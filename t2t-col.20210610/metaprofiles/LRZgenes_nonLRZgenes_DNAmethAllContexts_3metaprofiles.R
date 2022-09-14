#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 07.07.2022

# Calculate and plot metaprofiles of DeepSignal DNA methylation
# (feature windowed means and 95% confidence intervals, CIs)
# for LRZgenes, matched randomly positioned loci, all other genes (nonLRZgenes),
# and matched randomly positioned loci

# Usage:
# conda activate R-4.0.0
# ./LRZgenes_nonLRZgenes_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 2000 2000 2kb 10 10bp '0.02,0.96' 'Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CpG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH' 'nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth' 'CG,CHG,CHH' 'dodgerblue4,dodgerblue1,cyan2' t2t-col.20210610
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CpG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth",
#                                split = ","))
#ChIPNamesPlot <- unlist(strsplit("mCG,mCHG,mCHH",
#                                 split = ","))
#ChIPColours <- unlist(strsplit("dodgerblue4,dodgerblue1,cyan2",
#                               split = ","))
#refbase <- "t2t-col.20210610"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
bodyLength <- as.numeric(args[2])
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
binSize <- as.numeric(args[5])
binName <- args[6]
legendPos <- as.numeric(unlist(strsplit(args[7],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[8],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[9],
                                split = ","))
ChIPNamesPlot <- unlist(strsplit(args[10],
                                 split = ","))
ChIPColours <- unlist(strsplit(args[11],
                               split = ","))
refbase <- args[12]

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

if(length(chrName) == 5) {
  LRZgenesNamePlot <- "All LRZ genes"
  LRZranLocNamePlot <- "All LRZ ranLoc"
  nonLRZgenesNamePlot <- "All non-LRZ genes"
  nonLRZranLocNamePlot <- "All non-LRZ ranLoc"
} else {
  LRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ genes")
  LRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ ranLoc")
  nonLRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " non-LRZ genes")
  nonLRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " non-LRZ ranLoc")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

# Load feature matrices for each dataset
ChIPDirs <- sapply(seq_along(ChIPNamesDir), function(x) {
  if(grepl("deepsignal", ChIPNamesDir[x])) {
    paste0("/home/ajt200/analysis/",
           ChIPNamesDir[x],
           "/")
  } else {
    paste0("/home/ajt200/analysis/",
           ChIPNamesDir[x],
           "/coverage/")
  }
})

## ChIP
# LRZgene
ChIP_LRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_LRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If LRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding LRZgene coverage matrices
ChIP_LRZgeneMats <- mclapply(seq_along(ChIP_LRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_LRZgeneMats[[x]])
  } else {
    ChIP_LRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_LRZgeneMats))

## ChIP
# LRZranLoc
ChIP_LRZranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_LRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If LRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding LRZranLoc coverage matrices
ChIP_LRZranLocMats <- mclapply(seq_along(ChIP_LRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_LRZranLocMats[[x]])
  } else {
    ChIP_LRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_LRZranLocMats))

## ChIP
# nonLRZgene
ChIP_nonLRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_nonLRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If nonLRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding nonLRZgene coverage matrices
ChIP_nonLRZgeneMats <- mclapply(seq_along(ChIP_nonLRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_nonLRZgeneMats[[x]])
  } else {
    ChIP_nonLRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_nonLRZgeneMats))

## ChIP
# nonLRZranLoc
ChIP_nonLRZranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_nonLRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If nonLRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding nonLRZranLoc coverage matrices
ChIP_nonLRZranLocMats <- mclapply(seq_along(ChIP_nonLRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_nonLRZranLocMats[[x]])
  } else {
    ChIP_nonLRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_nonLRZranLocMats))


# ChIP
# Add column names
for(x in seq_along(ChIP_LRZgeneMats)) {
  colnames(ChIP_LRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_LRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                         paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_nonLRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_nonLRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                            paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                            paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_LRZgeneMats), function(x) {
  list(
       # LRZgenes
       ChIP_LRZgeneMats[[x]],
       # LRZranLoc
       ChIP_LRZranLocMats[[x]],
       # nonLRZgenes
       ChIP_nonLRZgeneMats[[x]],
       # nonLRZranLoc
       ChIP_nonLRZranLocMats[[x]]
      )
}, mc.cores = length(ChIP_LRZgeneMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(ChIP_mats[[x]][[y]]),
               t(ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    tidyDFfeature_list_ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]]$window,
                                                       levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    summaryDFfeature_list_ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]]$window,
                                                          levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
    summaryDFfeature_list_ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]])[1])
    summaryDFfeature_list_ChIP[[x]][[y]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_ChIP into
# a list of single data.frames containing all meta-profiles for plotting
LRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
LRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
nonLRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
nonLRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
names(LRZgenesTmp) <- ChIPNamesPlot
names(LRZranLocTmp) <- ChIPNamesPlot
names(nonLRZgenesTmp) <- ChIPNamesPlot
names(nonLRZranLocTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(LRZgenesTmp, .id = "libName"),
  bind_rows(LRZranLocTmp, .id = "libName"),
  bind_rows(nonLRZgenesTmp, .id = "libName"),
  bind_rows(nonLRZranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                               levels = ChIPNamesPlot)
}

# Define y-axis limits
#ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
#                   summaryDFfeature_ChIP[[2]]$CI_lower,
#                   summaryDFfeature_ChIP[[3]]$CI_lower,
#                   summaryDFfeature_ChIP[[4]]$CI_lower),
#                 na.rm = T)
#ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
#                   summaryDFfeature_ChIP[[2]]$CI_upper,
#                   summaryDFfeature_ChIP[[3]]$CI_upper,
#                   summaryDFfeature_ChIP[[4]]$CI_upper),
#                 na.rm = T)
ymin_ChIP <- 0
ymax_ChIP <- 100

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 22)))
})

# Plot average profiles with 95% CI ribbon

## LRZgenes
summaryDFfeature <- summaryDFfeature_ChIP[[1]]
ggObj1_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(LRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## LRZranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[2]]
ggObj2_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(LRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonLRZgenes
summaryDFfeature <- summaryDFfeature_ChIP[[3]]
ggObj3_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonLRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonLRZranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[4]]
ggObj4_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonLRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4
                                                      ))
ggsave(paste0(plotDir,
              "DNAmethAllContexts_",
              gsub("_MappedOn_.+", "", ChIPNames)[1],
              "_avgProfiles_around",
              "_LRZgenes_LRZranLoc_nonLRZgenes_nonLRZranLoc_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)
