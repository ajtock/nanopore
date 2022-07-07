#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 07.07.2022

# Calculate and plot metaprofiles of DeepSignal DNA methylation
# (feature windowed means and 95% confidence intervals, CIs)
# for NRZgenes, matched randomly positioned loci, all other genes (notNRZgenes),
# and matched randomly positioned loci

# Usage:
# conda activate R-4.0.0
# ./NRZgenes_notNRZgenes_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 2000 2000 2kb 10 10bp '0.02,0.96' 'Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CpG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHG,Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH' 'nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth,nanopore/t2t-col.20210610/deepsignal_DNAmeth' 'CG,CHG,CHH' 'dodgerblue4,dodgerblue1,cyan2' t2t-col.20210610
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
  NRZgenesNamePlot <- "All NRZ genes"
  NRZranLocNamePlot <- "All NRZ ranLoc"
  notNRZgenesNamePlot <- "All not NRZ genes"
  notNRZranLocNamePlot <- "All not NRZ ranLoc"
} else {
  NRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " NRZ genes")
  NRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " NRZ ranLoc")
  notNRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " not NRZ genes")
  notNRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " not NRZ ranLoc")
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
# NRZgene
ChIP_NRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_NRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If NRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding NRZgene coverage matrices
ChIP_NRZgeneMats <- mclapply(seq_along(ChIP_NRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_NRZgeneMats[[x]])
  } else {
    ChIP_NRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_NRZgeneMats))

## ChIP
# NRZranLoc
ChIP_NRZranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_NRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If NRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding NRZranLoc coverage matrices
ChIP_NRZranLocMats <- mclapply(seq_along(ChIP_NRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_NRZranLocMats[[x]])
  } else {
    ChIP_NRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_NRZranLocMats))

## ChIP
# notNRZgene
ChIP_notNRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_notNRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If notNRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding notNRZgene coverage matrices
ChIP_notNRZgeneMats <- mclapply(seq_along(ChIP_notNRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_notNRZgeneMats[[x]])
  } else {
    ChIP_notNRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_notNRZgeneMats))

## ChIP
# notNRZranLoc
ChIP_notNRZranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_notNRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If notNRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding notNRZranLoc coverage matrices
ChIP_notNRZranLocMats <- mclapply(seq_along(ChIP_notNRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_notNRZranLocMats[[x]])
  } else {
    ChIP_notNRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_notNRZranLocMats))


# ChIP
# Add column names
for(x in seq_along(ChIP_NRZgeneMats)) {
  colnames(ChIP_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                         paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_notNRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_notNRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                            paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                            paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_NRZgeneMats), function(x) {
  list(
       # NRZgenes
       ChIP_NRZgeneMats[[x]],
       # NRZranLoc
       ChIP_NRZranLocMats[[x]],
       # notNRZgenes
       ChIP_notNRZgeneMats[[x]],
       # notNRZranLoc
       ChIP_notNRZranLocMats[[x]]
      )
}, mc.cores = length(ChIP_NRZgeneMats))

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
NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
notNRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
notNRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
names(NRZgenesTmp) <- ChIPNamesPlot
names(NRZranLocTmp) <- ChIPNamesPlot
names(notNRZgenesTmp) <- ChIPNamesPlot
names(notNRZranLocTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(NRZgenesTmp, .id = "libName"),
  bind_rows(NRZranLocTmp, .id = "libName"),
  bind_rows(notNRZgenesTmp, .id = "libName"),
  bind_rows(notNRZranLocTmp, .id = "libName")
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

## NRZgenes
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
ggtitle(bquote(.(NRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## NRZranLoc
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
ggtitle(bquote(.(NRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## notNRZgenes
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
ggtitle(bquote(.(notNRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## notNRZranLoc
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
ggtitle(bquote(.(notNRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
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
              "_NRZgenes_NRZranLoc_notNRZgenes_notNRZranLoc_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)
