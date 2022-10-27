#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 27.10.2022

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for all CEN180 sequences, CENAthila, Ty3, and randomly positioned loci

# Usage:
# conda activate R-4.0.0
# ./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_H3K9me2_Rep1_ChIP' 'WT_H3K9me2_Rep1_input' '170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610' '170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610' 'H3K9me2' 'H3K9me2 input' 'dodgerblue1' 'dodgerblue1' 't2t-col.20210610'
# ./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_REC8_HA_Rep2_ChIP' 'WT_REC8_Myc_Rep1_input' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8' 'REC8 input' 'firebrick1' 'firebrick1' 't2t-col.20210610'
# ./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_ASY1_Rep1_ChIP' 'WT_REC8_Myc_Rep1_input' '20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'ASY1' 'REC8 input' 'darkgreen' 'firebrick1' 't2t-col.20210610'
# ./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_SPO11oligos_Rep1' 'WT_gDNA_Rep1_R1' '160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_t2t-col.20210610' '150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_t2t-col.20210610' 'SPO11-1-oligos' 'gDNA' 'magenta3' 'magenta3' 't2t-col.20210610'
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
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
#ChIPNames <- unlist(strsplit("WT_H3K9me2_Rep1_ChIP",
#                             split = ","))
#controlNames <- unlist(strsplit("WT_H3K9me2_Rep1_input",
#                                split = ","))
#ChIPNamesDir <- unlist(strsplit("170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610",
#                                split = ","))
#controlNamesDir <- unlist(strsplit("170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610",
#                                   split = ","))
#ChIPNamesPlot <- unlist(strsplit("H3K9me2",
#                                 split = ","))
#controlNamesPlot <- unlist(strsplit("H3K9me2 input",
#                                    split = ","))
#ChIPColours <- unlist(strsplit("dodgerblue1",
#                               split = ","))
#controlColours <- unlist(strsplit("dodgerblue1",
#                                  split = ","))
#refbase <- unlist(strsplit("t2t-col.20210610",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[4])
flankName <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
legendPos <- as.numeric(unlist(strsplit(args[8],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[9],
                             split = ","))
controlNames <- unlist(strsplit(args[10],
                                split = ","))
ChIPNamesDir <- unlist(strsplit(args[11],
                                split = ","))
controlNamesDir <- unlist(strsplit(args[12],
                                   split = ","))
ChIPNamesPlot <- unlist(strsplit(args[13],
                                 split = ","))
controlNamesPlot <- unlist(strsplit(args[14],
                                    split = ","))
ChIPColours <- unlist(strsplit(args[15],
                               split = ","))
controlColours <- unlist(strsplit(args[16],
                                  split = ","))
refbase <- unlist(strsplit(args[17],
                           split = ","))

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
  NRZgenesNamePlot <- "All NRZ genes"
  NRZranLocNamePlot <- "All NRZ ranLoc"
  LRZ_NRZgenesNamePlot <- "All LRZ+NRZ genes"
  LRZ_NRZranLocNamePlot <- "All LRZ+NRZ ranLoc"
  armgenesNamePlot <- "All arm genes"
  armranLocNamePlot <- "All arm ranLoc"
} else {
  LRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ genes")
  LRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ ranLoc")
  NRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " NRZ genes")
  NRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " NRZ ranLoc")
  LRZ_NRZgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ+NRZ genes")
  LRZ_NRZranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " LRZ+NRZ ranLoc")
  armgenesNamePlot <- paste0(paste0(chrName, collapse = ","), " arm genes")
  armranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " arm ranLoc")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/")
})


## ChIP
# LRZgene
ChIP_LRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_LRZgenes_in_",
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
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_LRZgenes_in_",
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
# NRZgene
ChIP_NRZgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_NRZgenes_in_",
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
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_NRZgenes_in_",
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
# LRZ_NRZgene
ChIP_LRZ_NRZgeneMats <- mclapply(seq_along(ChIP_LRZgeneMats), function(x) {
  rbind(ChIP_LRZgeneMats[[x]], ChIP_NRZgeneMats[[x]])
})

## ChIP
# LRZ_NRZranLoc
ChIP_LRZ_NRZranLocMats <- mclapply(seq_along(ChIP_LRZranLocMats), function(x) {
  rbind(ChIP_LRZranLocMats[[x]], ChIP_NRZranLocMats[[x]])
})


## ChIP
# armgene
ChIP_armgeneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_armgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If armgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding armgene coverage matrices
ChIP_armgeneMats <- mclapply(seq_along(ChIP_armgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_armgeneMats[[x]])
  } else {
    ChIP_armgeneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_armgeneMats))

## ChIP
# armranLoc
ChIP_armranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_armgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If armranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding armranLoc coverage matrices
ChIP_armranLocMats <- mclapply(seq_along(ChIP_armranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_armranLocMats[[x]])
  } else {
    ChIP_armranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_armranLocMats))


## control
# LRZgene
control_LRZgeneMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_LRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If LRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding LRZgene coverage matrices
control_LRZgeneMats <- mclapply(seq_along(control_LRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_LRZgeneMats[[x]])
  } else {
    control_LRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(control_LRZgeneMats))

## control
# LRZranLoc
control_LRZranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_LRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If LRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding LRZranLoc coverage matrices
control_LRZranLocMats <- mclapply(seq_along(control_LRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_LRZranLocMats[[x]])
  } else {
    control_LRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_LRZranLocMats))


## control
# NRZgene
control_NRZgeneMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_NRZgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If NRZgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding NRZgene coverage matrices
control_NRZgeneMats <- mclapply(seq_along(control_NRZgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_NRZgeneMats[[x]])
  } else {
    control_NRZgeneMats[[x]][[1]]
  }
}, mc.cores = length(control_NRZgeneMats))

## control
# NRZranLoc
control_NRZranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_NRZgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If NRZranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding NRZranLoc coverage matrices
control_NRZranLocMats <- mclapply(seq_along(control_NRZranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_NRZranLocMats[[x]])
  } else {
    control_NRZranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_NRZranLocMats))


## control
# LRZ_NRZgene
control_LRZ_NRZgeneMats <- mclapply(seq_along(control_LRZgeneMats), function(x) {
  rbind(control_LRZgeneMats[[x]], control_NRZgeneMats[[x]])
})

## control
# LRZ_NRZranLoc
control_LRZ_NRZranLocMats <- mclapply(seq_along(control_LRZranLocMats), function(x) {
  rbind(control_LRZranLocMats[[x]], control_NRZranLocMats[[x]])
})


## control
# armgene
control_armgeneMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_armgenes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If armgenes from multiple chromosomes are to be analysed,
# concatenate the corresponding armgene coverage matrices
control_armgeneMats <- mclapply(seq_along(control_armgeneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_armgeneMats[[x]])
  } else {
    control_armgeneMats[[x]][[1]]
  }
}, mc.cores = length(control_armgeneMats))

## control
# armranLoc
control_armranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "LRZ_NRZgeneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_armgenes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If armranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding armranLoc coverage matrices
control_armranLocMats <- mclapply(seq_along(control_armranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_armranLocMats[[x]])
  } else {
    control_armranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_armranLocMats))



# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# LRZgene
log2ChIP_LRZgeneMats <- mclapply(seq_along(ChIP_LRZgeneMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_LRZgeneMats[[x]]+1)/(control_LRZgeneMats[[x]]+1))
}, mc.cores = length(ChIP_LRZgeneMats))

# LRZranLoc
log2ChIP_LRZranLocMats <- mclapply(seq_along(ChIP_LRZranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_LRZranLocMats[[x]]+1)/(control_LRZranLocMats[[x]]+1))
}, mc.cores = length(ChIP_LRZranLocMats))

# NRZgene
log2ChIP_NRZgeneMats <- mclapply(seq_along(ChIP_NRZgeneMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_NRZgeneMats[[x]]+1)/(control_NRZgeneMats[[x]]+1))
}, mc.cores = length(ChIP_NRZgeneMats))

# NRZranLoc
log2ChIP_NRZranLocMats <- mclapply(seq_along(ChIP_NRZranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_NRZranLocMats[[x]]+1)/(control_NRZranLocMats[[x]]+1))
}, mc.cores = length(ChIP_NRZranLocMats))

# LRZ_NRZgene
log2ChIP_LRZ_NRZgeneMats <- mclapply(seq_along(ChIP_LRZ_NRZgeneMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_LRZ_NRZgeneMats[[x]]+1)/(control_LRZ_NRZgeneMats[[x]]+1))
}, mc.cores = length(ChIP_LRZ_NRZgeneMats))

# LRZ_NRZranLoc
log2ChIP_LRZ_NRZranLocMats <- mclapply(seq_along(ChIP_LRZ_NRZranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_LRZ_NRZranLocMats[[x]]+1)/(control_LRZ_NRZranLocMats[[x]]+1))
}, mc.cores = length(ChIP_LRZ_NRZranLocMats))

# armgene
log2ChIP_armgeneMats <- mclapply(seq_along(ChIP_armgeneMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_armgeneMats[[x]]+1)/(control_armgeneMats[[x]]+1))
}, mc.cores = length(ChIP_armgeneMats))

# armranLoc
log2ChIP_armranLocMats <- mclapply(seq_along(ChIP_armranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_armranLocMats[[x]]+1)/(control_armranLocMats[[x]]+1))
}, mc.cores = length(ChIP_armranLocMats))


# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_LRZgeneMats)) {
  colnames(log2ChIP_LRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_LRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_LRZ_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_LRZ_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                 paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                 paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_armgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_armranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIP_LRZgeneMats), function(x) {
  list(
       # LRZgenes
       log2ChIP_LRZgeneMats[[x]],
       # LRZranLocs
       log2ChIP_LRZranLocMats[[x]],
       # NRZgenes
       log2ChIP_NRZgeneMats[[x]],
       # NRZranLocs
       log2ChIP_NRZranLocMats[[x]],
       # LRZ_NRZgenes
       log2ChIP_LRZ_NRZgeneMats[[x]],
       # LRZ_NRZranLocs
       log2ChIP_LRZ_NRZranLocMats[[x]],
       # armgenes
       log2ChIP_armgeneMats[[x]],
       # armranLocs
       log2ChIP_armranLocMats[[x]]
      )
}, mc.cores = length(log2ChIP_LRZgeneMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting
LRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
LRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
LRZ_NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[5]]
})
LRZ_NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[6]]
})
armgenesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[7]]
})
armranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[8]]
})
names(LRZgenesTmp) <- log2ChIPNamesPlot
names(LRZranLocTmp) <- log2ChIPNamesPlot
names(NRZgenesTmp) <- log2ChIPNamesPlot
names(NRZranLocTmp) <- log2ChIPNamesPlot
names(LRZ_NRZgenesTmp) <- log2ChIPNamesPlot
names(LRZ_NRZranLocTmp) <- log2ChIPNamesPlot
names(armgenesTmp) <- log2ChIPNamesPlot
names(armranLocTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(LRZgenesTmp, .id = "libName"),
  bind_rows(LRZranLocTmp, .id = "libName"),
  bind_rows(NRZgenesTmp, .id = "libName"),
  bind_rows(NRZranLocTmp, .id = "libName"),
  bind_rows(LRZ_NRZgenesTmp, .id = "libName"),
  bind_rows(LRZ_NRZranLocTmp, .id = "libName"),
  bind_rows(armgenesTmp, .id = "libName"),
  bind_rows(armranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
                       summaryDFfeature_log2ChIP[[4]]$CI_lower,
                       summaryDFfeature_log2ChIP[[5]]$CI_lower,
                       summaryDFfeature_log2ChIP[[6]]$CI_lower,
                       summaryDFfeature_log2ChIP[[7]]$CI_lower,
                       summaryDFfeature_log2ChIP[[8]]$CI_lower),
                     na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
                       summaryDFfeature_log2ChIP[[4]]$CI_upper,
                       summaryDFfeature_log2ChIP[[5]]$CI_upper,
                       summaryDFfeature_log2ChIP[[6]]$CI_upper,
                       summaryDFfeature_log2ChIP[[7]]$CI_upper,
                       summaryDFfeature_log2ChIP[[8]]$CI_upper),
                     na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 22)))
})


# Plot average profiles with 95% CI ribbon

## LRZgenes
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
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

## NRZgenes
summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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

## LRZ_NRZgenes
summaryDFfeature <- summaryDFfeature_log2ChIP[[5]]
ggObj5_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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
ggtitle(bquote(.(LRZ_NRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## LRZ_NRZranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[6]]
ggObj6_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
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
ggtitle(bquote(.(LRZ_NRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armgenes
summaryDFfeature <- summaryDFfeature_log2ChIP[[7]]
ggObj7_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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
ggtitle(bquote(.(armgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[8]]
ggObj8_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
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
ggtitle(bquote(.(armranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP,
                                              ggObj5_combined_log2ChIP,
                                              ggObj6_combined_log2ChIP,
                                              ggObj7_combined_log2ChIP,
                                              ggObj8_combined_log2ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_",
              paste0(log2ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_LRZgenes_NRZgenes_armgenes_in_", refbase, "_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*8, limitsize = FALSE)


# ChIP
# Add column names
for(x in seq_along(ChIP_LRZgeneMats)) {
  colnames(ChIP_LRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_LRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_LRZ_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_LRZ_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                 paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                 paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_armgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(ChIP_armranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_LRZgeneMats), function(x) {
  list(
       # LRZgenes
       ChIP_LRZgeneMats[[x]],
       # LRZranLocs
       ChIP_LRZranLocMats[[x]],
       # NRZgenes
       ChIP_NRZgeneMats[[x]],
       # NRZranLocs
       ChIP_NRZranLocMats[[x]],
       # LRZ_NRZgenes
       ChIP_LRZ_NRZgeneMats[[x]],
       # LRZ_NRZranLocs
       ChIP_LRZ_NRZranLocMats[[x]],
       # armgenes
       ChIP_armgeneMats[[x]],
       # armranLocs
       ChIP_armranLocMats[[x]]
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
NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
LRZ_NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[5]]
})
LRZ_NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[6]]
})
armgenesTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[7]]
})
armranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[8]]
})
names(LRZgenesTmp) <- ChIPNamesPlot
names(LRZranLocTmp) <- ChIPNamesPlot
names(NRZgenesTmp) <- ChIPNamesPlot
names(NRZranLocTmp) <- ChIPNamesPlot
names(LRZ_NRZgenesTmp) <- ChIPNamesPlot
names(LRZ_NRZranLocTmp) <- ChIPNamesPlot
names(armgenesTmp) <- ChIPNamesPlot
names(armranLocTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(LRZgenesTmp, .id = "libName"),
  bind_rows(LRZranLocTmp, .id = "libName"),
  bind_rows(NRZgenesTmp, .id = "libName"),
  bind_rows(NRZranLocTmp, .id = "libName"),
  bind_rows(LRZ_NRZgenesTmp, .id = "libName"),
  bind_rows(LRZ_NRZranLocTmp, .id = "libName"),
  bind_rows(armgenesTmp, .id = "libName"),
  bind_rows(armranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                                   levels = ChIPNamesPlot)
}

# Define y-axis limits
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
                       summaryDFfeature_ChIP[[2]]$CI_lower,
                       summaryDFfeature_ChIP[[3]]$CI_lower,
                       summaryDFfeature_ChIP[[4]]$CI_lower,
                       summaryDFfeature_ChIP[[5]]$CI_lower,
                       summaryDFfeature_ChIP[[6]]$CI_lower,
                       summaryDFfeature_ChIP[[7]]$CI_lower,
                       summaryDFfeature_ChIP[[8]]$CI_lower),
                     na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
                       summaryDFfeature_ChIP[[2]]$CI_upper,
                       summaryDFfeature_ChIP[[3]]$CI_upper,
                       summaryDFfeature_ChIP[[4]]$CI_upper,
                       summaryDFfeature_ChIP[[5]]$CI_upper,
                       summaryDFfeature_ChIP[[6]]$CI_upper,
                       summaryDFfeature_ChIP[[7]]$CI_upper,
                       summaryDFfeature_ChIP[[8]]$CI_upper),
                     na.rm = T)

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
                   labels = function(x) sprintf("%4.2f", x)) +
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
     y = bquote("ChIP")) +
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
                   labels = function(x) sprintf("%4.2f", x)) +
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
     y = bquote("ChIP")) +
annotation_custom(legendLabs[[1]]) +
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

## NRZgenes
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
                   labels = function(x) sprintf("%4.2f", x)) +
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
     y = bquote("ChIP")) +
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
                   labels = function(x) sprintf("%4.2f", x)) +
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
     y = bquote("ChIP")) +
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

## LRZ_NRZgenes
summaryDFfeature <- summaryDFfeature_ChIP[[5]]
ggObj5_combined_ChIP <- ggplot(data = summaryDFfeature,
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
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
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
ggtitle(bquote(.(LRZ_NRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## LRZ_NRZranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[6]]
ggObj6_combined_ChIP <- ggplot(data = summaryDFfeature,
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
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
annotation_custom(legendLabs[[1]]) +
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
ggtitle(bquote(.(LRZ_NRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armgenes
summaryDFfeature <- summaryDFfeature_ChIP[[7]]
ggObj7_combined_ChIP <- ggplot(data = summaryDFfeature,
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
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
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
ggtitle(bquote(.(armgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[8]]
ggObj8_combined_ChIP <- ggplot(data = summaryDFfeature,
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
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
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
ggtitle(bquote(.(armranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP,
                                              ggObj5_combined_ChIP,
                                              ggObj6_combined_ChIP,
                                              ggObj7_combined_ChIP,
                                              ggObj8_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_",
              paste0(ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_LRZgenes_NRZgenes_armgenes_in_", refbase, "_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*8, limitsize = FALSE)


# control
# Add column names
for(x in seq_along(control_LRZgeneMats)) {
  colnames(control_LRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_LRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_LRZ_NRZgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_LRZ_NRZranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                 paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                 paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_armgeneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(control_armranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                             paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
control_mats <- mclapply(seq_along(control_LRZgeneMats), function(x) {
  list(
       # LRZgenes
       control_LRZgeneMats[[x]],
       # LRZranLocs
       control_LRZranLocMats[[x]],
       # NRZgenes
       control_NRZgeneMats[[x]],
       # NRZranLocs
       control_NRZranLocMats[[x]],
       # LRZ_NRZgenes
       control_LRZ_NRZgeneMats[[x]],
       # LRZ_NRZranLocs
       control_LRZ_NRZranLocMats[[x]],
       # armgenes
       control_armgeneMats[[x]],
       # armranLocs
       control_armranLocMats[[x]]
      )
}, mc.cores = length(control_LRZgeneMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = colnames(control_mats[[x]][[y]]),
               t(control_mats[[x]][[y]]))
  })
}, mc.cores = length(control_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_control[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    tidyDFfeature_list_control[[x]][[y]]$window <- factor(tidyDFfeature_list_control[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    summaryDFfeature_list_control[[x]][[y]]$window <- factor(summaryDFfeature_list_control[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
    summaryDFfeature_list_control[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]])[1])
    summaryDFfeature_list_control[[x]][[y]]$sem <- summaryDFfeature_list_control[[x]][[y]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]]$n-1)
    summaryDFfeature_list_control[[x]][[y]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
    summaryDFfeature_list_control[[x]][[y]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_control into
# a list of single data.frames containing all meta-profiles for plotting
LRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[1]]
})
LRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[2]]
})
NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[3]]
})
NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[4]]
})
LRZ_NRZgenesTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[5]]
})
LRZ_NRZranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[6]]
})
armgenesTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[7]]
})
armranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[8]]
})
names(LRZgenesTmp) <- controlNamesPlot
names(LRZranLocTmp) <- controlNamesPlot
names(NRZgenesTmp) <- controlNamesPlot
names(NRZranLocTmp) <- controlNamesPlot
names(LRZ_NRZgenesTmp) <- controlNamesPlot
names(LRZ_NRZranLocTmp) <- controlNamesPlot
names(armgenesTmp) <- controlNamesPlot
names(armranLocTmp) <- controlNamesPlot
summaryDFfeature_control <- list(
  bind_rows(LRZgenesTmp, .id = "libName"),
  bind_rows(LRZranLocTmp, .id = "libName"),
  bind_rows(NRZgenesTmp, .id = "libName"),
  bind_rows(NRZranLocTmp, .id = "libName"),
  bind_rows(LRZ_NRZgenesTmp, .id = "libName"),
  bind_rows(LRZ_NRZranLocTmp, .id = "libName"),
  bind_rows(armgenesTmp, .id = "libName"),
  bind_rows(armranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_control)) {
  summaryDFfeature_control[[x]]$libName <- factor(summaryDFfeature_control[[x]]$libName,
                                                   levels = controlNamesPlot)
}

# Define y-axis limits
ymin_control <- min(c(summaryDFfeature_control[[1]]$CI_lower,
                       summaryDFfeature_control[[2]]$CI_lower,
                       summaryDFfeature_control[[3]]$CI_lower,
                       summaryDFfeature_control[[4]]$CI_lower,
                       summaryDFfeature_control[[5]]$CI_lower,
                       summaryDFfeature_control[[6]]$CI_lower,
                       summaryDFfeature_control[[7]]$CI_lower,
                       summaryDFfeature_control[[8]]$CI_lower),
                     na.rm = T)
ymax_control <- max(c(summaryDFfeature_control[[1]]$CI_upper,
                       summaryDFfeature_control[[2]]$CI_upper,
                       summaryDFfeature_control[[3]]$CI_upper,
                       summaryDFfeature_control[[4]]$CI_upper,
                       summaryDFfeature_control[[5]]$CI_upper,
                       summaryDFfeature_control[[6]]$CI_upper,
                       summaryDFfeature_control[[7]]$CI_upper,
                       summaryDFfeature_control[[8]]$CI_upper),
                     na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(controlNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(controlNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = controlColours[x], fontsize = 22)))
})


# Plot average profiles with 95% CI ribbon

## LRZgenes
summaryDFfeature <- summaryDFfeature_control[[1]]
ggObj1_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[1]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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
summaryDFfeature <- summaryDFfeature_control[[2]]
ggObj2_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[2]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
annotation_custom(legendLabs[[1]]) +
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

## NRZgenes
summaryDFfeature <- summaryDFfeature_control[[3]]
ggObj3_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[3]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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
summaryDFfeature <- summaryDFfeature_control[[4]]
ggObj4_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[4]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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

## LRZ_NRZgenes
summaryDFfeature <- summaryDFfeature_control[[5]]
ggObj5_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[5]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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
ggtitle(bquote(.(LRZ_NRZgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## LRZ_NRZranLoc
summaryDFfeature <- summaryDFfeature_control[[6]]
ggObj6_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[6]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
annotation_custom(legendLabs[[1]]) +
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
ggtitle(bquote(.(LRZ_NRZranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armgenes
summaryDFfeature <- summaryDFfeature_control[[7]]
ggObj7_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[7]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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
ggtitle(bquote(.(armgenesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## armranLoc
summaryDFfeature <- summaryDFfeature_control[[8]]
ggObj8_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[8]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("control")) +
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
ggtitle(bquote(.(armranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_control,
                                              ggObj2_combined_control,
                                              ggObj3_combined_control,
                                              ggObj4_combined_control,
                                              ggObj5_combined_control,
                                              ggObj6_combined_control,
                                              ggObj7_combined_control,
                                              ggObj8_combined_control
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "control_",
              paste0(controlNames, collapse = "_"),
              "_avgProfiles_around",
              "_LRZgenes_NRZgenes_armgenes_in_", refbase, "_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*8, limitsize = FALSE)
