#!/applications/R/R-4.0.0/bin/Rscript

# Usage:
# /applications/R/R-4.0.0/bin/Rscript per_read_methyl_density_heatmap.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 1000 200000 CpG

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#readBinSize <- 1000
#genomeBinSize <- 200000
#context <- "CpG"

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
readBinSize <- as.integer(args[3])
genomeBinSize <- as.integer(args[4])
context <- args[5]

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

if(floor(log10(readBinSize)) + 1 < 4) {
  readBinName <- paste0(readBinSize, "bp")
} else if(floor(log10(readBinSize)) + 1 >= 4 &
          floor(log10(readBinSize)) + 1 <= 6) {
  readBinName <- paste0(readBinSize/1e3, "kb")
} else if(floor(log10(readBinSize)) + 1 >= 7) {
  readBinName <- paste0(readBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(segmentSeq)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[1:5]
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

CEN <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.centromeres"), header = T)
CENstart <- CEN$start
CENend <- CEN$end
CENGR <- GRanges(seqnames = CEN$chr,
                 ranges = IRanges(start = CEN$start,
                                  end = CEN$end),
                 strand = "*")

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# Load table of per-read methylation proportions, in which
# each coordinate corresponds to the midpoint between the
# first and last cytosine position with methylation info in the read
tab <- read.table(paste0(sampleName, "_MappedOn_", refbase, "_", context,
                         "_raw_winSize", readBinName, "_per_read_midpoint.tsv"),
                  header = T)
colnames(tab) <- c("chr", "midpoint", "per_read_mProp")
tab[,2] <- round(tab[,2])
tab <- tab[with(tab, order(chr, midpoint, decreasing = FALSE)),]

# Define heatmap colours
rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
rich8 <- function() {manual_pal(values = c("#000041","#0000CB","#0081FF","#02DA81","#80FE1A","#FDEE02","#FFAB00","#FF3300"))}
rich6 <- function() {manual_pal(values = c("#000043","#0033FF","#01CCA4","#BAFF12","#FFCC00","#FF3300"))}
revSpectralScale11 <- rev(brewer.pal(11, "Spectral"))
viridisScale <- viridis_pal()(6)
#htmpColour <- rich6()(6)
htmpColour <- revSpectralScale11

# For each genomeBinSize-bp adjacent window, calculate the mean of the genomeBinSize
# per-base DNA methylation proportions within that window
print(genomeBinName)
htmps <- NULL
for(i in seq_along(chrs)) {
  # Define adjacent windows
  winSeq <- seq(from = 1, to = chrLens[i], by = genomeBinSize)
  winCum <- winSeq + sumchr[i]
  winIR <- IRanges(start = winSeq,
                   width = genomeBinSize)
  winIR <- winIR[-length(winIR)]
  winIR <- append(winIR,
                  IRanges(start = winSeq[length(winSeq)],
                          end = chrLens[i]))
  winGR <- GRanges(seqnames = chrs[i],
                   ranges = winIR,
                   strand = "*")
  print(winGR)

  # Define per-read midpoint coordinates as GRanges objects
  # and get corresponding DNA methylation proportions that overlap genomic windows
  chr_tab <- tab[tab[,1] == chrs[i],]
  chr_tab_GR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chr_tab[,2],
                                         width = 1),
                        strand = "*",
                        per_read_mProp = chr_tab[,3])

#  ## Note: findOverlaps() approach does not work where a window does not overlap
#  ##       andy positions in chr_tab_GR, which can occur with smaller genomeBinSize
#  # Identify overlapping windows and midpoint coordinates
#  fOverlaps <- findOverlaps(query = winGR,
#                            subject = chr_tab_GR,
#                            type = "any",
#                            select = "all",
#                            ignore.strand = TRUE)
#  # Convert fOverlaps into list object equivalent to that
#  # generated by segmentSeq::getOverlaps(), in which each
#  # list element corresponds to a vector of sequentially numbered indices of
#  # read midpoint coordinates that overlap a given genomic window
#  fOverlapsList <- mclapply(seq_along(unique(queryHits(fOverlaps))),
#                            function(x) {
#                              subjectHits(fOverlaps[queryHits(fOverlaps) == x])
#                            }, mc.cores = detectCores(), mc.preschedule = TRUE)
  fOverlapsList <- getOverlaps(coordinates = winGR,
                               segments = chr_tab_GR,
                               overlapType = "overlapping",
                               whichOverlaps = TRUE,
                               ignoreStrand = TRUE)

  # Get per-read methylation proportion values overlapping each window
  win_mProp_list <- lapply(fOverlapsList, function(x) {
                      data.frame(matrix(data = chr_tab[,3][x], nrow = 1))
                    })
  # Convert into matrix in which each column corresponds to a genomic window
  win_mProp_matrix <- t(as.matrix(x = bind_rows(win_mProp_list)))
  colnames(win_mProp_matrix) <- round(start(winGR)/1e6, digits = 1)
  # Remove columns where less than 2 rows are not NA
  win_mProp_matrix <- win_mProp_matrix[,colSums(is.na(win_mProp_matrix)) < nrow(win_mProp_matrix) - 1]  

  # Define ylim depending on context
  if(context == "CpG") {
    ylimContext <- c(0, 1)
  } else if(context == "CHG") {
    ylimContext <- c(0, 0.6)
  } else if(context == "CHH") {
    ylimContext <- c(0, 0.2)
  }

  # Make heatmap for chromosome
  ha1 <- HeatmapAnnotation(Region = c(
                           rep("NonCEN", length(which(as.numeric(colnames(win_mProp_matrix))*1e6 < CENstart[i] &
                                                       as.numeric(colnames(win_mProp_matrix))*1e6 < CENend[i]))),
                           rep("CEN", length(which(as.numeric(colnames(win_mProp_matrix))*1e6 > CENstart[i] &
                                                   as.numeric(colnames(win_mProp_matrix))*1e6 < CENend[i]))),
                           rep("NonCEN", length(which(as.numeric(colnames(win_mProp_matrix))*1e6 > CENend[i])))),
                           col = list(Region = c("NonCEN" = "grey70", "CEN" = "red")),
                           annotation_legend_param = list(title = "Region",
                                                          title_position = "topleft",
                                                          title_gp = gpar(font = 2, fontsize = 12),
                                                          labels_gp = gpar(font = 3, fontize = 10),
                                                          legend_direction = "vertical",
                                                          nrow = 2,
                                                          labels_gp = gpar(fontsize = 10)),
                           show_annotation_name = FALSE)

  htmp <- densityHeatmap(data = win_mProp_matrix,
                         col = htmpColour,
                         density_param = list(na.rm = TRUE),
                         top_annotation = ha1,
                         column_title = paste0(chrs[i], " ",
                                               gsub(pattern = "[A-z]", replacement = "", x = genomeBinName),
                                               "-", gsub(pattern = "[0-9]", replacement = "", x = genomeBinName),
                                               " window (Mb)"),
                         column_title_side = "bottom",
                         column_title_rot = 0,
                         title_gp = gpar(fontsize = 12),
                         show_quantiles = FALSE,
                         ylab = paste0("Per-read m", context, " proportion"),
                         ylab_gp = gpar(fontsize = 12),
                         ylim = ylimContext, 
                         column_names_side = "bottom",
                         column_names_gp = gpar(fontsize = 6),
                         column_names_rot = 90,
                         column_names_centered = TRUE,
                         column_gap = unit(0, "mm"),
                         heatmap_legend_param = list(title = "Density",
                                                     title_position = "topleft",
                                                     title_gp = gpar(font = 2, fontsize = 12),
                                                     legend_direction = "vertical",
                                                     labels_gp = gpar(fontsize = 10)),
                         border = FALSE,
                         # If converting into png with pdfTotiffTopng.sh,
                         # set use_raster to FALSE
                         use_raster = FALSE)
                         #use_raster = TRUE, raster_device = "png", raster_quality = 4)
  htmps <- htmps + htmp
}

pdf(paste0(plotDir,
           sampleName, "_MappedOn_", refbase, "_", context,
           "_prop_raw_readBinSize", readBinName, "_per_read_midpoint_density_heatmap_genomeBinSize", genomeBinName, ".pdf"), 
    height = 4, width = 12 * length(htmps))
draw(htmps,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     gap = unit(c(1), "mm"))
dev.off()
