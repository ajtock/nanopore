#!/applications/R/R-3.5.0/bin/Rscript

# Generate plots of chromosome profiles

# Rscript ./log2_ChIPinput_perwin_PLOTONLY_v080121.R 10kb T2T_Col both 'WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input,WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me3_ChIP14_WT_REC8_Myc_Rep1_input,WT_H3K4me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K27me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,H2AW6_ChIP_SRR5298545_H2AW_input_SRR5298544,H2AW7_ChIP_SRR5298546_H2AW_input_SRR5298544,WT_MNase_Rep1_WT_gDNA_Rep1,WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep2_ChIP_WT_REC8_Myc_Rep1_input,kss_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,cmt3_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,ColLerF1_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep1_WT_gDNA_Rep1_R1' 'CENH3,H3K9me2,H3K4me3,H3K4me1,H3K4me2,H3K27me1,H2A.W.6,H2A.W.7,MNase,REC8-HA,ASY1,MTOPVIB-HA Rep1,MTOPVIB-HA Rep2,wt DMC1 Rep1,wt DMC1 Rep2,kss DMC1 Rep1,cmt3 DMC1 Rep1,ColLerF1 DMC1 Rep1,wt SPO11-1-oligos,met1 SPO11-1-oligos,kss SPO11-1-oligos' 'red,dodgerblue2'

#genomeBinName <- "10kb"
#refbase <- "T2T_Col"
#align <- "both"
#geno1Names <- unlist(strsplit("WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input,WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me3_ChIP14_WT_REC8_Myc_Rep1_input,WT_H3K4me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K27me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,H2AW6_ChIP_SRR5298545_H2AW_input_SRR5298544,H2AW7_ChIP_SRR5298546_H2AW_input_SRR5298544,WT_MNase_Rep1_WT_gDNA_Rep1,WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep2_ChIP_WT_REC8_Myc_Rep1_input,kss_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,cmt3_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,ColLerF1_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep1_WT_gDNA_Rep1_R1",
#                              split = ","))
#geno1NamesPlot <- unlist(strsplit("CENH3,H3K9me2,H3K4me3,H3K4me1,H3K4me2,H3K27me1,H2A.W.6,H2A.W.7,MNase,REC8-HA,ASY1,MTOPVIB-HA Rep1,MTOPVIB-HA Rep2,wt DMC1 Rep1,wt DMC1 Rep2,kss DMC1 Rep1,cmt3 DMC1 Rep1,ColLerF1 DMC1 Rep1,wt SPO11-1-oligos,met1 SPO11-1-oligos,kss SPO11-1-oligos",
#                                  split = ","))
#geno1Colours <- unlist(strsplit("red,dodgerblue2",
#                                split = ","))

args <- commandArgs(trailingOnly = T)
genomeBinName <- as.character(args[1])
refbase <- args[2]
align <- args[3]
geno1Names <- unlist(strsplit(args[4],
                              split = ","))
#geno2Names <- unlist(strsplit(args[5],
#                              split = ","))
geno1NamesPlot <- unlist(strsplit(args[5],
                                  split = ","))
#geno2NamesPlot <- unlist(strsplit(args[6],
#                                  split = ","))
geno1Colours <- unlist(strsplit(args[6],
                                split = ","))
#geno2Colours <- unlist(strsplit(args[7],
#                                split = ","))

options(stringsAsFactors = FALSE)
library(parallel)

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1][1:5])
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

geno1_filt_winDF_list <- mclapply(seq_along(geno1Names), function(x) {
  read.table(paste0(geno1Names[x],
                    "_MappedOn_", refbase, "_lowXM_", align,
                    "_sort_norm_binSize",
                     genomeBinName, "_smoothed.tsv"),
             header = T)
}, mc.cores = length(geno1Names))

## smallRNAs
#sRNAsizes <- c("all", "21nt", "22nt", "23nt", "24nt")
#sRNAsizesPlot <- c("All", "21-nt", "22-nt", "23-nt", "24-nt")
#sRNAsizesColours <- c("grey40", "darkorange2", "purple3", "green", "blue")
#
#sRNA1_filt_winDF_list <- mclapply(seq_along(sRNAsizes), function(x) {
#  read.table(paste0("filt_WT_sRNAseq_Rep1_SRR6376831",
#                    "_genome_norm_coverage_",
#                    genomeBinName, "_", sRNAsizes[x], ".tsv"),
#             header = T)
#}, mc.cores = length(sRNAsizes))
#sRNA2_filt_winDF_list <- mclapply(seq_along(sRNAsizes), function(x) {
#  read.table(paste0("filt_WT_sRNAseq_Rep2_SRR6376832",
#                    "_genome_norm_coverage_",
#                    genomeBinName, "_", sRNAsizes[x], ".tsv"),
#             header = T)
#}, mc.cores = length(sRNAsizes))

# DNA methylation
filt_DNAmeth <- read.table(paste0("DNAmeth_Col0_BSseq_Rep1_MappedOn_", refbase, "_dedup_binSize200kb_smoothed.tsv"),
                           header = T)

# Function to plot genome-scale coverage overlaid with other datasets
dat1to2diffY <- function(xplot,
                         dat1, dat2,
                         dat1Lab, dat2Lab,
                         dat1Colour, dat2Colour) {
  plot(xplot, dat1, type = "l", lwd = 2, col = dat1Colour,
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        text = dat1Lab, col = dat1Colour)
  par(new = T)
  plot(xplot, dat2, type = "l", lwd = 2, col = dat2Colour,
       ylim = c(min(dat2),
                max(dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.25), labels = dat2Lab, xpd = NA, srt = -90,
       col = dat2Colour)
  axis(side = 4, at = pretty(dat2), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8, 1.4e8),
       labels = c("0", "20", "40", "60", "80", "100", "120", "140"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  #abline(v = centromeres, lty = 5, lwd = 2)
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
}

# Function to plot genome-scale coverage overlaid with DNA methylation (each context)
ChIPvsDNAmethSepGenomePlot <- function(xplot1,
                                       xplot2,
                                       dat1A,
                                       dat2,
                                       dat1ALab,
                                       legendNames, plotColours) {
  plot(xplot1, dat1A, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1A), max(dat1A)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        #text = expression("Log"[2]*"(ChIP/input)"))
        text = dat1ALab, col = plotColours[1])
  par(new = T)
  plot(xplot2, dat2[,4], type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat2[,4], dat2[,5], dat2[,6]), max(dat2[,4], dat2[,5], dat2[,6])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.95), labels = "DNA methylation", xpd = NA, srt = -90,
       col = "blue")
  lines(xplot2, dat2[,5], type = "l", lwd = 2, col = plotColours[3])
  lines(xplot2, dat2[,6], type = "l", lwd = 2, col = plotColours[4])
  axis(side = 4, at = pretty(c(dat2[,4], dat2[,5], dat2[,6])), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8, 1.4e8),
       labels = c("0", "20", "40", "60", "80", "100", "120", "140"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  #abline(v = centromeres, lty = 5, lwd = 2)
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours[2:4],
         ncol = 1, cex = 1.4, lwd = 1.5, bty = "n")
}

# Function to plot genome-scale coverage overlaid with other datasets 
dat1GenomePlot <- function(xplot,
                           dat1,
                           Ylab,
                           plotColours) {
  plot(xplot, dat1, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8, 1.4e8),
       labels = c("0", "20", "40", "60", "80", "100", "120", "140"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2.0,
        text = Ylab, col = plotColours[1])
  abline(v = sumchr, lty = 1, lwd = 2)
  #abline(v = centromeres, lty = 5, lwd = 2)
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
}

# Function to plot genome-scale coverage overlaid with other datasets 
dat1to2GenomePlot <- function(xplot,
                              dat1, dat2,
                              Ylab,
                              legendNames,
                              plotColours) {
  plot(xplot, dat2, type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8, 1.4e8),
       labels = c("0", "20", "40", "60", "80", "100", "120", "140"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2.0,
        text = Ylab)
  abline(v = sumchr, lty = 1, lwd = 2)
  #abline(v = centromeres, lty = 5, lwd = 2)
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  #rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.4, lwd = 1.5, bty = "n")
}

# Plot
pdf(paste0(plotDir,
           "log2ChIPcontrol_genomeProfiles_", genomeBinName, "_" , refbase, "_v080121.pdf"),
    height = 70, width = 32)
par(mfcol = c(14, 2))
par(mar = c(5.1, 7.1, 2.1, 7.1))
par(mgp = c(3, 1.5, 0))
# DNA methylation vs CENH3 (log2)
ChIPvsDNAmethSepGenomePlot(xplot1 = geno1_filt_winDF_list[[1]]$cumwindow,
                           xplot2 = filt_DNAmeth$cumwindow,
                           dat1A = geno1_filt_winDF_list[[1]][,6],
                           dat2 = filt_DNAmeth,
                           dat1ALab = bquote("Log"[2]*"("*.(geno1NamesPlot[1])*"/control)"),
                           legendNames = c("mCG", "mCHG", "mCHH"),
                           plotColours = c(geno1Colours[1], "navy", "blue", "deepskyblue1"))
# histone modifications and variants vs CENH3
for(x in 2:(length(geno1Names)-3)) {
dat1to2diffY(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
             dat1 = geno1_filt_winDF_list[[x]][,6],
             dat2 = geno1_filt_winDF_list[[1]][,6],
             dat1Lab = bquote("Log"[2]*"("*.(geno1NamesPlot[x])*"/control)"), 
             dat2Lab = bquote("Log"[2]*"("*.(geno1NamesPlot[1])*"/control)"), 
             dat1Colour = geno1Colours[2],
             dat2Colour = geno1Colours[1])
}
for(x in (length(geno1Names)-2):length(geno1Names)) {
dat1to2diffY(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
             dat1 = geno1_filt_winDF_list[[x]][,6],
             dat2 = geno1_filt_winDF_list[[1]][,6],
             dat1Lab = bquote("Log"[2]*"("*.(geno1NamesPlot[x])*"/gDNA)"), 
             dat2Lab = bquote("Log"[2]*"("*.(geno1NamesPlot[1])*"/control)"), 
             dat1Colour = geno1Colours[2],
             dat2Colour = geno1Colours[1])
}

dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-2]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)-1]][,6],
                  Ylab = bquote("Log"[2]*"(SPO11-1-oligos/gDNA)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-2], geno1NamesPlot[length(geno1Names)-1]),
                  plotColours = c("navy", "green2"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-2]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)]][,6],
                  Ylab = bquote("Log"[2]*"(SPO11-1-oligos/gDNA)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-2], geno1NamesPlot[length(geno1Names)]),
                  plotColours = c("navy", "magenta3"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-1]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)]][,6],
                  Ylab = bquote("Log"[2]*"(SPO11-1-oligos/gDNA)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-1], geno1NamesPlot[length(geno1Names)]),
                  plotColours = c("green2", "magenta3"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-6]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-6]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)-5]][,6],
                  Ylab = bquote("Log"[2]*"(ChIP/control)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-6], geno1NamesPlot[length(geno1Names)-5]),
                  plotColours = c("navy", "magenta3"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-6]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-6]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)-4]][,6],
                  Ylab = bquote("Log"[2]*"(ChIP/control)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-6], geno1NamesPlot[length(geno1Names)-4]),
                  plotColours = c("navy", "darkgreen"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-5]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-5]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)-4]][,6],
                  Ylab = bquote("Log"[2]*"(ChIP/control)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-5], geno1NamesPlot[length(geno1Names)-4]),
                  plotColours = c("magenta3", "darkgreen"))
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[length(geno1Names)-6]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[length(geno1Names)-6]][,6],
                  dat2 = geno1_filt_winDF_list[[length(geno1Names)-3]][,6],
                  Ylab = bquote("Log"[2]*"(ChIP/control)"),
                  legendNames = c(geno1NamesPlot[length(geno1Names)-6], geno1NamesPlot[length(geno1Names)-3]),
                  plotColours = c("navy", "darkorange"))
dev.off()
