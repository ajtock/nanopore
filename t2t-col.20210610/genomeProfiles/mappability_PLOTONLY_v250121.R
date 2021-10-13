#!/applications/R/R-4.0.0/bin/Rscript

# Generate plots of chromosome profiles

# Rscript ./mappability_PLOTONLY_v250121.R 10kb T2T_Col both 'WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input' 'CENH3' 'red'

#genomeBinName <- "10kb"
#refbase <- "T2T_Col"
#align <- "both"
#geno1Names <- unlist(strsplit("WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input",
#                              split = ","))
#geno1NamesPlot <- unlist(strsplit("CENH3",
#                                  split = ","))
#geno1Colours <- unlist(strsplit("red",
#                                split = ","))

args <- commandArgs(trailingOnly = T)
genomeBinName <- as.character(args[1])
refbase <- args[2]
align <- args[3]
geno1Names <- unlist(strsplit(args[4],
                              split = ","))
geno1NamesPlot <- unlist(strsplit(args[5],
                                  split = ","))
geno1Colours <- unlist(strsplit(args[6],
                                split = ","))

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
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

CENstart <- sapply(seq_along(CENstart), function(x) {
  CENstart[x] + sumchr[x]
})
print(CENstart)
CENend <- sapply(seq_along(CENend), function(x) {
  CENend[x] + sumchr[x]
})
print(CENend)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

geno1_filt_winDF_list <- mclapply(seq_along(geno1Names), function(x) {
  read.table(paste0(geno1Names[x],
                    "_MappedOn_", refbase, "_lowXM_", align,
                    "_sort_norm_binSize",
                     genomeBinName, "_smoothed.tsv"),
             header = T)
}, mc.cores = length(geno1Names))

mapNames <- c("map_K40_E2", "map_K45_E2", "map_K50_E2",
              "map_K150_E4", "map_K200_E4", "map_K300_E4")
mapNamesPlot <- c("k=40 e=2", "k=45 e=2", "k=50 e=2",
                  "k=150 e=4", "k=200 e=4", "k=300 e=4")
mapWins <- read.table(paste0(mapNames[1],
                             "_genmap_MappedOn_", refbase,
                             "_binSize",
                              genomeBinName, "_smoothed.tsv"),
                      header = T, colClasses = c(rep(NA, 3), "NULL"))
filt_MAP_list <- mclapply(seq_along(mapNames), function(x) {
  read.table(paste0(mapNames[x],
                    "_genmap_MappedOn_", refbase,
                    "_binSize",
                     genomeBinName, "_smoothed.tsv"),
             header = T, colClasses = c(rep("NULL", 3), NA))
}, mc.cores = length(mapNames))
filt_map <- data.frame(mapWins,
                       do.call(cbind, filt_MAP_list))
colnames(filt_map) <- c(colnames(mapWins), mapNames)

# Function to plot genome-scale coverage overlaid with mappability (each in mapNames)
ChIPvsMapGenomePlot <- function(xplot1,
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
       ylim = c(min(dat2[,4], dat2[,5], dat2[,6],
                    dat2[,7], dat2[,8], dat2[,9]),
                max(dat2[,4], dat2[,5], dat2[,6],
                    dat2[,7], dat2[,8], dat2[,9])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 2.00, adj = c(0.5, -2.95), labels = "Mappability", xpd = NA, srt = -90,
       col = "black")
  lines(xplot2, dat2[,5], type = "l", lwd = 2, col = plotColours[3])
  lines(xplot2, dat2[,6], type = "l", lwd = 2, col = plotColours[4])
  lines(xplot2, dat2[,7], type = "l", lwd = 2, col = plotColours[5])
  lines(xplot2, dat2[,8], type = "l", lwd = 2, col = plotColours[6])
  lines(xplot2, dat2[,9], type = "l", lwd = 2, col = plotColours[7])
  axis(side = 4, at = pretty(c(dat2[,4], dat2[,5], dat2[,6])), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8, 1.4e8),
       labels = c("0", "20", "40", "60", "80", "100", "120", "140"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = c(CENstart, CENend), lty = 5, lwd = 2, col = "grey20")
  #rug(x = c(CENstart, CENend), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  #rug(x = c(CENstart, CENend), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("left",
         legend = legendNames,
         col = "white",
         text.col = plotColours[2:7],
         ncol = 1, cex = 1.0, lwd = 1.5, bty = "n")
}

# Plot
pdf(paste0(plotDir,
           "log2ChIPcontrol_mappability_genomeProfiles_", genomeBinName, "_" , refbase, "_v250121.pdf"),
    height = 5, width = 16)
par(mfcol = c(1, 1))
par(mar = c(5.1, 7.1, 2.1, 7.1))
par(mgp = c(3, 1.5, 0))
# Mappability vs CENH3 (log2[ChIP/input])
ChIPvsMapGenomePlot(xplot1 = geno1_filt_winDF_list[[1]]$cumwindow,
                    xplot2 = filt_map$cumwindow,
                    dat1A = geno1_filt_winDF_list[[1]][,6],
                    dat2 = filt_map,
                    dat1ALab = bquote("Log"[2]*"("*.(geno1NamesPlot[1])*"/control)"),
                    legendNames = mapNamesPlot,
                    plotColours = c(geno1Colours[1],
                                    "navy", "blue", "deepskyblue1",
                                    "darkgreen", "seagreen", "springgreen3"))
dev.off()
