#!/applications/R/R-3.5.0/bin/Rscript

#
# author: Andy Tock (adapted from original script logistic2.R by Tom Hardcastle)
# contact: ajt200@cam.ac.uk
# date: 07.10.2020
# 

# Build GLM with GBS-derived Col/Ler crossovers as the response variable and
# Col/Ler SNPs and interval width as pedictor variables

# Usage:
# Rscript ./logistic2_ajtEdit_v180520.R 1000 1

#winSize <- 1000
#stepSize <- 1

args <- commandArgs(trailingOnly = T)
winSize <- as.numeric(args[1])
stepSize <- as.numeric(args[2])

library(segmentSeq)
library(glm2)
library(MASS)

outDir <- paste0(winSize/1e3, "kbWin_", stepSize, "bpStep/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
genomeFAI <- read.table("/projects/ajt200/TAIR10/TAIR10_chr_all.fa.fai",
                        colClasses = c(rep(NA, 2), rep("NULL", 3)), header = F)[1:5,]
chrs <- paste0("Chr", genomeFAI$V1)
chrLens <- genomeFAI$V2

# Load SNPs (519834 SNPs in collerF2.complete.tiger.txt and 481252 SNPs in BC.complete.tiger.txt)
# and create a GRanges object of 519839 intervals within which a crossover can be detected
# There are 5 more intervals than SNPs due to the interval
# between the last SNP in each chromosome and the end of that chromosome
SNPs1 <- read.table("collerF2.complete.tiger.txt",
                    header = T)
print(dim(SNPs1))
#[1] 519834      6
SNPs1 <- SNPs1[,1:2]
colnames(SNPs1) <- c("chr", "pos")
SNPs2 <- read.table("BC.complete.tiger.txt",
                    header = T)
print(dim(SNPs2))
#[1] 481252	6
SNPs2 <- SNPs2[,1:2]
colnames(SNPs2) <- c("chr", "pos")
SNPs <- rbind(SNPs1[,1:2], SNPs2[,1:2])
print(dim(SNPs))
#[1] 1001086       2
SNPs <- unique(SNPs)
SNPs <- SNPs[with(SNPs, order(chr, pos)),]
print(dim(SNPs))
#[1] 534775
SNPsGR <- do.call("c", lapply(1:5, function(chr) {
  chrSNPs <- SNPs[SNPs[,1] == chr,]
  GRanges(seqnames = chrs[chr],
          ranges = IRanges(start = c(1, chrSNPs[,2]),
                           end = c(chrSNPs[,2], chrLens[chr])))
}))
print(SNPsGR)

# Load 3320 wild type crossovers
COsBED <- read.table("/projects/ajt200/GBS_CO/HS_CU_080617/wt/GBS_crossovers.bed",
                     header = F)
# Note addition of 1 to start coordinates as these are in BED format
# (i.e., 0-based start coordinates)
# NOTE: Tom merged overlapping COs (i.e., reduce(COsGR), which gives 3115 crossovers) 
COsGR <- GRanges(seqnames = COsBED$V1,
                 ranges = IRanges(start = COsBED$V2+1,
                                  end = COsBED$V3))
COsGR <- reduce(COsGR)

# Identify which SNP intervals overlap at least one crossover interval
# Tom set overlapType to "within", meaning that called overlaps correspond to
# SNP intervals that are completely contained by at least one crossover interval
COevents <- getOverlaps(coordinates = SNPsGR,
                        segments = COsGR,
                        overlapType = "within",
                        whichOverlaps = FALSE,
                        ignoreStrand = TRUE)
# Number of SNP intervals completely contained by at least one crossover interval
print(sum(COevents))
#[1] 4865
# Number of SNP intervals not completely contained by at least one crossover interval
print(sum(!COevents))
#[1] 529915 

# Define regionsGR object as SNP intervals that are
# not completely contained by at least one crossover interval
# + all crossover intervals

regionsGR <- c(SNPsGR[!COevents], COsGR)
regionsGR <- sort(regionsGR)
print(length(regionsGR))
#[1] 533030

# Load SPO11-1-oligos, nuclesome occupancy, H3K4me3 ChIP-seq and DNA methylation data
# and calculate sum and mean in each region

SPO11 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
                           "log2wtSPO11oligoMeanAllRepsNakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
colnames(SPO11) <- c("chr", "signal")
#SPO11_Rep1 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_Rep2 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_Rep3 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_df <- data.frame(Rep1 = SPO11_Rep1$V4,
#                       Rep2 = SPO11_Rep2$V4,
#                       Rep3 = SPO11_Rep3$V4,
#                       stringsAsFactors = F)
#SPO11 <- data.frame(chr = SPO11_Rep1$V1,
#                    signal = rowMeans(SPO11_df),
#                    stringsAsFactors = F)
splitSPO11 <- list()
for(x in 1:5) {
  chr_SPO11 <- SPO11[SPO11$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitSPO11 <- c(splitSPO11,
                  split(x = chr_SPO11,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanSPO11 <- sapply(splitSPO11, mean)
sumSPO11 <- sapply(splitSPO11, sum)

MNase <- read.table(paste0("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/",
                           "log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
colnames(MNase) <- c("chr", "signal")
splitMNase <- list()
for(x in 1:5) {
  chr_MNase <- MNase[MNase$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitMNase <- c(splitMNase,
                  split(x = chr_MNase,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanMNase <- sapply(splitMNase, mean)
sumMNase <- sapply(splitMNase, sum)

H3K4me3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/WT/coverage/log2ChIPinput/noZscore/",
                             "log2wtH3K4me3ChIPwtH3K9me2input_noZscore_norm_allchrs_coverage_coord_tab.bed"),
                      colClasses = c(NA, rep("NULL", 2), NA),
                      header = F)
colnames(H3K4me3) <- c("chr", "signal")
#H3K4me3_Rep1 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_Rep2 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_Rep3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_df <- data.frame(
#                         Rep1 = H3K4me3_Rep1$V4,
#                         Rep2 = H3K4me3_Rep2$V4,
#                         Rep3 = H3K4me3_Rep3$V4,
#                         stringsAsFactors = F)
#H3K4me3 <- data.frame(chr = H3K4me3_Rep1$V1,
#                      signal = rowMeans(H3K4me3_df),
#                      stringsAsFactors = F)
splitH3K4me3 <- list()
for(x in 1:5) {
  chr_H3K4me3 <- H3K4me3[H3K4me3$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitH3K4me3 <- c(splitH3K4me3,
                  split(x = chr_H3K4me3,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanH3K4me3 <- sapply(splitH3K4me3, mean)
sumH3K4me3 <- sapply(splitH3K4me3, sum)

# DNA methylation
CG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                        "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"),
                 colClasses = c(rep(NA, 2), "NULL", NA),
                 header = F)
CHG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
CHH <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
meanCG <- NULL
meanCHG <- NULL
meanCHH <- NULL
meanDNAmeth <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  windowsGR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = start(windowsGR),
                                        end = end(windowsGR)-1))
  end(windowsGR[length(windowsGR)]) <- chrLens[x]
  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  # CG
  chr_CG <- CG[CG[,1] == paste0("chr", x),]
  chr_CG_GR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = chr_CG[,2],
                                        width = 1),
                       strand = "*")
  CGoverlaps <- getOverlaps(coordinates = windowsGR,
                            segments = chr_CG_GR,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE,
                            ignoreStrand = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chr_CG[,3][x])))
  # CHG
  chr_CHG <- CHG[CHG[,1] == paste0("chr", x),]
  chr_CHG_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHG[,2],
                                         width = 1),
                        strand = "*")
  CHGoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHG_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHGwinVals <- sapply(CHGoverlaps, function(x) mean(as.numeric(chr_CHG[,3][x])))
  # CHH
  chr_CHH <- CHH[CHH[,1] == paste0("chr", x),]
  chr_CHH_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHH[,2],
                                         width = 1),
                        strand = "*")
  CHHoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHH_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHHwinVals <- sapply(CHHoverlaps, function(x) mean(as.numeric(chr_CHH[,3][x])))
  # Mean of all 3 contexts
  CwinVals <- sapply(seq_along(CGwinVals), function(x) {
                mean(c(CGwinVals[x], CHGwinVals[x], CHHwinVals[x]))
              })
  # Combine values from all chromosomes
  meanCG <- c(meanCG, CGwinVals)
  meanCHG <- c(meanCHG, CHGwinVals)
  meanCHH <- c(meanCHH, CHHwinVals)
  meanDNAmeth <- c(meanDNAmeth, CwinVals)
}

# Calculate SNP frequency within each interval in regionsGR
allSNPsGR <- GRanges(seqnames = paste0("Chr", SNPs[,1]),
                     ranges = IRanges(start = SNPs[,2],
                                      width = 1))
# Count SNPs within each interval in which a crossover could potentially have been detected
# and divide by the width of the corresponding interval
meanSNPs <- NULL
sumSNPs <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  chr_allSNPsGR <- allSNPsGR[seqnames(allSNPsGR) == chrs[x]]
  winSNPs <- countOverlaps(query = windowsGR,
                           subject = chr_allSNPsGR,
                           type = "any",
                           ignore.strand = TRUE)
  meanSNPs <- c(meanSNPs, winSNPs/width(windowsGR))
  sumSNPs <- c(sumSNPs, winSNPs)
} 

# Count SNPs in windows of width winSize nt, with a step of stepSize nt
SNPsPB <- data.frame()
for(x in 1:5) {
  # Define sliding windows of width winSize nt,
  # with a step of stepSize nt
  ## Note: the commented-out code would create windows of winSize nt only,
  ## whereas the active code creates windows decreasing from winSize nt to stepSize nt
  ## at the right-hand end of each chromosome ( from chrLens[x]-winSize to chrLens[x] ),
  winStarts <- seq(from = 1,
                   to = chrLens[x],
#                   to = chrLens[x]-winSize,
                   by = stepSize)
  stopifnot(winStarts[length(winStarts)] == chrLens[x])
#  if(chrLens[x] - winStarts[length(winStarts)] >= winSize) {
#    winStarts <- c(winStarts,
#                   winStarts[length(winStarts)]+stepSize)
#  }
  winEnds <- seq(from = winStarts[1]+winSize-1,
                 to = chrLens[x],
                 by = stepSize)
  stopifnot(winEnds[length(winEnds)] == chrLens[x])
  winEnds <- c(winEnds,
               rep(chrLens[x], times = length(winStarts)-length(winEnds)))
  stopifnot(length(winStarts) == length(winEnds))

  windowsGR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = winStarts,
                                        end = winEnds),
                       strand = "*")

  # Count SNPs in sliding windows
  chr_allSNPsGR <- allSNPsGR[seqnames(allSNPsGR) == chrs[x]]
  winSNPs <- countOverlaps(query = windowsGR,
                           subject = chr_allSNPsGR,
                           type = "any",
                           ignore.strand = TRUE)
  SNPsPBDF <- data.frame(chr = chrs[x],
                         signal = winSNPs/width(windowsGR)) 
  SNPsPB <- rbind(SNPsPB, SNPsPBDF) 
}

splitSNPsPB <- list()
for(x in 1:5) {
  chr_SNPsPB <- SNPsPB[SNPsPB$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitSNPsPB <- c(splitSNPsPB,
                   split(x = chr_SNPsPB,
                         f = c(rep(1:length(chr_regionsGR),
                                   times = width(chr_regionsGR)-1),
                               length(chr_regionsGR))))
}
meanSNPsPB <- sapply(splitSNPsPB, mean)
sumSNPsPB <- sapply(splitSNPsPB, sum)

# Create data object for model
dat <- cbind.data.frame(CO = regionsGR %in% COsGR,
                        meanSPO11 = meanSPO11,
                        meanMNase = meanMNase,
                        meanH3K4me3 = meanH3K4me3,
                        meanDNAmeth = meanDNAmeth,
                        meanSNPs = meanSNPs,
                        meanSNPsPB = meanSNPsPB,
                        sumSPO11 = sumSPO11,
                        sumMNase = sumMNase,
                        sumH3K4me3 = sumH3K4me3,
                        sumSNPs = sumSNPs,
                        sumSNPsPB = sumSNPsPB,
                        width = width(regionsGR))
save(dat, file = paste0(outDir, "df_for_GLM_binomial_logit_",
                        winSize/1e3, "kbWin_", stepSize, "bpStep.RData"))

## Discard rows with missing data
dat <- dat[!is.na(dat$meanDNAmeth),]

# Build binomial GLM with "logit" link function
glmCO <- glm2(formula = CO ~ (meanSPO11 + meanMNase + meanH3K4me3 + meanDNAmeth + meanSNPsPB + width)^2,
              family = binomial(link="logit"),
              control = glm.control(maxit = 100000),
              data = dat)

glm_stepAIC <- stepAIC(object = glmCO, direction = "both")
print("stepAIC-selected model formula:")
print(glm_stepAIC$formula)
glm_select <- glm2(formula = glm_stepAIC$formula, 
                   family = binomial(link="logit"),
                   control = glm.control(maxit = 100000),
                   data = dat)
glm_summary <- summary(glm_select)
glm_coeffs <- glm_summary$coefficients
glm_predict <- predict(glm_select, type = "response")
glm_formula <- glm_select$formula
save(glm_stepAIC, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_stepAIC.RData"))
save(glm_select, file = paste0(outDir, "GLM_binomial_logit_",
                               winSize/1e3, "kbWin_", stepSize, "bpStep.RData"))
save(glm_summary, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_summary.RData"))
write.csv(glm_coeffs, file = paste0(outDir, "GLM_binomial_logit_",
                                    winSize/1e3, "kbWin_", stepSize, "bpStep_coeff.csv"))
save(glm_predict, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_predict.RData"))
save(glm_formula, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_stepAIC_selected_formula.RData"))

# Plot observed and predicted crossovers for regionsGR grouped into hexiles

# ylims 0 to max
# meanSNPsPB hexiles
sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))]
#[1] 0.002005495 0.004004255 0.006650165 0.010117647 0.015663462 0.062791667
levels(cut(x = dat$meanSNPsPB,
           breaks = c(min(dat$meanSNPsPB, na.rm = T),
                      sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
#[1] "(1.04e-07,0.00201]" "(0.00201,0.004]"    "(0.004,0.00665]"
#[4] "(0.00665,0.0101]"   "(0.0101,0.0157]"    "(0.0157,0.0628]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSNPsPB,
                      breaks = c(min(dat$meanSNPsPB, na.rm = T),
                                 sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_",
           winSize/1e3, "kbWin_", stepSize, "bpStep_SNPsPerBase_hexiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SNPs per base hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# width hexiles
sort(dat$width)[round(1:6*(nrow(dat)/6))]
#[1]     53     99    168    292    595 127250
levels(cut(x = dat$width,
           breaks = c(min(dat$width, na.rm = T),
                      sort(dat$width)[round(1:6*(nrow(dat)/6))])))
#[1] "(4,53]"         "(53,99]"        "(99,168]"       "(168,292]"
#[5] "(292,595]"      "(595,1.27e+05]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$width,
                      breaks = c(min(dat$width, na.rm = T),
                                 sort(dat$width)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_width_hexiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Width hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanSPO11 hexiles
sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))]
#[1] -1.231059665 -0.740196379 -0.383552660 -0.004741569  0.465551989 [6]  4.399179154
levels(cut(x = dat$meanSPO11,
           breaks = c(min(dat$meanSPO11, na.rm = T),
                      sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.47,-1.23]"     "(-1.23,-0.74]"     "(-0.74,-0.384]"
#[4] "(-0.384,-0.00474]" "(-0.00474,0.466]"  "(0.466,4.4]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSPO11,
                      breaks = c(min(dat$meanSPO11, na.rm = T),
                                 sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_SPO11_hexiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SPO11-1-oligo hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanMNase hexiles
sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))]
#[1] -1.16445198 -0.50423235 -0.01477396  0.39903072  0.82942352  3.26628143
levels(cut(x = dat$meanMNase,
           breaks = c(min(dat$meanMNase, na.rm = T),
                      sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.76,-1.16]"    "(-1.16,-0.504]"   "(-0.504,-0.0148]" "(-0.0148,0.399]"
#[5] "(0.399,0.829]"    "(0.829,3.27]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanMNase,
                      breaks = c(min(dat$meanMNase, na.rm = T),
                                 sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_MNase_hexiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Nucleosomes hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanH3K4me3 hexiles
sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))]
#[1] -1.8795789 -1.3227856 -0.7459590 -0.1397163  0.5624519  4.1235247
levels(cut(x = dat$meanH3K4me3,
           breaks = c(min(dat$meanH3K4me3, na.rm = T),
                      sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
#[1] "(-5.03,-1.88]"  "(-1.88,-1.32]"  "(-1.32,-0.746]" "(-0.746,-0.14]"
#[5] "(-0.14,0.562]"  "(0.562,4.12]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanH3K4me3,
                      breaks = c(min(dat$meanH3K4me3, na.rm = T),
                                 sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_H3K4me3_hexiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "H3K4me3 hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanDNAmeth quintiles
sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))]
#[1] 0.000000000 0.002380967 0.014949459 0.204084944 0.450664333 1.000000000
levels(cut(x = dat$meanDNAmeth,
           breaks = c(
                      sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
#[1] "(0,0.00238]"      "(0.00238,0.0149]" "(0.0149,0.204]"   "(0.204,0.451]"
#[5] "(0.451,1]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanDNAmeth,
                      breaks = c(
                                 sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_DNAmeth_quintiles_CO_boxplot_ymin0.pdf"))
boxplot(
  pco,
  ylim = c(0, max(c(unlist(pco), tco), na.rm = T)),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "DNA methylation quintiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# ylims 0 to 60
# meanSNPsPB hexiles
sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))]
#[1] 0.002005495 0.004004255 0.006650165 0.010117647 0.015663462 0.062791667
levels(cut(x = dat$meanSNPsPB,
           breaks = c(min(dat$meanSNPsPB, na.rm = T),
                      sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
#[1] "(1.04e-07,0.00201]" "(0.00201,0.004]"    "(0.004,0.00665]"
#[4] "(0.00665,0.0101]"   "(0.0101,0.0157]"    "(0.0157,0.0628]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSNPsPB,
                      breaks = c(min(dat$meanSNPsPB, na.rm = T),
                                 sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_",
           winSize/1e3, "kbWin_", stepSize, "bpStep_SNPsPerBase_hexiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SNPs per base hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# width hexiles
sort(dat$width)[round(1:6*(nrow(dat)/6))]
#[1]     53     99    168    292    595 127250
levels(cut(x = dat$width,
           breaks = c(min(dat$width, na.rm = T),
                      sort(dat$width)[round(1:6*(nrow(dat)/6))])))
#[1] "(4,53]"         "(53,99]"        "(99,168]"       "(168,292]"
#[5] "(292,595]"      "(595,1.27e+05]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$width,
                      breaks = c(min(dat$width, na.rm = T),
                                 sort(dat$width)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_width_hexiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Width hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanSPO11 hexiles
sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))]
#[1] -1.231059665 -0.740196379 -0.383552660 -0.004741569  0.465551989 [6]  4.399179154
levels(cut(x = dat$meanSPO11,
           breaks = c(min(dat$meanSPO11, na.rm = T),
                      sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.47,-1.23]"     "(-1.23,-0.74]"     "(-0.74,-0.384]"
#[4] "(-0.384,-0.00474]" "(-0.00474,0.466]"  "(0.466,4.4]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSPO11,
                      breaks = c(min(dat$meanSPO11, na.rm = T),
                                 sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_SPO11_hexiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SPO11-1-oligo hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanMNase hexiles
sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))]
#[1] -1.16445198 -0.50423235 -0.01477396  0.39903072  0.82942352  3.26628143
levels(cut(x = dat$meanMNase,
           breaks = c(min(dat$meanMNase, na.rm = T),
                      sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.76,-1.16]"    "(-1.16,-0.504]"   "(-0.504,-0.0148]" "(-0.0148,0.399]"
#[5] "(0.399,0.829]"    "(0.829,3.27]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanMNase,
                      breaks = c(min(dat$meanMNase, na.rm = T),
                                 sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_MNase_hexiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Nucleosomes hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanH3K4me3 hexiles
sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))]
#[1] -1.8795789 -1.3227856 -0.7459590 -0.1397163  0.5624519  4.1235247
levels(cut(x = dat$meanH3K4me3,
           breaks = c(min(dat$meanH3K4me3, na.rm = T),
                      sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
#[1] "(-5.03,-1.88]"  "(-1.88,-1.32]"  "(-1.32,-0.746]" "(-0.746,-0.14]"
#[5] "(-0.14,0.562]"  "(0.562,4.12]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanH3K4me3,
                      breaks = c(min(dat$meanH3K4me3, na.rm = T),
                                 sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_H3K4me3_hexiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "H3K4me3 hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanDNAmeth quintiles
sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))]
#[1] 0.000000000 0.002380967 0.014949459 0.204084944 0.450664333 1.000000000
levels(cut(x = dat$meanDNAmeth,
           breaks = c(
                      sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
#[1] "(0,0.00238]"      "(0.00238,0.0149]" "(0.0149,0.204]"   "(0.204,0.451]"
#[5] "(0.451,1]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanDNAmeth,
                      breaks = c(
                                 sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
pco <- lapply(ssID, function(x) {
  sapply(1:100, function(ii) {
    sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
    sum( dat$width[x] ) * 1e6
  })
})
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
pdf(paste0(plotDir, "GLM_binomial_logit_", winSize/1e3, "kbWin_", stepSize, "bpStep_DNAmeth_quintiles_CO_boxplot_ymin0_ymax60.pdf"))
boxplot(
  pco,
  ylim = c(0, 60),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "DNA methylation quintiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()
