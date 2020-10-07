#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock (adapted from original script logistic2.R by Tom Hardcastle)
# contact: ajt200@cam.ac.uk
# date: 07.10.2020

# Build GLM for the CTL3.9 interval within RaGOO v2.0 (assembled by Matt Naish),
# with recombination rate (cM/Mb; obtained by Joiselle) as the response variable
# and chromatin and recombination ChIP-seq signals as predictor variables

# Usage:
# Rscript ./GLM_v071020.R

options(stringsAsFactors = FALSE)

library(fitdistrplus)
library(glm2)
library(MASS) # glm.nb() included; MASS is also loaded as a dependency of fitdistrplus
library(pscl) # zeroinfl() included
library(vcd) # goodfit included
#library(segmentSeq)
#library(GenomicRanges)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
#library(segmentSeq)
#library(GenomicRanges)

# Load cM/Mb data within CTL3.9 interval
dat <- read.csv("ragoo2.cMMb.final.csv", header = T)
# Convert into data.frame of marker intervals
winDF <- data.frame()
for(x in 1:(dim(dat)[1]-1)) {
  rowDF <- data.frame(chr = "Chr3",
                      start = dat[x,]$ragoo2.markers,
                      end = dat[x+1,]$ragoo2.markers,
                      width = (dat[x+1,]$ragoo2.markers-dat[x,]$ragoo2.markers)+1,
                      wt_COs = dat[x,]$wt.cos,
                      wt_cMMb = dat[x,]$wt.cMMb,
                      cmt3_COs = dat[x,]$cmt3.cos,
                      cmt3_cMMb = dat[x,]$cmt3.cMMb,
                      met1_COs = dat[x,]$met1.cos,
                      met1_cMMb = dat[x,]$met1.cMMb)
  winDF <- rbind(winDF, rowDF)
}
#winGR <- GRanges(seqnames = winDF$chr,
#                 ranges = IRanges(winDF$start,
#                                  winDF$end),
#                 strand = "*",
#                 winDF[,5:10])

# Inspect distribution of wt_cMMb by plotting a skewness-kurtosis (Cullen and Frey) graph
# using descdist from the fitdistrplus package, with 1000 
# See https://stats.stackexchange.com/questions/58220/what-distribution-does-my-data-follow
# and https://stackoverflow.com/questions/31741742/how-to-identify-the-distribution-of-the-given-data-using-r
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_descdistPlot.pdf"))
CFplot <- descdist(winDF$wt_cMMb, discrete = T,
                   obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()

# Test goodness-of-fit to the negative binomial distribution using goodfit from vcd
# A P-value > 0.05 indicates that the model fits the data
print(summary(goodfit(winDF$wt_cMMb, type = "nbinomial", method = "MinChisq")))
print(summary(goodfit(winDF$wt_cMMb, type = "nbinomial", method = "ML")))
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_goodness_of_fit_nbinomial.pdf"))
plot(goodfit(winDF$wt_cMMb, type = "nbinomial", method = "ML"))
dev.off()

# Load log2(ChIP/input) coverage within marker intervals
covMat <- read.table(paste0("/home/ajt200/analysis/CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/",
                            "xiaohui_pipeline_RaGOO_v2.0/coverage/log2ChIPinput/CTL3.9profiles/",
                            "log2ChIPinput_RaGOO_v2.0_CTL3.9_marker_intervals_norm_coverage_matrix_unsmoothed.tsv"),
                            header = T)
covMat_colnames <- colnames(covMat)
covMat_colnames <- gsub("log2_WT_", "", covMat_colnames)
covMat_colnames <- gsub("log2_", "", covMat_colnames)
covMat_colnames <- gsub("ChIP_set(\\d)", "Set\\1_ChIP", covMat_colnames, perl = T)
covMat_colnames <- gsub("input_set(\\d)", "Set\\1_input", covMat_colnames, perl = T)
covMat_colnames <- gsub("_ChIP.*", "", covMat_colnames, perl = T) 
covMat_colnames <- gsub("_WT_gDNA_Rep1_R1", "", covMat_colnames) 
colnames(covMat) <- covMat_colnames

# Create data object for model
dat <- data.frame(covMat,
                  cMMb = winDF$wt_cMMb,
                  COs = winDF$wt_COs)

# cM/Mb and COs data show overdispersion with unequal means and variances
# Therefore, build negative binomial GLM with "log" link function
nbinom <- glm.nb(formula = cMMb ~ (CENH3_Rep1 + H2A_Set1 + H2AW6_Set1 + H2AX_Set1 +
                                   H2AZ_Set1 + H3_Set1 + H2AW7_Set2 + H3_Set2 +
                                   H3K27me3_Set2 + H3K4me1_Set2 + H3K4me3_Set2 + H4K20me1_Set2 +
                                   H1_Set3 + H3_Set3 + H3_Set5 + H3K27me1_Set5 +
                                   H3K36me3_Set5 + H3K9me1_Set5 + H3K9me2_Set5 + H3_Set7 +
                                   HTB1_Set7 + HTB2_Set7 + HTB3_Set7 + HTB4_Set7 +
                                   H3K9me2_Rep1 + REC8_HA_Rep2 + ASY1_Rep1 + MSH4_HA_Rep1 +
                                   DMC1_V5_Rep1 + MTOPVIB_HA_Rep1 + MTOPVIB_HA_Rep2 + SPO11oligos_Rep1 +
                                   width)^2,
                 data = dat,
                 control = glm.control(maxit = 100000),
                 link = log)
print(warnings())
nbinom_stepAIC <- stepAIC(object = nbinom, direction = "both")
print(warnings())
print("stepAIC-selected model formula:")
print(formula(nbinom_stepAIC)) # ????
#print(nbinom_stepAIC$formula)
nbinom_select <- glm.nb(formula = formula(nbinom_stepAIC),
                        data = dat,
                        control = glm.control(maxit = 100000),
                        link = log)
nbinom_summary <- summary(nbinom_select)
nbinom_coeffs <- nbinom_summary$coefficients
nbinom_predict <- predict(nbinom_select, type = "response")
nbinom_formula <- nbinom_select$formula
save(nbinom_stepAIC, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_stepAIC.RData"))
save(nbinom_select, file = paste0(outDir, "GLM_binomial_logit_",
                               winSize/1e3, "kbWin_", stepSize, "bpStep.RData"))
save(nbinom_summary, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_summary.RData"))
write.csv(nbinom_coeffs, file = paste0(outDir, "GLM_binomial_logit_",
                                    winSize/1e3, "kbWin_", stepSize, "bpStep_coeff.csv"))
save(nbinom_predict, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_predict.RData"))
save(nbinom_formula, file = paste0(outDir, "GLM_binomial_logit_",
                                winSize/1e3, "kbWin_", stepSize, "bpStep_stepAIC_selected_formula.RData"))


## Zero-inflated negative binomial regression:
#ZINB <- lapply(seq_along(pops_winIndCOs_list), function(x) {
#  zeroinfl(formula = pops_winIndCOs_list[[x]] ~ genotype | 1,
#           dist = "negbin")
#})

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
