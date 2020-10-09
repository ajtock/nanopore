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

library(fitdistrplus) # descdist, plotdist, fitdist included
library(glm2)
library(MASS) # glm.nb included; MASS is also loaded as a dependency of fitdistrplus
library(pscl) # zeroinfl included
library(vcd) # goodfit included
library(qualityTools) # qqPlot included
library(stats4) # mle (for estimating parameters by maximum likelihood) included
#library(segmentSeq)
#library(GenomicRanges)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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

# Inspect distribution of wt_cMMb
# By plotting the empirical density and the empirical cumulative distribution function (ECDF)
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plotdist.pdf"))
fitdistrplus::plotdist(winDF$wt_cMMb, histo = T, demp = T)
dev.off()

# Plot with offset of +1
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_plotdist.pdf"))
fitdistrplus::plotdist(winDF$wt_cMMb+1, histo = T, demp = T)
dev.off()

# Plot with offset of +1e-06
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_plotdist.pdf"))
fitdistrplus::plotdist(winDF$wt_cMMb+1e-06, histo = T, demp = T)
dev.off()

# Plot random gamma distribution
set.seed(9382)
pdf(paste0(plotDir, "random_gamma_n95_plotdist.pdf"))
fitdistrplus::plotdist(rgamma(n = length(winDF$wt_cMMb), shape = 1, rate = 0.25),
                       histo = T, demp = T)
dev.off()

# Plot random exponential distribution
set.seed(9382)
pdf(paste0(plotDir, "random_exponential_n95_plotdist.pdf"))
fitdistrplus::plotdist(rexp(n = length(winDF$wt_cMMb), rate = 0.25),
                       histo = T, demp = T)
dev.off()

# Plot random weibull distribution
set.seed(9382)
pdf(paste0(plotDir, "random_weibull_n95_plotdist.pdf"))
fitdistrplus::plotdist(rweibull(n = length(winDF$wt_cMMb), shape = 1, scale = 4),
                       histo = T, demp = T)
dev.off()

# By plotting a skewness-kurtosis (Cullen and Frey) graph
# using descdist from the fitdistrplus package, with 1000 
# See https://stats.stackexchange.com/questions/58220/what-distribution-does-my-data-follow
# and https://stackoverflow.com/questions/31741742/how-to-identify-the-distribution-of-the-given-data-using-r
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_descdistPlot.pdf"))
fitdistrplus::descdist(winDF$wt_cMMb, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  0   max:  19.46348 
#median:  2.897937 
#mean:  4.516331 
#estimated sd:  4.750273 
#estimated skewness:  1.301378 
#estimated kurtosis:  4.201431
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_descdistPlot.pdf"))
fitdistrplus::descdist(winDF$wt_cMMb+1e-06, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  1e-06   max:  19.46348 
#median:  2.897938 
#mean:  4.516332 
#estimated sd:  4.750273 
#estimated skewness:  1.301378 
#estimated kurtosis:  4.201431 
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_descdistPlot.pdf"))
fitdistrplus::descdist(winDF$wt_cMMb+1, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  1   max:  20.46348
#median:  3.897937
#mean:  5.516331
#estimated sd:  4.750273
#estimated skewness:  1.301378
#estimated kurtosis:  4.201431

# Based on CTL3.9_wt_cMMb_descdistPlot.pdf, inspect fit to distributions
# Quantile-quantile (qq) plot for gamma distribution, with offset of +1
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_qqPlot_gamma.pdf"))
qualityTools::qqPlot(x = winDF$wt_cMMb+1, y = "gamma")
dev.off()
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_qqPlot_gamma.pdf"))
qualityTools::qqPlot(x = winDF$wt_cMMb+1e-06, y = "gamma")
dev.off()
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_qqPlot_weibull.pdf"))
qualityTools::qqPlot(x = winDF$wt_cMMb+1, y = "weibull")
dev.off()
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_qqPlot_plus1_exponential.pdf"))
qualityTools::qqPlot(x = winDF$wt_cMMb+1, y = "exponential")
dev.off()
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_qqPlot_exponential.pdf"))
qualityTools::qqPlot(x = winDF$wt_cMMb, y = "exponential")
dev.off()

# On the gamma distribution, see https://bookdown.org/probability/beta/beta-and-gamma.html

# Inspect fit of distribution using maximum likelihood estimation (MLE) within fitdistrplus, with offset of +1
# See https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
wt_cMMb_plus1_fitdist_gamma <- fitdistrplus::fitdist(data = winDF$wt_cMMb+1,
                                                          distr = "gamma",
                                                          method = "mle")
wt_cMMb_plus1_fitdist_lnorm <- fitdistrplus::fitdist(data = winDF$wt_cMMb+1,
                                                          distr = "lnorm",
                                                          method = "mle")
wt_cMMb_plus1_fitdist_weibull <- fitdistrplus::fitdist(data = winDF$wt_cMMb+1,
                                                            distr = "weibull",
                                                            method = "mle")
wt_cMMb_plus1_fitdist_exp <- fitdistrplus::fitdist(data = winDF$wt_cMMb+1,
                                                        distr = "exp",
                                                        method = "mle")
wt_cMMb_fitdist_exp <- fitdistrplus::fitdist(data = winDF$wt_cMMb,
                                             distr = "exp",
                                             method = "mle")
print(summary(wt_cMMb_plus1_fitdist_gamma))
print(summary(wt_cMMb_plus1_fitdist_lnorm))
print(summary(wt_cMMb_plus1_fitdist_weibull))
print(summary(wt_cMMb_plus1_fitdist_exp))
print(summary(wt_cMMb_fitdist_exp))

# Generate goodness-of-fit plots:
# See https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
# 1. "a density plot representing the density function of the fitted distribution
#     along with the histogram of the empirical distribution"
# 2. "a cumulative distribution function (CDF) plot of both the
#     empirical distribution and the fitted distribution"
# 3. "a Q-Q plot representing the empirical quantiles (y-axis)
#     against the theoretical quantiles (x-axis)"
# 4. "a P-P plot representing the empirical distribution function evaluated at each data point (y-axis)
#     against the fitted distribution function (x-axis)"

# Load modified denscomp and cdfcomp functions (allowing line widths to be changed)
load("fitdistrplus_denscomp_mod.RData")
load("fitdistrplus_cdfcomp_mod.RData")
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Log-normal", "Weibull", "Exponential", "Gamma")
denscomp(list(wt_cMMb_plus1_fitdist_lnorm, wt_cMMb_plus1_fitdist_weibull,
              wt_cMMb_plus1_fitdist_exp, wt_cMMb_plus1_fitdist_gamma),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(wt_cMMb_plus1_fitdist_lnorm, wt_cMMb_plus1_fitdist_weibull,
            wt_cMMb_plus1_fitdist_exp, wt_cMMb_plus1_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(wt_cMMb_plus1_fitdist_lnorm, wt_cMMb_plus1_fitdist_weibull,
             wt_cMMb_plus1_fitdist_exp, wt_cMMb_plus1_fitdist_gamma),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(wt_cMMb_plus1_fitdist_lnorm, wt_cMMb_plus1_fitdist_weibull,
            wt_cMMb_plus1_fitdist_exp, wt_cMMb_plus1_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
dev.off()
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_plus1eNeg6_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Log-normal", "Weibull", "Exponential", "Gamma")
denscomp(list(wt_cMMb_plus1eNeg6_fitdist_lnorm, wt_cMMb_plus1eNeg6_fitdist_weibull,
              wt_cMMb_plus1eNeg6_fitdist_exp, wt_cMMb_plus1eNeg6_fitdist_gamma),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(wt_cMMb_plus1eNeg6_fitdist_lnorm, wt_cMMb_plus1eNeg6_fitdist_weibull,
            wt_cMMb_plus1eNeg6_fitdist_exp, wt_cMMb_plus1eNeg6_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(wt_cMMb_plus1eNeg6_fitdist_lnorm, wt_cMMb_plus1eNeg6_fitdist_weibull,
             wt_cMMb_plus1eNeg6_fitdist_exp, wt_cMMb_plus1eNeg6_fitdist_gamma),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(wt_cMMb_plus1eNeg6_fitdist_lnorm, wt_cMMb_plus1eNeg6_fitdist_weibull,
            wt_cMMb_plus1eNeg6_fitdist_exp, wt_cMMb_plus1eNeg6_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
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
                  cMMb = winDF$wt_cMMb)

# First fit model using the gamma distribution,
# which will enable model fitting using the exponential distribution
#glmGamma <- glm2(formula = cMMb+1 ~ (CENH3_Rep1 + H2A_Set1 + H2AW6_Set1 + H2AX_Set1 +
#                                     H2AZ_Set1 + H3_Set1 + H2AW7_Set2 + H3_Set2 +
#                                     H3K27me3_Set2 + H3K4me1_Set2 + H3K4me3_Set2 + H4K20me1_Set2 +
#                                     H1_Set3 + H3_Set3 + H3_Set5 + H3K27me1_Set5 +
#                                     H3K36me3_Set5 + H3K9me1_Set5 + H3K9me2_Set5 + H3_Set7 +
#                                     HTB1_Set7 + HTB2_Set7 + HTB3_Set7 + HTB4_Set7 +
#                                     H3K9me2_Rep1 + REC8_HA_Rep2 + ASY1_Rep1 + MSH4_HA_Rep1 +
#                                     DMC1_V5_Rep1 + MTOPVIB_HA_Rep1 + MTOPVIB_HA_Rep2 + SPO11oligos_Rep1 +
#                                     width)^2,
#                 data = dat,
#                 family = Gamma(link="inverse"),
#                 control = glm.control(maxit = 100000))
glmGamma <- glm2(formula = cMMb+1 ~ (CENH3_Rep1 + H2AW6_Set1 + H2AZ_Set1 + H3K4me1_Set2 + H3K4me3_Set2 +
                                     H3K9me2_Rep1 + REC8_HA_Rep2 + ASY1_Rep1 + MSH4_HA_Rep1 +
                                     DMC1_V5_Rep1 + MTOPVIB_HA_Rep1 + MTOPVIB_HA_Rep2 + SPO11oligos_Rep1)^2,
                 data = dat,
                 family = Gamma(link="inverse"),
                 control = glm.control(maxit = 100000))
# See https://stat.ethz.ch/pipermail/r-help/2003-June/034852.html
# and https://stats.stackexchange.com/questions/240455/fitting-exponential-regression-model-by-mle
# and https://stats.stackexchange.com/questions/250077/how-do-you-specify-exponential-distribution-in-glm-in-r?noredirect=1&lq=1
# "The Gamma family is parametrised in glm() by two parameters: 
#  mean and dispersion; the "dispersion" regulates the shape.
#  "So must fit a GLM with the Gamma family, and then produce a 'summary'
#  with dispersion parameter set equal to 1, since this value 
#  corresponds to the exponential distribution in the Gamma family.
#  
#  "In practice:
#  
#  fit <- glm(formula =...,  family = Gamma)
#  summary(fit,dispersion=1)
#
#  "
summary(glmGamma)
# Summary for model using exponential distribution
summary(glmGamma, dispersion = 1)

glmGamma_stepAIC <- stepAIC(object = glmGamma, direction = "both")
print("stepAIC-selected model formula:")
print(formula(glmGamma_stepAIC))
#print(glmGamma_stepAIC$formula)
glmGamma_formula <- formula(glmGamma_stepAIC)
glmGamma_select <- glm2(formula = formula(glmGamma_stepAIC),
                        data = dat,
                        family = Gamma(link="inverse"),
                        control = glm.control(maxit = 100000))
# Estimate the shape parameter and adjust GLM coefficient estimates and predictions
# See https://stats.stackexchange.com/questions/58497/using-r-for-glm-with-gamma-distribution
glmGamma_shape <- gamma.shape(glmGamma_select)
glmGamma_predict <- predict(glmGamma_select, type = "response",
                            se = T, dispersion = 1/glmGamma_shape$alpha)
glmGamma_summary <- summary(glmGamma_select, dispersion = 1/glmGamma_shape$alpha)
glmGamma_coeffs <- glmGamma_summary$coefficients
glmExp_predict <- predict(glmGamma_select, type = "response",
                          se = T, dispersion = 1)
glmExp_summary <- summary(glmGamma_select, dispersion = 1)
glmExp_coeffs <- glmExp_summary$coefficients

save(glmGamma_stepAIC, file = "CTL3.9_wt_cMMb_glmGamma_stepAIC.RData")
save(glmGamma_formula, file = "CTL3.9_wt_cMMb_glmGamma_formula.RData")
save(glmGamma_select, file = "CTL3.9_wt_cMMb_glmGamma_select.RData")
save(glmGamma_shape, file = "CTL3.9_wt_cMMb_glmGamma_shape.RData")
save(glmGamma_predict, file = "CTL3.9_wt_cMMb_glmGamma_predict.RData")
save(glmGamma_summary, file = "CTL3.9_wt_cMMb_glmGamma_summary.RData")
save(glmGamma_coeffs, file = "CTL3.9_wt_cMMb_glmGamma_coeffs.RData")
save(glmExp_predict, file = "CTL3.9_wt_cMMb_glmExp_predict.RData")
save(glmExp_summary, file = "CTL3.9_wt_cMMb_glmExp_summary.RData")
save(glmExp_coeffs, file = "CTL3.9_wt_cMMb_glmExp_coeffs.RData")

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# P-value > 0.05 indicates failure to reject the null hypothesis
# that the model fits the data
1 - pchisq(q = summary(glmGamma_select, dispersion = 1/glmGamma_shape$alpha)$deviance,
           df = summary(glmGamma_select, dispersion = 1/glmGamma_shape$alpha)$df.residual)
#[1] 1
1 - pchisq(q = summary(glmGamma_select, dispersion = 1)$deviance,
           df = summary(glmGamma_select, dispersion = 1)$df.residual)
#[1] 1

# Disable scientific notation for plotting
options("scipen"=100)

# Plot observed cM/Mb and predicted cM/Mb by GLM (glmGamma_select)
pdf(paste0(plotDir, "CTL3.9_wt_cMMb_observed_predicted_gamma_GLM.pdf"),
    height = 5, width = 10)
plot(x = round((winDF$start+winDF$end)/2),
     y = glmGamma_select$y-1,
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     type = "l", lwd = 2, col = "red")
lines(x = round((winDF$start+winDF$end)/2),
      y = glmGamma_predict$fit-1,
      lty = 2, lwd = 2, col = "black")
axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
     at = round((winDF$start+winDF$end)/2),
     labels = round(round((winDF$start+winDF$end)/2)/1e6, digits = 3))
axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
mtext(side = 1, line = 2, cex = 1, text = "Coordinates (Mb)")
mtext(side = 2, line = 2, cex = 1, text = "cM/Mb")
mtext(side = 3, line = 2, cex = 0.75,
      text = "Gamma GLM: cM/Mb~(CENH3+H2A.W6+H2A.Z+H3K4me1+H3K4me3+H3K9me2+REC8+ASY1+MSH4+DMC1+MTOPVIBr1+MTOPVIBr2+SPO11-1)^2")
box(lwd = 1.5)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"),
       col = "white",
       ncol = 1, cex = 1, lwd = 1.5, bty = "n")
dev.off()

## For modelling CO counts, we would need to use log2(sum ChIP/sum input) coverage,
## rather than log2(mean ChIP/mean input) coverage in each windows
## COs data show overdispersion with unequal means and variances
## Therefore, evaluate fit of negative binomial GLM with "log" link function
#nbinom <- glm.nb(formula = COs ~ (CENH3_Rep1 +
#                                  H2AZ_Set1 +
#                                  H3K9me2_Rep1 + REC8_HA_Rep2 + ASY1_Rep1 + MSH4_HA_Rep1 +
#                                  DMC1_V5_Rep1 +
##                                  MTOPVIB_HA_Rep1 +
#                                  MTOPVIB_HA_Rep2 +
#                                  SPO11oligos_Rep1)^2,
#                 data = dat,
#                 control = glm.control(maxit = 100000),
#                 link = log)
