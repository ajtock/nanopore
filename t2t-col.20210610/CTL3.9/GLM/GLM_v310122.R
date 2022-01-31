#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 31.01.2022

# Build GLM with recombination rate (cM/Mb; obtained by Joiselle) as the response variable
# and chromatin and recombination ChIP-seq signals (and others) as predictor variables

# Usage:
# Rscript ./GLM_v310122.R

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

# Load cM/Mb data
dat <- read.csv("glm.dat.csv", header = T)

# Inspect distribution of cMMb:
# 1. by plotting the empirical density and the empirical cumulative distribution function (ECDF)
pdf(paste0(plotDir, "cMMb_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb, histo = T, demp = T)
dev.off()

# Plot with offset of +1
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "cMMb_plus1_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb+1, histo = T, demp = T)
dev.off()

# Plot with offset of +1e-06
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "cMMb_plus1e-06_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb+1e-06, histo = T, demp = T)
dev.off()

# Plot random gamma distribution
set.seed(9382)
pdf(paste0(plotDir, "random_gamma_n", length(dat$cMMb), "_plotdist.pdf"))
fitdistrplus::plotdist(rgamma(n = length(dat$cMMb), shape = 1, rate = 0.25),
                       histo = T, demp = T)
dev.off()

# Plot random exponential distribution
set.seed(9382)
pdf(paste0(plotDir, "random_exponential_n", length(dat$cMMb), "_plotdist.pdf"))
fitdistrplus::plotdist(rexp(n = length(dat$cMMb), rate = 0.25),
                       histo = T, demp = T)
dev.off()

# Plot random weibull distribution
set.seed(9382)
pdf(paste0(plotDir, "random_weibull_n", length(dat$cMMb), "_plotdist.pdf"))
fitdistrplus::plotdist(rweibull(n = length(dat$cMMb), shape = 1, scale = 4),
                       histo = T, demp = T)
dev.off()

# 2. by plotting a skewness-kurtosis (Cullen and Frey) graph
# using descdist from the fitdistrplus package, with 1000 bootstraps 
# See https://stats.stackexchange.com/questions/58220/what-distribution-does-my-data-follow
# and https://stackoverflow.com/questions/31741742/how-to-identify-the-distribution-of-the-given-data-using-r
pdf(paste0(plotDir, "cMMb_descdistPlot.pdf"))
fitdistrplus::descdist(dat$cMMb, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  0   max:  19.47513
#median:  2.900039
#mean:  4.510157
#estimated sd:  4.742442
#estimated skewness:  1.303429
#estimated kurtosis:  4.211326

pdf(paste0(plotDir, "cMMb_descdistPlot_discrete.pdf"))
fitdistrplus::descdist(dat$cMMb, discrete = T,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  0   max:  19.47513 
#median:  2.900039 
#mean:  4.510157 
#estimated sd:  4.742442 
#estimated skewness:  1.303429 
#estimated kurtosis:  4.211326

pdf(paste0(plotDir, "cMMb_plus1e-06_descdistPlot.pdf"))
fitdistrplus::descdist(dat$cMMb+1e-06, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  1e-06   max:  19.47513
#median:  2.90004
#mean:  4.510158
#estimated sd:  4.742442
#estimated skewness:  1.303429
#estimated kurtosis:  4.211326

pdf(paste0(plotDir, "cMMb_plus1_descdistPlot.pdf"))
fitdistrplus::descdist(dat$cMMb+1, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  1   max:  20.47513
#median:  3.900039
#mean:  5.510157
#estimated sd:  4.742442
#estimated skewness:  1.303429
#estimated kurtosis:  4.211326


# Based on cMMb_descdistPlot.pdf, inspect fit to distributions
# Quantile-quantile (qq) plot for gamma distribution, with offset of +1 or +1e-06
pdf(paste0(plotDir, "cMMb_plus1_qqPlot_gamma.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1, y = "gamma")
dev.off()
pdf(paste0(plotDir, "cMMb_plus1e-06_qqPlot_gamma.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1e-06, y = "gamma")
dev.off()
#There were 12 warnings (use warnings() to see them)
#Warning messages:
#1: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#2: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#3: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#4: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#5: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#6: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#7: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#8: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#9: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#10: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#11: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#12: In densfun(x, parm[1], parm[2], ...) : NaNs produced

# Quantile-quantile (qq) plot for Weibull distribution, with offset of +1 or +1e-06
pdf(paste0(plotDir, "cMMb_plus1_qqPlot_weibull.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1, y = "weibull")
dev.off()
#Warning messages:
#1: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#2: In densfun(x, parm[1], parm[2], ...) : NaNs produced
pdf(paste0(plotDir, "cMMb_plus1e-06_qqPlot_weibull.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1e-06, y = "weibull")
dev.off()
#Warning messages:
#1: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#pdf(paste0(plotDir, "cMMb_qqPlot_weibull.pdf"))
#qualityTools::qqPlot(x = dat$cMMb, y = "weibull")
#dev.off()
#Error in (function (x, densfun, start, ...)  : Weibull values must be > 0 

# Quantile-quantile (qq) plot for exponential distribution, with offset of +1, or +1e-06 or no offset
pdf(paste0(plotDir, "cMMb_qqPlot_plus1_exponential.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1, y = "exponential")
dev.off()
pdf(paste0(plotDir, "cMMb_qqPlot_plus1e-06_exponential.pdf"))
qualityTools::qqPlot(x = dat$cMMb+1e-06, y = "exponential")
dev.off()
pdf(paste0(plotDir, "cMMb_qqPlot_exponential.pdf"))
qualityTools::qqPlot(x = dat$cMMb, y = "exponential")
dev.off()

# On the gamma distribution, see https://bookdown.org/probability/beta/beta-and-gamma.html

# Inspect fit of distribution using maximum likelihood estimation (MLE) within fitdistrplus, with offset of +1
# See https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
cMMb_plus1_fitdist_gamma <- fitdistrplus::fitdist(data = dat$cMMb+1,
                                                  distr = "gamma",
                                                  method = "mle")
cMMb_plus1_fitdist_lnorm <- fitdistrplus::fitdist(data = dat$cMMb+1,
                                                  distr = "lnorm",
                                                  method = "mle")
cMMb_plus1_fitdist_weibull <- fitdistrplus::fitdist(data = dat$cMMb+1,
                                                    distr = "weibull",
                                                    method = "mle")
cMMb_plus1_fitdist_exp <- fitdistrplus::fitdist(data = dat$cMMb+1,
                                                distr = "exp",
                                                method = "mle")
cMMb_fitdist_exp <- fitdistrplus::fitdist(data = dat$cMMb,
                                          distr = "exp",
                                          method = "mle")

cMMb_plus1eNeg06_fitdist_gamma <- fitdistrplus::fitdist(data = dat$cMMb+1e-06,
                                                        distr = "gamma",
                                                        method = "mle")
cMMb_plus1eNeg06_fitdist_lnorm <- fitdistrplus::fitdist(data = dat$cMMb+1e-06,
                                                        distr = "lnorm",
                                                        method = "mle")
cMMb_plus1eNeg06_fitdist_weibull <- fitdistrplus::fitdist(data = dat$cMMb+1e-06,
                                                          distr = "weibull",
                                                          method = "mle")
cMMb_plus1eNeg06_fitdist_exp <- fitdistrplus::fitdist(data = dat$cMMb+1e-06,
                                                      distr = "exp",
                                                      method = "mle")
print(summary(cMMb_plus1_fitdist_gamma))
print(summary(cMMb_plus1_fitdist_lnorm))
print(summary(cMMb_plus1_fitdist_weibull))
print(summary(cMMb_plus1_fitdist_exp))
print(summary(cMMb_fitdist_exp))

print(summary(cMMb_plus1eNeg06_fitdist_gamma))
print(summary(cMMb_plus1eNeg06_fitdist_lnorm))
print(summary(cMMb_plus1eNeg06_fitdist_weibull))
print(summary(cMMb_plus1eNeg06_fitdist_exp))

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
pdf(paste0(plotDir, "cMMb_plus1_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Log-normal", "Weibull", "Exponential", "Gamma")
denscomp(list(cMMb_plus1_fitdist_lnorm, cMMb_plus1_fitdist_weibull,
              cMMb_plus1_fitdist_exp, cMMb_plus1_fitdist_gamma),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(cMMb_plus1_fitdist_lnorm, cMMb_plus1_fitdist_weibull,
            cMMb_plus1_fitdist_exp, cMMb_plus1_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(cMMb_plus1_fitdist_lnorm, cMMb_plus1_fitdist_weibull,
             cMMb_plus1_fitdist_exp, cMMb_plus1_fitdist_gamma),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(cMMb_plus1_fitdist_lnorm, cMMb_plus1_fitdist_weibull,
            cMMb_plus1_fitdist_exp, cMMb_plus1_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
dev.off()

pdf(paste0(plotDir, "cMMb_plus1e-06_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Log-normal", "Weibull", "Exponential", "Gamma")
denscomp(list(cMMb_plus1eNeg06_fitdist_lnorm, cMMb_plus1eNeg06_fitdist_weibull,
              cMMb_plus1eNeg06_fitdist_exp, cMMb_plus1eNeg06_fitdist_gamma),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(cMMb_plus1eNeg06_fitdist_lnorm, cMMb_plus1eNeg06_fitdist_weibull,
            cMMb_plus1eNeg06_fitdist_exp, cMMb_plus1eNeg06_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(cMMb_plus1eNeg06_fitdist_lnorm, cMMb_plus1eNeg06_fitdist_weibull,
             cMMb_plus1eNeg06_fitdist_exp, cMMb_plus1eNeg06_fitdist_gamma),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(cMMb_plus1eNeg06_fitdist_lnorm, cMMb_plus1eNeg06_fitdist_weibull,
            cMMb_plus1eNeg06_fitdist_exp, cMMb_plus1eNeg06_fitdist_gamma),
       legendtext = plot.legend, fitpch = 20)
dev.off()


## Load log2(ChIP/input) coverage within marker intervals
#covMat <- read.table(paste0("/home/ajt200/analysis/CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/",
#                            "xiaohui_pipeline_RaGOO_v2.0/coverage/log2ChIPinput/CTL3.9profiles/",
#                            "log2ChIPinput_RaGOO_v2.0_marker_intervals_norm_coverage_matrix_unsmoothed.tsv"),
#                            header = T)
#covMat_colnames <- colnames(covMat)
#covMat_colnames <- gsub("log2_WT_", "", covMat_colnames)
#covMat_colnames <- gsub("log2_", "", covMat_colnames)
#covMat_colnames <- gsub("ChIP_set(\\d)", "Set\\1_ChIP", covMat_colnames, perl = T)
#covMat_colnames <- gsub("input_set(\\d)", "Set\\1_input", covMat_colnames, perl = T)
#covMat_colnames <- gsub("_ChIP.*", "", covMat_colnames, perl = T) 
#covMat_colnames <- gsub("_WT_gDNA_Rep1_R1", "", covMat_colnames) 
#colnames(covMat) <- covMat_colnames
#
## Create data object for model
#dat <- data.frame(covMat,
#                  cMMb = dat$cMMb)

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
glmGamma <- glm2(formula = cMMb+1 ~ (CG + CHG + CHH + genes + H2AZ + H3K4me3 + H3K9me2 +
                                     REC8 + ASY1 + CENH3 + SPO11 + SNV)^2,
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
stopifnot(identical(formula(glmGamma_select), formula(glmGamma_stepAIC)))
glmGamma_formula <- formula(glmGamma_select)
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

save(glmGamma_stepAIC, file = "cMMb_glmGamma_stepAIC.RData")
save(glmGamma_formula, file = "cMMb_glmGamma_formula.RData")
save(glmGamma_select, file = "cMMb_glmGamma_select.RData")
save(glmGamma_shape, file = "cMMb_glmGamma_shape.RData")
save(glmGamma_predict, file = "cMMb_glmGamma_predict.RData")
save(glmGamma_summary, file = "cMMb_glmGamma_summary.RData")
save(glmGamma_coeffs, file = "cMMb_glmGamma_coeffs.RData")
glmGamma_coeffs_df <- data.frame(glmGamma_coeffs)
colnames(glmGamma_coeffs_df) <- c("Estimate", "StdError", "z-value", "P-value")
write.table(glmGamma_coeffs_df,
            file = "cMMb_glmGamma_coeffs.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

save(glmExp_predict, file = "cMMb_glmExp_predict.RData")
save(glmExp_summary, file = "cMMb_glmExp_summary.RData")
save(glmExp_coeffs, file = "cMMb_glmExp_coeffs.RData")

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
pdf(paste0(plotDir, "cMMb_observed_predicted_gamma_GLM.pdf"),
    height = 5, width = 10)
#plot(x = round((dat$start+dat$end)/2),
plot(x = dat$X,
     y = glmGamma_select$y-1,
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     type = "l", lwd = 2, col = "red")
#lines(x = round((dat$start+dat$end)/2),
lines(x = dat$X,
      y = glmGamma_predict$fit-1,
      lty = 2, lwd = 2, col = "black")
axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
#     at = round((dat$start+dat$end)/2),
#     labels = round(round((dat$start+dat$end)/2)/1e6, digits = 3))
     at = dat$X,
     labels = dat$X)
axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
mtext(side = 1, line = 2, cex = 1, text = "Genomic window")
mtext(side = 2, line = 2, cex = 1, text = "cM/Mb")
mtext(side = 3, line = 2, cex = 0.75,
      text = "Gamma GLM: cM/Mb~(mCG+mCHG+mCHH+Genes+H2A.Z+H3K4me3+H3K9me2+REC8+ASY1+CENH3+SPO11+SNV)^2")
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
