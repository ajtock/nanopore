################################################################################
# The idea is to 
# analyse the stochasticity of methylation in nanopore data
# this should give a good idea of regions in the genomes that 
# are succeptible to mutibility during an organisms life 
# therefore cause ratios of methylation != [0|.5|1]
# and therefore suggesting mosaicism within the tissue
################################################################################

################################################################################
# Load in the libraries ----
################################################################################

#sampleName <- "Col_0_deepsignalDNAmeth_30kb_90pc"
#refbase <- "t2t-col.20210610"
#genomeBinSize <- 10000
#genomeStepSize <- 10000
#context <- "CpG"
#NAmax <- 0.50
#CPUpc <- 1.00
#chrName <- unlist(strsplit("Chr1", split = ","))

library(reshape2)
library(pheatmap)
library(ggplot2)
library(irr)
library(ggrepel)

print(paste0("Proportion of CPUs:", CPUpc))
options(stringsAsFactors = F)
library(parallel)
library(dplyr)
library(tidyr)


################################################################################
# Load in the files----
################################################################################


#table of methylatio
# Read in the raw output .tsv file from Deepsignal methylation model
tab_list <- mclapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/deepsignal_DNAmeth/",
                    sampleName, "_MappedOn_", refbase, "_", context, "_raw_", chrName[x], ".tsv"),
             header = F)
}, mc.cores = length(chrName), mc.preschedule = F)

if(length(chrName) > 1) {
  tab <- dplyr::bind_rows(tab_list)
} else {
  tab <- tab_list[[1]]
}
rm(tab_list); gc()


meCtable <- tab[ with(tab, order(V1, V2)), ]
rm(tab); gc()
head(meCtable)
max(meCtable$V2)
nrow(meCtable)



################################################################################
# Get a subset
################################################################################

meCtableHead = meCtable[1:10000000,]
colnames(meCtableHead ) = c("chr", "pos", "strand", "posStrand", "readName", "readStrand",
                            "prob_0", "prob_1", "MEC", "k_mer")
maxPos = max(meCtableHead$pos)

#subset annotation


################################################################################
# TE Annotation
################################################################################
teAnnotation = read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                                 "/annotation/TEs_EDTA/", refbase, "_TEs_All_",
                                 paste0(chrName, collapse = "_"), ".bed"),
                          stringsAsFactors = F)



teStarts = as.data.frame(table(signif(teAnnotation$V2, 3)))
colnames(teStarts) = c("pos", "Freq")
teStarts$bp = 0
for(i in 1:nrow(teStarts)){
  teStartstmp = teAnnotation[signif(teAnnotation$V2, 3) == teStarts$pos[i] ,  ]
  teStarts$bp[i] = sum(teStartstmp$V3 - teStartstmp$V2)
}
head(teStarts)

teStarts = teStarts[as.numeric(as.character(teStarts$pos)) < maxPos ,]

################################################################################
# Subset just for positive strand ( simplify the data for ease )
################################################################################
meCtableHead = meCtableHead[meCtableHead$strand == "+",]



################################################################################
# Raw coverage trying to find the right bin size
################################################################################

coverageDF = data.frame()
for(window in c(100,500,1000,2000,3000,4000,5000,6000,10000,20000,25000,30000)){
  
  windowSize = window
  maxPos = max(meCtableHead$pos)
  sigFig =  nchar(maxPos) - nchar(windowSize) + 1
  maxBinStart = signif(max(meCtableHead$pos),sigFig )
  
  numberOfBins = maxBinStart / windowSize
  
  
  cov = vector("numeric", numberOfBins)
  for(i in 1:numberOfBins){
    binStart = (i * windowSize) - (windowSize-1)
    binEnd = binStart + (windowSize-1)
    
    tmpMectable  = meCtableHead[meCtableHead$pos > binStart &
                                meCtableHead$pos < binEnd , ]
    cov[i] = length(unique(tmpMectable$readName))
    
  }
  
  tmpDf = data.frame("windowSize" = windowSize,
                     "coverage" = cov)
  
  
  coverageDF = rbind.data.frame(coverageDF,tmpDf, stringsAsFactors = F)
}


ggplot(coverageDF) + 
  geom_violin(aes(x = as.factor(windowSize), y = coverage))+
  theme_classic()

################################################################################
# full Read coverage trying to find the right bin size
################################################################################



coverageDF = data.frame()
for(window in c(100,500,1000,2000,3000,4000,5000,6000,10000,20000,25000,30000)){
  
  windowSize = window
  maxPos = max(meCtableHead$pos)
  sigFig =  nchar(maxPos) - nchar(windowSize) + 1
  maxBinStart = signif(max(meCtableHead$pos),sigFig )
  
  numberOfBins = maxBinStart / windowSize
  
  
  cov = vector("numeric", numberOfBins)
  for(i in 1:numberOfBins){
    binStart = (i * windowSize) - (windowSize-1)
    binEnd = binStart + (windowSize-1)
    
    tmpMectable  = meCtableHead[meCtableHead$pos > binStart &
                                meCtableHead$pos < binEnd, ]
    if(nrow(tmpMectable) == 0){
      cov[i] = 0
    } else {
    tmpMatrix = dcast(tmpMectable, pos ~ readName, value.var = "MEC")
    
    
    tmpMatrix = tmpMatrix[,(complete.cases(t(tmpMatrix)))]
    
    cov[i] = ncol(tmpMatrix )
    }
    
  }
  
  tmpDf = data.frame("windowSize" = windowSize,
                     "coverage" = cov)
  
  
  coverageDF = rbind.data.frame(coverageDF,tmpDf, stringsAsFactors = F)
}


ggplot(coverageDF) + 
  geom_violin(aes(x = as.factor(windowSize), y = coverage))+
  theme_classic()





################################################################################
# Now get a score and other metrics for each bin
################################################################################

# with a window size of 5kb check the fliess kappa of each window
# to find the good ones then plot them. Also grab other metrics for 
# each region whilst we are there

windowSize = 5000
maxPos = max(meCtableHead$pos)
sigFig =  nchar(maxPos) - nchar(windowSize) + 1
maxBinStart = signif(maxPos, sigFig)

numberOfBins = maxBinStart / windowSize


binCov = vector("numeric", numberOfBins)
binCyt = vector("numeric", numberOfBins)
meanBinMeC = vector("numeric", numberOfBins)
cov = vector("numeric", numberOfBins)
notMec = vector("numeric", numberOfBins)
wholeBinKappa =  vector("numeric", numberOfBins)
binStarts = vector("numeric", numberOfBins)
meCbinSize = 0.3
meCBins = seq(meCbinSize,1,meCbinSize)
binnedKappa = vector("numeric", (1/meCbinSize))

binScores = data.frame()




for(i in 1:numberOfBins){
  binStart = (i * windowSize) - (windowSize-1)
  binEnd = binStart + (windowSize-1)
  binStarts[i] = binStart
  # DeepSignal raw and methylation frequency file cytosine coordinates are 0-based
  tmpMectable  = meCtableHead[meCtableHead$pos >= binStart-1 &
                              meCtableHead$pos <= binEnd-1, ]
  if(nrow(tmpMectable) == 0){
    cov[i] = 0
    binScores = rbind.data.frame(binScores,binnedKappa )
  } else {
    
    ##Body
    
    #get a matrix of the region
    tmpMatrix = dcast(tmpMectable, pos ~ readName, value.var = "MEC")
    tmpMatrix = tmpMatrix[,(complete.cases(t(tmpMatrix)))]
    row.names(tmpMatrix) = tmpMatrix$pos
    tmpMatrix = tmpMatrix[,-1]
    
    ratio = rowSums(tmpMatrix) / ncol(tmpMatrix)
    
    # get the number of reads and cytosines for this region
    binCov[i] = ncol(tmpMatrix )
    binCyt[i] = nrow(tmpMatrix)
    
    # FleissKappa for whole region
    fleissKappa = kappam.fleiss(tmpMatrix)
    wholeBinKappa[i] = fleissKappa$value
    
    
    #mean bin methylation

    meanBinMeC[i] = mean(as.matrix(tmpMatrix))
    
    #### scores
    
    #number of non-methylated C's
    notMec[i] = length(which( rowSums(tmpMatrix) == 0 ))
    
    #scores in binned methylation values
    for(bin in 1:length(meCBins) ){
      higherMecScore = meCBins[bin]
      lowerMecScore = higherMecScore - meCbinSize
      binCytosines = which(ratio > lowerMecScore & ratio <= higherMecScore)
      if( length(binCytosines) ==0 ){       binnedKappa[bin] = NA
      } else {
      binMat = tmpMatrix[binCytosines,]
      fleissKappa = kappam.fleiss(binMat)
      testStatistic = fleissKappa$value
      
      binnedKappa[bin] = testStatistic
      }
    }
    binnedKappa = c(binnedKappa, i)

    binScores = rbind.data.frame(binScores,binnedKappa )
    colnames(binScores) = c(paste("bin_",1:(1/meCbinSize)), "window")
    
    
    
    ##
    
  }
  
}


# make the dataFrame

regionMetrics = data.frame("binCoverage" = binCov,
                           "numberOfCytosines" = binCyt,
                           "numberNoMethylation" = notMec,
                           "wholeBinKappa" = wholeBinKappa,
                           "meanBinMec" =  meanBinMeC,
                           "binStart" = binStarts)

regionMetrics = cbind.data.frame(regionMetrics, binScores)


regionMetrics = regionMetrics[regionMetrics$numberOfCytosines > 10,]




################################################################################
# identify exciting regions and plot them
################################################################################

#kappaScores = regionMetrics[,5:((5+(1/meCbinSize))-1)]
#kappaScores = melt(kappaScores)




which(regionMetrics$meanBinMec < 0.5 & regionMetrics$wholeBinKappa < 10)

plotDir <- "per_read_analysis/among_read_variation/"

pdf(paste0(plotDir, "wholeBinKappa_hist.pdf"))
hist(regionMetrics$wholeBinKappa,breaks = 30)
dev.off()

regionMetrics <- data.frame(regionMetrics, window = 1:nrow(regionMetrics))

meanBinMec_vs_wholeBinKappa_gg <- ggplot(regionMetrics) + 
  geom_point(aes(x = meanBinMec, y = wholeBinKappa))+
  geom_hline(yintercept = 0.6)+
  geom_hline(yintercept = 0.4)+
  theme_classic()
#  geom_text_repel(aes(label = window,x = meanBinMec, y = (wholeBinKappa)),
#                  color = "gray20",
#                  data = regionMetrics,
#                  force = 10)
ggsave(paste0(plotDir, "meanBinMec_vs_wholeBinKappa.pdf"),
       plot = meanBinMec_vs_wholeBinKappa_gg)

# heatmapWindow() function not defined
#heatmapWindow(tmpMectable,5000, 88)







################################################################################
# Here one could subset the region by some statistics and then
# find the those region with a "bi-modal" distribution
################################################################################














binStart_vs_wholeBinKappa_gg <- ggplot() + 
  geom_line(data = regionMetrics,aes(x = binStart, y = log10(wholeBinKappa)))+
  geom_line(data =  teStarts,  aes(x = as.numeric(as.character(pos)),colour = "red", y = -log10(bp         ))) +
  geom_hline(yintercept = log10(0.6))+
  geom_hline(yintercept = log10(0.4))+
  theme_classic()
ggsave(paste0(plotDir, "binStart_vs_wholeBinKappa.pdf"),
       plot = binStart_vs_wholeBinKappa_gg)

meanBinMec_vs_log10_wholeBinKappa_gg <- ggplot(regionMetrics) + 
  geom_point(aes(x = meanBinMec, y = log10(wholeBinKappa)))+
  geom_hline(yintercept = log10(0.6))+
  geom_hline(yintercept = log10(0.4))+
  theme_classic()
ggsave(paste0(plotDir, "meanBinMec_vs_log10_wholeBinKappa.pdf"),
       plot = meanBinMec_vs_log10_wholeBinKappa_gg)

