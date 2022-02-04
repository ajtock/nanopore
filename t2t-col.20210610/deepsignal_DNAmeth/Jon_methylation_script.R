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


library(reshape2)
library(pheatmap)
library(ggplot2)
library(irr)
library(ggrepel)


################################################################################
# Load in the files----
################################################################################


#table of methylatio
meCtable = read.table("./data/chr1_fast5s.CG.call_mods.tsv")
meCtable = meCtable[order(meCtable$V2) , ]
head(meCtable)

max(meCtable$V2)
nrow(meCtable)



################################################################################
# Get a subset
################################################################################

meCtableHead = meCtable[1:10000000,]
colnames(meCtableHead ) = c("chr", "pos", "strand", "posStrand", "readName", "readStrand",
                            "probMec", "probNot", "MEC", "context")
maxPos = max(meCtableHead$pos)

#subset annotation


################################################################################
# TE Annotation
################################################################################
teAnnotation = read.table("./data/t2t.tes.bed", stringsAsFactors = F)



teStarts = as.data.frame(table(signif(teAnnotation$V2, 3)))
colnames(teStarts) = c("pos", "Freq")
teStarts$bp = 0
for(i in 1:nrow(teStarts)){
  teStartstmp = teAnnotation[signif(teAnnotation$V2, 3) == teStarts$pos[i] ,  ]
  teStarts$bp[i] = sum(teStartstmp$V3 - teStartstmp$V2)
}
head(teStarts)


teStarts = teStarts[as.numeric(as.character(teStarts$pos)) < maxPos ,]
#teStarts = as.data.frame(table(signif(teAnnotation$V2, 3)))
#colnames(teStarts) = c("pos", "Freq")







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
    
    tmpMectable  = meCtableHead[meCtableHead$pos > binStart & meCtableHead$pos < binEnd , ]
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
    
    tmpMectable  = meCtableHead[meCtableHead$pos > binStart & meCtableHead$pos < binEnd , ]
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
maxBinStart = signif(max(meCtableHead$pos),sigFig )

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

binScores = data.frame()




for(i in 1:numberOfBins){
  binStart = (i * windowSize) - (windowSize-1)
  binEnd = binStart + (windowSize-1)
  binStarts[i] = binStart
  tmpMectable  = meCtableHead[meCtableHead$pos > binStart & meCtableHead$pos <= binEnd , ]
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
    
    binnedKappa = vector("numeric", (1/meCbinSize))
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


hist(regionMetrics$wholeBinKappa,breaks = 30)

ggplot(regionMetrics) + 
  geom_point(aes(x = meanBinMec, y = wholeBinKappa))+
  geom_hline(yintercept = 0.6)+
  geom_hline(yintercept = 0.4)+
  theme_classic() +
  geom_text_repel(aes(label = window,x = meanBinMec, y = (wholeBinKappa)),
                  color = "gray20",
                  data = regionMetrics,
                  force = 10)



heatmapWindow(tmpMectable,5000, 88)







################################################################################
# Here one could subset the region by some statistics and then
# find the those region with a "bi-modal" distribution
################################################################################














ggplot() + 
  geom_line(data = regionMetrics,aes(x = binStart, y = log10(wholeBinKappa)))+
  geom_line(data =  teStarts,  aes(x = as.numeric(as.character(pos)),colour = "red", y = -log10(bp         ))) +
  geom_hline(yintercept = log10(0.6))+
  geom_hline(yintercept = log10(0.4))+
  theme_classic()


ggplot(regionMetrics) + 
  geom_point(aes(x = meanBinMec, y = log10(wholeBinKappa)))+
  geom_hline(yintercept = log10(0.6))+
  geom_hline(yintercept = log10(0.4))+
  theme_classic()

hist(regionMetrics$wholeBinKappa,breaks = 30)



which(regionMetrics$meanBinMec < 0.5 & regionMetrics$wholeBinKappa < 10)





