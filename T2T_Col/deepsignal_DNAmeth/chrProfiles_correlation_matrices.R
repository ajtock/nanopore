#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 10.03.2021

# Create and plot correlation matrices of chromosome profiles for
# DNA methylation proportions derived from BS-seq, nanopolish or deepsignal

# Usage:
# ./chrProfiles_correlation_matrices.R T2T_Col 10000 chromwide 'Chr1,Chr2,Chr3,Chr4,Chr5' unsmoothed

#genomeBinSize <- 10000
#refbase <- "T2T_Col"
#region <- "chromwide"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#smoothing <- "unsmoothed"

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
genomeBinSize <- as.numeric(args[2])
region <- args[3]
chrName <- unlist(strsplit(args[4],
                           split = ","))
smoothing <- args[5]

plotDir <- paste0("correlation_matrices/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

options(stringsAsFactors = F)
library(GenomicRanges)
library(Hmisc) # includes rcorr() function which computes significance levels for Pearson and Spearman correlations
library(reshape)
library(ggplot2)
library(ggcorrplot)

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]

# Define region to be analysed
if(region == "nonCEN") {
  regionGR <- GRanges(seqnames = rep(chrs, 2),
                      ranges = IRanges(start = c(rep(1, length(chrs)), CENend+1),
                                       end = c(CENstart-1, chrLens)),
                      strand = "*")
} else if(region == "CEN") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = CENstart,
                                       end = CENend),
                      strand = "*")
} else if(region == "chromwide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
} else {
  stop("region is not nonCEN, CEN, or chromwide")
}

# Define region to be masked out of analysis
if(region == "nonCEN") {
  maskGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = CENstart,
                                     end = CENend),
                    strand = "*")
} else if(region == "CEN") {
  maskGR <- GRanges(seqnames = rep(chrs, 2),
                    ranges = IRanges(start = c(rep(1, length(chrs)), CENend+1),
                                     end = c(CENstart-1, chrLens)),
                    strand = "*")
} else if(region == "chromwide") {
  maskGR <- GRanges()
} else {
  stop("region is not nonCEN, CEN, or chromwide")
}

paths <- c(
           paste0("DNAmeth_WT_nanopolishDNAmeth_95_10kb_MappedOn_", refbase, "_dedup_binSize", genomeBinName, "_", smoothing, ".tsv"),
           paste0("DNAmeth_WT_deepsignalDNAmeth_95_30kb_MappedOn_", refbase, "_dedup_binSize", genomeBinName, "_", smoothing, ".tsv"),
           paste0("DNAmeth_Col_0_BSseq_Rep1_ERR965674_MappedOn_", refbase, "_dedup_binSize", genomeBinName, "_", smoothing, ".tsv"),
           paste0("DNAmeth_Col_0_BSseq_Rep2_ERR965675_MappedOn_", refbase, "_dedup_binSize", genomeBinName, "_", smoothing, ".tsv")
          )

profileNames <- c(
                  "Nanopolish mCG",
                  "DeepSignal mCG",
                  "DeepSignal mCHG",
#                  "DeepSignal mCHH",
                  "BS-seq Rep1 mCG",
                  "BS-seq Rep1 mCHG",
                  "BS-seq Rep1 mCHH",
                  "BS-seq Rep2 mCG",
                  "BS-seq Rep2 mCHG",
                  "BS-seq Rep2 mCHH"
)

profiles <- lapply(seq_along(paths), function(x) {
  read.table(paths[x], header = T)
})
profilesNew <- lapply(seq_along(profiles), function(x) {
  # Make columns consistent across data sets
  if(dim(profiles[[x]])[2] == 4) {
    data.frame(chr = profiles[[x]][,1],
               window = profiles[[x]][,2],
               value = profiles[[x]][,dim(profiles[[x]])[2]],
               stringsAsFactors = F)
  # Separate DNA methylation contexts and make columns consistent as above
  } else if(dim(profiles[[x]])[2] == 6) {
    DNAmethList <- list()
    # Exlcude average over all contexts
    for(y in 4:5) {
      print(y)
      DNAmethList[[y]] <- data.frame(chr = profiles[[x]][,1],
                                     window = profiles[[x]][,2],
                                     value = profiles[[x]][,y],
                                     stringsAsFactors = F)
    }
    # Remove empty ("NULL") list elements
    DNAmethList[-1:-3]
  # Separate DNA methylation contexts and make columns consistent as above
  } else if(dim(profiles[[x]])[2] == 7) {
    DNAmethList <- list()
    # Exlcude average over all contexts
    for(y in 4:6) {
      print(y)
      DNAmethList[[y]] <- data.frame(chr = profiles[[x]][,1],
                                     window = profiles[[x]][,2],
                                     value = profiles[[x]][,y],
                                     stringsAsFactors = F)
    }
    # Remove empty ("NULL") list elements
    DNAmethList[-1:-3]
  }
})

# DNAmethList in each case is one list of 2 or 3 elements within the 4-element list profilesNew
# To make each of these 2 or 3 elements its own element within a 9-element list profilesNew2:
profilesNew2 <- NULL
for(x in seq_along(profilesNew)) {
  if(class(profilesNew[[x]]) == "list") {
    profilesNew2 <- c(profilesNew2, lapply(profilesNew[[x]], function(y) { return(y) }))
  } else {
    profilesNew2 <- c(profilesNew2, list(profilesNew[[x]]))
  }
}

profiles <- profilesNew2 

profilesGR <- lapply(seq_along(profiles), function(x) {
  GRanges(seqnames = profiles[[x]]$chr,
          ranges = IRanges(start = profiles[[x]]$window,
                           end = profiles[[x]]$window+genomeBinSize-1),
          strand = "*",
          value = profiles[[x]]$value)
}) 
# Redefine end coordinate of last window in each chromosome
# to be equal to the length of that chromosome
profilesGR_chrs <- lapply(seq_along(profilesGR), function(x) {
  profilesGRx <- lapply(seq_along(chrs), function(i) {
    profilesGRx_chr <- profilesGR[[x]][seqnames(profilesGR[[x]]) == chrs[i]]
    end(profilesGRx_chr)[length(profilesGRx_chr)] <- chrLens[i]
    profilesGRx_chr
  })
  profilesGRx
})

# Combine GRanges list elements (1 element for each chromosome)
# in one GRanges object for each profile
profilesGR <- lapply(seq_along(profilesGR_chrs), function(x) {
  do.call(c, profilesGR_chrs[[x]])
})

# Subset to include only those windows not overlapping masked region (e.g., "nonCEN")
profilesGR <- lapply(seq_along(profilesGR), function(x) {
  mask_profilesGRx_overlap <- findOverlaps(query = maskGR,
                                           subject = profilesGR[[x]],
                                           type = "any",
                                           select = "all",
                                           ignore.strand = TRUE)
  if(length(mask_profilesGRx_overlap) != 0) {
    profilesGR[[x]][-subjectHits(mask_profilesGRx_overlap)]
  } else {
    profilesGR[[x]]
  }
})

# Combine profiles into one data.frame in which each profile is a column
profilesVal <- lapply(seq_along(profilesGR), function(x) {
  profilesGR[[x]]$value
})

profilesDF <- as.data.frame(do.call(cbind, profilesVal),
                            stringsAsFactors = F)
colnames(profilesDF) <- profileNames

# Create correlation matrix
corMat <- round(cor(profilesDF,
                    method = "spearman",
                    use = "pairwise.complete.obs"),
                digits = 2)

# Set duplicates to NA
for(x in 1:dim(corMat)[1]) {
  corMat[x, x] <- NA
  if(x > 1) {
    corMat[x, 1:x-1] <- NA
  }
}
corMat <- corMat[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDat <- melt(corMat)
corDat <- corDat[-which(is.na(corDat[,3])),]

# Order the data.frame for plotting
profileNamesList <- as.list(profileNames)
names(profileNamesList) <- profileNames
levels(corDat$X1) <- rev(profileNamesList) 
levels(corDat$X2) <- profileNamesList[-1]

# Get P-values for correlation matrix
corMatSig <- rcorr(as.matrix(profilesDF),
                   type = "spearman")$P
# Set duplicates to NA
for(x in 1:dim(corMatSig)[1]) {
  corMatSig[x, x] <- NA
  if(x > 1) {
    corMatSig[x, 1:x-1] <- NA
  }
}
corMatSig <- corMatSig[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDatSig <- melt(corMatSig)
corDatSig <- corDatSig[-which(is.na(corDatSig[,3])),]

# Standardise P-values to a sample size of 100 (q-values) as proposed by
# Good (1982) Standardized tail-area probabilities. Journal of Computation and Simulation 16: 65-66
# and summarised by Woolley (2003):
# https://stats.stackexchange.com/questions/22233/how-to-choose-significance-level-for-a-large-data-set
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.518.5341&rep=rep1&type=pdf
# Woolley (2003): "Clearly, the meaningfulness of the p-value diminishes as the sample size increases";
# Anne Z. (2012, Pearson eCollege, Denver): "In the real world, there are unlikely to be semi-partial correlations
# that are exactly zero, which is the null hypothesis in testing significance of a regression coefficient."
# Formally, the standardised p-value is defined as:
# q = min(0.5, p * sqrt( (n/100) ))
# Woolley (2003): "The value of 0.5 is somewhat arbitrary, though its purpose is to avoid q-values of greater than 1."
n <- dim(profilesDF)[1]
corDatSig$value <- sapply(corDatSig$value, function(x) {
  round(min(0.5, x * sqrt( (n/100) )),
        digits = 2)
})

# Order the data.frame for plotting
levels(corDatSig$X1) <- rev(profileNamesList) 
levels(corDatSig$X2) <- profileNamesList[-1]

# Plot
ggObj <- ggplot(data = corDat,
                mapping = aes(X2, X1, fill = value)) +
  geom_tile() +
#  geom_text(mapping = aes(X2, X1, label = value), size = 5) +
  geom_text(data = corDatSig,
            mapping = aes(X2, X1, label = value), size = 8) +
  scale_fill_gradient2(name = bquote("Spearman's" ~ italic(r[s])),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, breaks = seq(-1, 1, by = 0.4), limits = c(-1, 1)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 6,
                                title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0, size = 39, colour = "black"),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 39, colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.justification = c(1, 0),
        legend.position = c(0.65, 0.05),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(5.5, 200.5, 5.5, 5.5), "pt"),
        plot.title = element_text(hjust = 0.5, size = 30, colour = "black")) +
  ggtitle(bquote(.(genomeBinSize/1e3) * "-kb Spearman's" ~ italic(r[s]) ~ "for" ~
          .(paste0(chrName, collapse = ",")) ~
          .(region) ~ "regions (" * .(smoothing) * ")"))
ggsave(paste0(plotDir,
              "Spearman_correlation_matrix_", refbase, "_", genomeBinName,
              "_DNAmeth_Nanopolish_DeepSignal_Rigal2016PNASpeBSseq_",
              paste0(chrName, collapse = "_"), "_", region, "_", smoothing, "_qVals.pdf"),
       plot = ggObj, height = 20, width = 22)
