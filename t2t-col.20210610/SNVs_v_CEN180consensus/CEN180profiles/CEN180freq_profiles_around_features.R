#!/applications/R/R-4.0.0/bin/Rscript

# Profile CEN180 frequency around a given feature set

# Usage:
# /applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 180 2000 2kb 10 10bp genomewide CEN180
# /applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENgap
# /applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENAthila
# /applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENsoloLTR

#chrName <- "Chr1"
#regionBodyLength <- 180
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
#interval <- "genomewide"
#featureName <- "CEN180"

args <- commandArgs(trailingOnly = T)
chrName <- args[1]
regionBodyLength <- as.numeric(args[2])
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
binSize <- as.numeric(args[5])
binName <- args[6]
interval <- args[7]
featureName <- args[8]

options(stringsAsFactors = F)
library(EnrichedHeatmap)
library(parallel)

matDir <- paste0("matrices_unsmoothed/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14841110,3823792,13597188,4203902,11784131)
CENend <- c(17559778,6045243,15733925,6977949,14551809)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")
genomeGR <- genomeGR[grep(chrName,
                          seqnames(genomeGR))@values]
CENGR <- CENGR[grep(chrName,
                    seqnames(CENGR))@values]

# Define interval to be analysed
if(interval == "centromeres") {
  intervalGR <- GRanges(seqnames = chrs,
                        ranges = IRanges(start = CENstart,
                                         end = CENend),
                        strand = "*")
  intervalGR <- intervalGR[grep(chrName,
                                seqnames(intervalGR))@values]
} else if(interval == "genomewide") {
  intervalGR <- GRanges(seqnames = chrs,
                        ranges = IRanges(start = 1,
                                         end = chrLens),
                        strand = "*")
  intervalGR <- intervalGR[grep(chrName,
                                seqnames(intervalGR))@values]
} else {
  stop("interval is not centromeres or genomewide")
}

# Define interval to be masked out of analysis
if(interval == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrs, 2),
                    ranges = IRanges(start = c(rep(1, length(chrs)),
                                               CENend+1),
                                     end = c(CENstart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(chrName,
                        seqnames(maskGR))@values]
} else if(interval == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(chrName,
                        seqnames(maskGR))@values]
} else {
  stop("interval is not centromeres or genomewide")
}

# Load features in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
features <- read.table(paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/", featureName, "/",
                              featureName, "_in_t2t-col.20210610_", chrName, ".bed"),
                       header = F)
colnames(features) <- c("chr", "start", "end", "name", "score", "strand")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start+1,
                                       end = features$end),
                      strand = features$strand,
                      number = features$name)
print(length(featuresGR))
# Subset to include only those not overlapping masked interval
mask_features_overlap <- findOverlaps(query = maskGR,
                                      subject = featuresGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
if(length(mask_features_overlap) > 0) {
  featuresGR <- featuresGR[-subjectHits(mask_features_overlap)]
}
print(length(featuresGR))

# Load ranLoc in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
if(featureName %in% c("CEN180", "CENgapAll")) {
  ranLoc <- read.table(paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/", featureName, "/",
                              featureName, "_in_t2t-col.20210610_", chrName, "_CENrandomLoci.bed"),
                       header = F)
  colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
  ranLocGR <- GRanges(seqnames = ranLoc$chr,
                      ranges = IRanges(start = ranLoc$start+1,
                                       end = ranLoc$end),
                      strand = ranLoc$strand,
                      number = ranLoc$name)
  print(length(ranLocGR))
#  # Subset to include only those not overlapping masked interval
#  mask_ranLoc_overlap <- findOverlaps(query = maskGR,
#                                      subject = ranLocGR,
#                                      type = "any",
#                                      select = "all",
#                                      ignore.strand = TRUE)
#  if(length(mask_ranLoc_overlap) > 0) {
#    ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
#  }
#  print(length(ranLocGR))
}

# Load genomic physical coordinates of CEN180 sequences
# Note addition of 1 to 0-based BED start coordinates
CEN180 <- read.table(paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/CEN180/",
                            "CEN180_in_t2t-col.20210610_", chrName, ".bed"),
                       header = F)
colnames(CEN180) <- c("chr", "start", "end", "name", "score", "strand", "HORlengthsSum", "HORcount", "percentageIdentity")
CEN180GR <- GRanges(seqnames = CEN180$chr,
                    ranges = IRanges(start = CEN180$start+1,
                                     end = CEN180$end),
                    strand = "*",
                    coverage = rep(1, dim(CEN180)[1]))

if(featureName %in% c("CEN180", "CENgapAll")) {
  # Define matrix and column mean outfiles
  outDF <- list(paste0(matDir,
                       "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                       "_matrix_bin", binName, "_flank", flankName, ".tab"),
                paste0(matDir,
                       "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                       "_CENranLoc_matrix_bin", binName, "_flank", flankName, ".tab"))
  outDFcolMeans <- list(paste0(matDir,
                               "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                               "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
                        paste0(matDir,
                               "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                               "_CENranLoc_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
  
  # Function to create CEN180 frequency matrices for
  # feature loci and random loci (incl. flanking regions)
  # and to calculate mean profiles across all feature loci and random loci
  covMatrix <- function(signal,
                        feature,
                        ranLoc,
                        featureSize,
                        flankSize,
                        winSize,
                        outDF,
                        outDFcolMeans) {
    # feature loci
    set.seed(2840)
    feature_smoothed <- normalizeToMatrix(signal = signal,
                                          target = feature,
                                          value_column = "coverage",
                                          extend = flankSize,
                                          mean_mode = "w0",
                                          w = winSize,
                                          background = 0,
                                          smooth = FALSE,
                                          include_target = TRUE,
                                          target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("feature_smoothed")
    print(feature_smoothed)
    print("feature_smoothed rows = ")
    print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    feature_smoothed_DF <- data.frame(feature_smoothed)
    feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                       na.rm = T))
    write.table(feature_smoothed_DF,
                file = outDF[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(feature_smoothed_DF_colMeans,
                file = outDFcolMeans[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  
    # random loci
    set.seed(8472)
    ranLoc_smoothed <- normalizeToMatrix(signal = signal,
                                         target = ranLoc,
                                         value_column = "coverage",
                                         extend = flankSize,
                                         mean_mode = "w0",
                                         w = winSize,
                                         background = 0,
                                         smooth = FALSE,
                                         include_target = TRUE,
                                         target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("ranLoc_smoothed")
    print(ranLoc_smoothed)
    print("ranLoc_smoothed rows = ")
    print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    ranLoc_smoothed_DF <- data.frame(ranLoc_smoothed)
    ranLoc_smoothed_DF_colMeans <- as.vector(colMeans(ranLoc_smoothed_DF,
                                                      na.rm = T))
    write.table(ranLoc_smoothed_DF,
                file = outDF[[2]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(ranLoc_smoothed_DF_colMeans,
                file = outDFcolMeans[[2]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
  
  # Run covMatrix() function on each feature GRanges object to obtain matrices
  # containing normalised feature density values around target and random loci
  covMatrix(signal = CEN180GR,
            feature = featuresGR,
            ranLoc = ranLocGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            winSize = binSize,
            outDF = outDF,
            outDFcolMeans = outDFcolMeans)
  print(paste0(
               "CEN180 frequency around ", featureName, " in ", chrName,
               " profile calculation complete"))
} else {
  # Define matrix and column mean outfiles
  outDF <- list(paste0(matDir,
                       "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                       "_matrix_bin", binName, "_flank", flankName, ".tab"))
  outDFcolMeans <- list(paste0(matDir,
                               "CEN180freq_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                               "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
  
  # Function to create CEN180 frequency matrices for
  # feature loci and random loci (incl. flanking regions)
  # and to calculate mean profiles across all feature loci and random loci
  covMatrix <- function(signal,
                        feature,
                        featureSize,
                        flankSize,
                        winSize,
                        outDF,
                        outDFcolMeans) {
    # feature loci
    set.seed(2840)
    feature_smoothed <- normalizeToMatrix(signal = signal,
                                          target = feature,
                                          value_column = "coverage",
                                          extend = flankSize,
                                          mean_mode = "w0",
                                          w = winSize,
                                          background = 0,
                                          smooth = FALSE,
                                          include_target = TRUE,
                                          target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("feature_smoothed")
    print(feature_smoothed)
    print("feature_smoothed rows = ")
    print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    feature_smoothed_DF <- data.frame(feature_smoothed)
    feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                       na.rm = T))
    write.table(feature_smoothed_DF,
                file = outDF[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(feature_smoothed_DF_colMeans,
                file = outDFcolMeans[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
  
  # Run covMatrix() function on each feature GRanges object to obtain matrices
  # containing normalised feature density values around target and random loci
  covMatrix(signal = CEN180GR,
            feature = featuresGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            winSize = binSize,
            outDF = outDF,
            outDFcolMeans = outDFcolMeans)
  print(paste0(
               "CEN180 frequency around ", featureName, " in ", chrName,
               " profile calculation complete"))
}
