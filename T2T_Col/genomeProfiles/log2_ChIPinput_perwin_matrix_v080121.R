#!/applications/R/R-3.5.0/bin/Rscript

######################################
######################################
# Generate matrix of chromosome profiles
######################################
######################################

# Rscript ./log2_ChIPinput_perwin_matrix_v080121.R 10kb T2T_Col both 'WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input,WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me3_ChIP14_WT_REC8_Myc_Rep1_input,WT_H3K4me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K27me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,H2AW6_ChIP_SRR5298545_H2AW_input_SRR5298544,H2AW7_ChIP_SRR5298546_H2AW_input_SRR5298544,WT_MNase_Rep1_WT_gDNA_Rep1,WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep2_ChIP_WT_REC8_Myc_Rep1_input,kss_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,cmt3_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,ColLerF1_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,WT_SPO11oligos_Rep2_WT_gDNA_Rep1_R1,WT_SPO11oligos_Rep3_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep2_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep3_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep2_WT_gDNA_Rep1_R1'

#genomeBinName <- "10kb"
#refbase <- "T2T_Col"
#align <- "both"
#sampleNames <- unlist(strsplit("WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input,WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me3_ChIP14_WT_REC8_Myc_Rep1_input,WT_H3K4me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K4me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_H3K27me1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,H2AW6_ChIP_SRR5298545_H2AW_input_SRR5298544,H2AW7_ChIP_SRR5298546_H2AW_input_SRR5298544,WT_MNase_Rep1_WT_gDNA_Rep1,WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_MTOPVIB_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_DMC1_V5_Rep2_ChIP_WT_REC8_Myc_Rep1_input,kss_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,cmt3_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,ColLerF1_DMC1_V5_Rep1_ChIP_WT_REC8_Myc_Rep1_input,WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,WT_SPO11oligos_Rep2_WT_gDNA_Rep1_R1,WT_SPO11oligos_Rep3_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep2_WT_gDNA_Rep1_R1,met1_SPO11oligos_Rep3_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,kss_SPO11oligos_Rep2_WT_gDNA_Rep1_R1",
#                              split = ","))

args <- commandArgs(trailingOnly = T)
genomeBinName <- as.character(args[1])
refbase <- args[2]
align <- args[3]
sampleNames <- unlist(strsplit(args[4],
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

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

chrWindows <- read.table(paste0(sampleNames[1],
                                "_MappedOn_", refbase, "_lowXM_", align,
                                "_sort_norm_binSize",
                                 genomeBinName, "_unsmoothed.tsv"),
                         header = T, colClasses = c(rep(NA, 3), rep("NULL", 3)))

sample_winDF_list <- mclapply(seq_along(sampleNames), function(x) {
  read.table(paste0(sampleNames[x],
                    "_MappedOn_", refbase, "_lowXM_", align,
                    "_sort_norm_binSize",
                     genomeBinName, "_unsmoothed.tsv"),
             header = T, colClasses = c(rep("NULL", 5), rep(NA, 1)))
}, mc.cores = length(sampleNames))

sample_filt_winDF_list <- mclapply(seq_along(sampleNames), function(x) {
  read.table(paste0(sampleNames[x],
                    "_MappedOn_", refbase, "_lowXM_", align,
                    "_sort_norm_binSize",
                     genomeBinName, "_smoothed.tsv"),
             header = T, colClasses = c(rep("NULL", 5), rep(NA, 1)))
}, mc.cores = length(sampleNames))

sample_winDF_matrix <- cbind(chrWindows,
                             do.call(cbind, sample_winDF_list))
colnames(sample_winDF_matrix) <- c("chr", "window", "cumwindow",
                                   paste0("log2_", sampleNames))
sample_filt_winDF_matrix <- cbind(chrWindows,
                                  do.call(cbind, sample_filt_winDF_list))
colnames(sample_filt_winDF_matrix) <- c("chr", "window", "cumwindow",
                                        paste0("log2_", sampleNames))

write.table(sample_winDF_matrix,
            file = paste0("log2ChIPcontrol_MappedOn_", refbase, "_",
                          "genome_norm_coverage_matrix_binSize",
                          genomeBinName, "_unsmoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(sample_filt_winDF_matrix,
            file = paste0("log2ChIPcontrol_MappedOn_", refbase, "_",
                          "genome_norm_coverage_matrix_binSize",
                          genomeBinName, "_smoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
