#!/usr/bin/env python

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 19/07/22

# Build supervised learning model to predict lethal-phenotype genes,
# with features (predictors) including among-read DNA methylation variation (Fleiss' kappa or Krippendorff's alpha),
# within-read (site-to-site) DNA methylation variation (stochasticity), and
# others from 


# ==== Import libraries ====
import argparse
import pandas as pd
import sys 
import os

from sklearn.neighbors import KNeighborsClassifier


# ==== Capture user input as command-line arguments ====
def parse_args(args):
    parser = argparse.ArgumentParser(description = 'DeepSignal-derived filename variables.')
    #### Define command-line arguments
    parser.add_argument('-s', '--sampleName', type=str, default='Col_0_deepsignalDNAmeth_30kb_90pc',
                        help='The sample prefix of the DeepSignal-derived file containing DNA methylation variation values. Default: Col_0_deepsignalDNAmeth_30kb_90pc')
    parser.add_argument('-r', '--refbase', type=str, default='t2t-col.20210610',
                        help='The base name (excluding ".fa") of the reference assembly to which data were mapped. Default: t2t-col.20210610')
    parser.add_argument('-c', '--context', type=str, default='CpG',
                        help='The DNA methylation sequence context to be analysed. Default: CpG')
    parser.add_argument('-n', '--NAmax', type=float, default=0.50,
                        help='The maximum number of missing mC values across all sites in a read (set by prior filtering and included in filename). Default: 0.50')
    parser.add_argument('-cN', '--chrName', type=str, default='Chr1 Chr2 Chr3 Chr4 Chr5'.split(),
                        help='The chromosomes to be analysed. Default: Chr1 Chr2 Chr3 Chr4 Chr5')
#    parser.add_argument('-lT', '--locusType', type=str, default='gene',
#                        help='The genomic locus type to be analysed. Default: gene')
#    parser.add_argument('-lR', '--locusRegion', type=str, default='regions',
#                        help='The region type of the genomic loci to be analysed. Default: regions (i.e., 1 kb upstream of TSS to 1 kb downstream of TTS)')
    #### Create parser
    return parser.parse_args(args)
 
parser = parse_args(sys.argv[1:])
#parser = parse_args([])

class Parser_test:
    def __init__(self):
        print("in init")
    def test_parser(self):
        parser = parse_args(['-s Col_0_deepsignalDNAmeth_30kb_90pc', '-r t2t-col.20210610', '-c CpG', '-n 0.50', '-cN Chr1 Chr2 Chr3 Chr4 Chr5'])
        self.assertTrue(parser.long)

t = Parser_test()        
t.test_parser()

outDir = str(parser.locusType) + "_" + str(parser.locusRegion) + "/" + str('_'.join(parser.chrName)) + "/"
plotDir = str(outDir + "plots/")

try:
    os.mkdir(outDir)
except OSError as error:
    print(error)

try:
    os.mkdir(plotDir)
except OSError as error:
    print(error)

def load_target_DF():



def main():



outDir <- paste0(featName, "_", featRegion, "/", paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))


# Read in Lloyd et al. (2015) Plant Cell lethal-phenotype status for each gene (target, a.k.a. response variable)
ds1_DF = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/TPC2015-00051-LSBR3_Supplemental_Data_set_1_Sheet1.csv")
ds1_DF.columns = ["gene", "phenotype", "predicted_lethal"]

# Read in Lloyd et al. (2015) Plant Cell gene features, a.k.a predictor variables
ds3_DF = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv")
ds3_DF.rename(columns = {"Locus number":"gene"}, inplace = True)

# Read in gene DNA methylation features 
mC_DF <- pd.read_csv("

# Load among-read and within-read mC data for featName featRegion
featDF <- read.table(paste0(outDir,
                            featName, "_", featRegion, "_", sampleName, "_MappedOn_", refbase,
                            "_", context,
                            "_NAmax", NAmax,
                            "_filt_df_fk_kappa_all_mean_mC_all_complete_",
                            paste0(chrName, collapse = "_"), ".tsv"),
                     header = T)
colnames(featDF)[which(colnames(featDF) == "fk_kappa_all")] <- "kappa"
colnames(featDF)[which(colnames(featDF) == "ka_alpha_all")] <- "alpha"
colnames(featDF)[which(colnames(featDF) == "mean_stocha_all")] <- "stocha"
colnames(featDF)[which(colnames(featDF) == "mean_min_acf_all")] <- "min_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mean_acf_all")] <- "mean_ACF"
colnames(featDF)[which(colnames(featDF) == "mean_mC_all")] <- paste0("mean_m", context)
featDF$kappa_C_density <- featDF$fk_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$stocha_C_density <- featDF$stocha_Cs_all / ( (featDF$end - featDF$start + 1) / 1e3)
featDF$parent <- sub(pattern = "\\.\\d+", replacement = "", x = featDF$name)
featDF$parent <- sub(pattern = "_\\d+", replacement = "", x = featDF$parent)


df = pd.merge(ds1_DF, ds3_DF, how = "inner", on = "gene")
df.replace(to_replace = "?", value = NA, inplace = True)
df = df[df["phenotype"].isin(["Lethal", "Non-Lethal"])]
df.phenotype.replace(["Non-Lethal", "Lethal"], [0, 1], inplace = True)
# Drop rows containing missing values across ANY of the features (predictors) ??
df.dropna(axis = 0, inplace = True)

df_features = df.loc[:, ~df.columns.isin(["gene", "phenotype", "predicted_lethal"])]
df_target = df["phenotype"]

X = df_features.values
y = df_target.values

#print(X.shape, y.shape)
## With NaNs
## (3443, 57) (3443,)

print(X.shape, y.shape)
# Without NaNs
# (3443, 57) (3443,)


knn = KNeighborsClassifier(n_neighbors = 15)
knn.fit(X, y)

