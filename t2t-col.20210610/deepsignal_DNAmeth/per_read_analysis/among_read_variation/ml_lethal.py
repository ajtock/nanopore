#!/usr/bin/env python

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 19/07/22

# Build supervised learning model to predict lethal-phenotype genes,
# with features (predictors) including among-read DNA methylation variation (Fleiss' kappa or Krippendorff's alpha),
# within-read (site-to-site) DNA methylation variation (stochasticity), and
# others from 


# ==== Import libraries
import argparse
import unittest
import sys 
import os
import re
import pandas as pd
import numpy as np

from time import time, sleep
from functools import reduce
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.compose import make_column_transformer
from sklearn.compose import make_column_selector
from sklearn.preprocessing import OrdinalEncoder
from sklearn.datasets import load_diabetes
from sklearn.datasets import fetch_california_housing
from sklearn.datasets import fetch_openml
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LinearRegression


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="DeepSignal-derived filename variables.")
    #### Define command-line arguments
    parser.add_argument("-s", "--sampleName", type=str, default="Col_0_deepsignalDNAmeth_30kb_90pc",
                        help="The sample prefix of the DeepSignal-derived file containing DNA methylation variation values. Default: Col_0_deepsignalDNAmeth_30kb_90pc")
    parser.add_argument("-r", "--refbase", type=str, default="t2t-col.20210610",
                        help="The base name (excluding '.fa') of the reference assembly to which data were mapped. Default: t2t-col.20210610")
    parser.add_argument("-c", "--context", type=str, default="CpG",
                        help="The DNA methylation sequence context to be analysed. Default: CpG")
    parser.add_argument("-n", "--NAmax", type=float, default=0.50,
                        help="The maximum number of missing mC values across all sites in a read (set by prior filtering and included in filename). Default: 0.50")
    parser.add_argument("-cN", "--chrName", type=str, default="Chr1 Chr2 Chr3 Chr4 Chr5".split(),
                        help="The chromosomes to be analysed. Default: Chr1 Chr2 Chr3 Chr4 Chr5")
    parser.add_argument("-lT", "--locusType", type=str, default="gene",
                        help="The genomic locus type to be analysed. Default: gene")
    parser.add_argument("-lR", "--locusRegion", type=str, default="regions",
                        help="The region type of the genomic loci to be analysed. Default: regions (i.e., 1 kb upstream of TSS to 1 kb downstream of TTS)")
    #### Create parser
    return parser
 
parser = create_parser().parse_args()
#parser = parse_args([])

#class TestParser(unittest.TestCase):
#    def setUp(self):
#        self.parser = create_parser()
#    def parser_test(self):
#        parsed = self.parser.parse_args(["sampleName", "Col_0_deepsignalDNAmeth_30kb_90pc"])
#        self.assertEqual(parsed.sampleName, "Col_0_deepsignalDNAmeth_30kb_90pc")
#
#if __name__ == "__main__":
#    unittest.main()

#def test_parser(self):
#    parser = parse_args(["-s Col_0_deepsignalDNAmeth_30kb_90pc", "-r t2t-col.20210610", "-c CpG", "-n 0.50", "-cN Chr1 Chr2 Chr3 Chr4 Chr5"])
#    self.assertTrue(parser.sampleName)
#
#t = Parser_test()        
#t.test_parser()

# ==== Define and create output directories
outDir = str(parser.locusType) + "_" + str(parser.locusRegion) + "/" + str("_".join(parser.chrName)) + "/"
plotDir = str(outDir + "plots/")

try:
    os.mkdir(outDir)
except OSError as error:
    print(error)

try:
    os.mkdir(plotDir)
except OSError as error:
    print(error)


# ==== Read in Lloyd et al. (2015) Plant Cell lethal-phenotype status for each gene (target, a.k.a. response variable)
def load_target_DF():
    try:
        target = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/" +
                             "TPC2015-00051-LSBR3_Supplemental_Data_set_1_Sheet1.csv",
                             na_values="-")
        return target
    except OSError as error:
        print(error)

ds1_DF = load_target_DF()
ds1_DF.columns = ["gene", "phenotype", "predicted_lethal"]


# ==== Read in Lloyd et al. (2015) Plant Cell gene features, a.k.a predictor variables
def load_features_DF():
    try:
        features = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/" +
                               "TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv",
                               na_values="?")
        return features
    except OSError as error:
        print(error)

ds3_DF = load_features_DF()
ds3_DF.rename(columns={"Locus number":"gene"}, inplace=True)

# ==== Read in gene DNA methylation features 
def load_mC_DF():
    try:
        mC = pd.read_csv(outDir +
                         parser.locusType +
                         "_" + parser.locusRegion +
                         "_" + parser.sampleName +
                         "_MappedOn_" + parser.refbase +
                         "_" + parser.context +
                         "_NAmax" + str(parser.NAmax) +
                         "_filt_df_fk_kappa_all_mean_mC_all_complete_" +
                         str("_".join(parser.chrName)) + ".tsv",
                         sep="\t")
        return mC
    except OSError as error:
        print(error)

mC_DF = load_mC_DF()
mC_DF.rename(columns={
                      "fk_kappa_all":"kappa",
                      "ka_alpha_all":"alpha",
                      "mean_stocha_all":"stocha",
                      "mean_min_acf_all":"min_ACF",
                      "mean_mean_acf_all":"mean_ACF",
                      "mean_mC_all":str("".join("mean_m" + parser.context))
                     },
             inplace=True)
mC_DF["gene"] = mC_DF["name"].str.replace("\.\d+", "", regex=True)
mC_DF["kappa_C_density"] = mC_DF["fk_Cs_all"] / ( ( mC_DF["end"] - mC_DF["start"] + 1 ) / 1e3 )
mC_DF["C_density"] = mC_DF["stocha_Cs_all"] / ( ( mC_DF["end"] - mC_DF["start"] + 1 ) / 1e3 )

kappa_stocha_DF = mC_DF[["gene", "kappa", "stocha", "C_density"]]

# Combine dataframes
merged_DF = pd.merge(ds1_DF, ds3_DF, how="inner", on="gene")
merged_DF = merged_DF[merged_DF["phenotype"].isin(["Lethal", "Non-Lethal"])]
merged_DF.phenotype.replace(["Non-Lethal", "Lethal"], [0, 1], inplace=True)
## NOT DONE AS DROPS ALL ROWS: Drop rows containing missing values across ANY of the features (predictors) ??
#merged_DF.dropna(axis=0, inplace=True)

merged_DF_features = merged_DF.loc[:, ~merged_DF.columns.isin(["gene", "phenotype", "predicted_lethal"])]
merged_DF_target = merged_DF["phenotype"]


# ==== Re-encode categorical features
# See https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_categorical.html#sphx-glr-auto-examples-ensemble-plot-gradient-boosting-categorical-py

# Categorical features
merged_DF_features.iloc[:, 23:57].head()
# Continuous features
merged_DF_features.iloc[:,0:23].head()

## Get column indices of categorical features
#col_idx_categorical_features = list(range(23, 57))

# Get column names of categorical features
categorical_columns = list(merged_DF_features.iloc[:, 23:57].columns)

# Get column names of continuous features
continuous_columns = list(merged_DF_features.iloc[:, 0:23].columns) 


# ==== Remove features that show multi-collinearity, as indicated by high variance-inflation factors (VIFs)
# This occurs when a predictor variable (feature) exhibits a linear relationship
# with two or more other predictors (features)
# "VIF calculates how much the variance of a coefficient is inflated because of its
# linear dependencies with other predictors."
# See https://towardsdatascience.com/statistics-in-python-collinearity-and-multicollinearity-4cc4dcd82b3f
def calculate_vif(df, datatype):
    vif, tolerance = {}, {}
    #
    if datatype=="categorical":
        features = categorical_columns
    elif datatype=="continuous":
        features = continuous_columns
    else:
        print("datatype must be either 'categorical' or 'continuous'")
    #
    # All the features to be evaluated
    for feature in features:
        df = df.dropna(subset=[feature], inplace=False)
        #
        # Extract all the other features to be regressed against
        X = [f for f in features if f != feature]
        X, y = df[X], df[feature]
        #
        # Extract R-squared from the fit
        if datatype=="categorical":
            r2 = HistGradientBoostingClassifier().fit(X, y).score(X, y)
        elif datatype=="continuous":
            r2 = HistGradientBoostingRegressor().fit(X, y).score(X, y)
        else:
            print("datatype must be either 'categorical' or 'continuous'")
        #
        # Calculate tolerance
        tolerance[feature] = 1 - r2
        #
        # Calculate VIF
        vif[feature] = 1/(tolerance[feature])
    #
    # Return VIF DataFrame
    return pd.DataFrame({"VIF": vif, "Tolerance": tolerance})

vif_categorical = calculate_vif(df=merged_DF_features, datatype="categorical") 
vif_continuous = calculate_vif(df=merged_DF_features, datatype="continuous")

vif_categorical_pass = list(vif_categorical[vif_categorical["VIF"] < 10].index)
vif_continuous_pass = list(vif_continuous[vif_continuous["VIF"] < 10].index)

# Append categorical_columns (left) and continuous_columns (right)
merged_DF_features = merged_DF_features[vif_categorical_pass + vif_continuous_pass]

# Set categorical_column type to "category" 
merged_DF_features[vif_categorical_pass] = merged_DF_features[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features.select_dtypes(include="category").shape[1]
n_continuous_features = merged_DF_features.select_dtypes(include="number").shape[1]

print(f"Number of rows (loci): {merged_DF_features.shape[0]}")
# Number of rows (loci): 3443
print(f"Number of features (predictors): {merged_DF_features.shape[1]}")
# Number of features (predictors): 57
# Number of features (predictors): 25
print(f"Number of categorical features: {n_categorical_features}")
# Number of categorical features: 34
# Number of categorical features: 12
print(f"Number of continuous features: {n_continuous_features}")
# Number of continuous features: 23
# Number of continuous features: 13

seed=42

X, y = merged_DF_features, merged_DF_target
print(X.shape, y.shape)
# including rows containing NaNs
# (3443, 25) (3443,)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

# === Gradient boosting estimator with native categorical support, which requires ordinal encoding

# "Create a [HistGradientBoostingClassifier] estimator that will natively handle categorical features.
# "This estimator will not treat categorical features as ordered quantiles.
# "Since the HistGradientBoostingRegressor [or HistGradientBoostingClassifier] requires
# category values to be encoded in [0, n_unique_categories - 1], we still rely on an
# OrdinalEncoder ( https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OrdinalEncoder.html#sklearn.preprocessing.OrdinalEncoder ) to pre-process the data.
# "The main difference between this pipeline and the [one with ordinal encoding only] is that in this one,
# we let the [HistGradientBoostingClassifier] know which features are categorical."
ordinal_encoder = make_column_transformer(
    (
        OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=np.nan),
        make_column_selector(dtype_include="category"),
    ),
    remainder="passthrough",
)

# "The ordinal encoder will first output the categorical features, and then the
# continuous (passed-through) features
categorical_mask = [True] * n_categorical_features + [False] * n_continuous_features

# Make the HistGradientBoostingClassifier estimator
est_native = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask
    ),
)

# K-fold cross validation
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=2,
                             random_state=36851234)

X_np = X.to_numpy()
y_np = y.to_numpy()

for train_index, test_index in cv.split(X_np, y_np):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X_np[train_index], X_np[test_index]
    y_train, y_test = y_np[train_index], y_np[test_index]

n_scores = cross_val_score(est_native, X, y, scoring="accuracy",
                           cv=cv, n_jobs=-1, error_score="raise")


# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native.score(X_test, y_test):.2f}")



## Examples
#XX, yy = load_diabetes(return_X_y=True)
#HistGradientBoostingRegressor().fit(XX, yy).score(XX, yy)
#
#cal_housing = fetch_california_housing()




#knn = KNeighborsClassifier(n_neighbors=15)
## Cannot handle NAs
#knn.fit(X, y)
##ValueError: Input X contains NaN.
##KNeighborsClassifier does not accept missing values encoded as NaN natively. For supervised learning, you might want to consider sklearn.ensemble.HistGradientBoostingClassifier and Regressor which accept missing values encoded as NaNs natively. Alternatively, it is possible to preprocess the data, for instance by using an imputer transformer in a pipeline or drop samples with missing values. See https://scikit-learn.org/stable/modules/impute.html You can find a list of all estimators that handle NaN values at the following page: https://scikit-learn.org/stable/modules/impute.html#estimators-that-handle-nan-values

# ==== ######
def main():



