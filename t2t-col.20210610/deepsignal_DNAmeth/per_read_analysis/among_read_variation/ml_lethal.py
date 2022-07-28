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
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot as plt

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
mC_DF["gene"] = mC_DF["gene"].str.replace("_\d+", "", regex=True)

mC_DF["kappa_C_density"] = mC_DF["fk_Cs_all"] / ( ( mC_DF["end"] - mC_DF["start"] + 1 ) / 1e3 )
mC_DF["C_density"] = mC_DF["stocha_Cs_all"] / ( ( mC_DF["end"] - mC_DF["start"] + 1 ) / 1e3 )

kappa_DF = mC_DF[["gene", str("".join("mean_m" + parser.context)), "kappa"]]
#, "kappa_C_density"]]
alpha_DF = mC_DF[["gene", str("".join("mean_m" + parser.context)), "alpha"]]
#, "C_density"]]
stocha_DF = mC_DF[["gene", str("".join("mean_m" + parser.context)), "stocha"]]
#, "C_density"]]


# Combine dataframes
list_DF = [ds1_DF, kappa_DF, ds3_DF]

merged_DF = reduce(lambda  left, right: pd.merge(left, right,
                                                 on=["gene"], how="outer"),
                   list_DF)

merged_DF = merged_DF[merged_DF["phenotype"].isin(["Lethal", "Non-Lethal"])]
merged_DF.phenotype.replace(["Non-Lethal", "Lethal"], [0, 1], inplace=True)
## NOT DONE AS DROPS ALL ROWS: Drop rows containing missing values across ANY of the features (predictors) ??
#merged_DF.dropna(axis=0, inplace=True)

merged_DF_features = merged_DF.loc[:, ~merged_DF.columns.isin(["gene", "phenotype", "predicted_lethal"])]
merged_DF_target = merged_DF["phenotype"]


# ==== Re-encode categorical features
# See https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_categorical.html#sphx-glr-auto-examples-ensemble-plot-gradient-boosting-categorical-py

# Categorical features
merged_DF_features.iloc[:, 25:59].head()
# Continuous features
merged_DF_features.iloc[:,0:25].head()

## Get column indices of categorical features
#col_idx_categorical_features = list(range(23, 57))

# Get column names of categorical features
categorical_columns = list(merged_DF_features.iloc[:, 25:59].columns)

# Get column names of continuous features
continuous_columns = list(merged_DF_features.iloc[:, 0:25].columns) 


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

## Make single-feature prediction models and get ROC AUC for each
#def feature_pred(feat, targ, datatype):
#    #
#    if datatype=="categorical":
#        features = vif_categorical_pass 
#    elif datatype=="continuous":
#        features = vif_continuous_pass
#    else:
#        print("datatype must be either 'categorical' or 'continuous'")
#    #
#    # All the features to be evaluated
#    for feature in features:
#        print(feat[feature].name)
#        X, y = feat[feature], targ
#        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)
#        # Generate a no-skill classifier prediction (majority class)
#        ns_probs = [0 for _ in range(len(y_test))]
#        # Set categorical mask for the HistGradientBoostingClassifier
#        if datatype=="categorical":
#            # Make the HistGradientBoostingClassifier estimator
#            est_native = make_pipeline(
#                ordinal_encoder,
#                HistGradientBoostingClassifier(
#                    random_state=seed,
#                    categorical_features=[True]
#                ),
#            )
#        elif datatype=="continuous":
#            est_native = HistGradientBoostingClassifier(random_state=seed) 
#        else:
#            print("datatype must be either 'categorical' or 'continuous'")
#
## Train the estimator on the training set (90% of loci)
#print("Training HistGradientBoostingClassifier...")
#tic = time()
#est_native.fit(X_train, y_train)
#print(f"Training done in {time() - tic:.3f}s")
## Evaluate performance on the test set (the other 10% of loci)
#print(f"Test R2 score: {est_native.score(X_test, y_test):.2f}")

### !!! ### !!!
vif_categorical_pass = ['Core eukaryotic gene', 'Gene body methylated', 'GOslim C nucleus', 'alpha WGD paralog retained', 'GOslim F protein binding', 'GOslim C plastid', 'GOslim P cellular component organization', 'No homolog in rice', 'GOslim P nucleobase-containing compound metabolic process', 'GOslim C plasma membrane', 'GOslim P response to abiotic stimulus', 'GOslim P reproduction'] 
vif_continuous_pass = ['mean_mCpG', 'kappa', 'Median expression', 'PPIs - AIMC', 'Expression correlation - Ks below 2', 'No. of amino acids in protein', 'Ks with putative paralog', 'Nucleotide diversity', 'Sequence conservation in plants (% ID)', 'Expression variation', 'Co-expression module size', 'Percent identity with putative paralog', 'Ka/Ks with putative paralog', 'Expression breadth', 'OrthoMCL paralog cluster size']
### !!! ### !!!

#vif_categorical = calculate_vif(df=merged_DF_features, datatype="categorical") 
#vif_continuous = calculate_vif(df=merged_DF_features, datatype="continuous")
#
#vif_categorical_pass = list(vif_categorical[vif_categorical["VIF"] < 10].index)
#vif_continuous_pass = list(vif_continuous[vif_continuous["VIF"] < 10].index)

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
# Number of features (predictors): 27
print(f"Number of categorical features: {n_categorical_features}")
# Number of categorical features: 34
# Number of categorical features: 12
print(f"Number of continuous features: {n_continuous_features}")
# Number of continuous features: 23
# Number of continuous features: 15

seed=42

X, y = merged_DF_features, merged_DF_target
print(X.shape, y.shape)
# including rows containing NaNs
# (3443, 27) (3443,)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

# Generate a no-skill classifier prediction (majority class)
# See https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
ns_probs = [0 for _ in range(len(y_test))]


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
# continuous (passed-through) features"
categorical_mask = [True] * n_categorical_features + [False] * n_continuous_features

# Make the HistGradientBoostingClassifier estimator
est_native = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask
    ),
)

## K-fold cross validation
#cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=2,
#                             random_state=36851234)
#
#X_np = X.to_numpy()
#y_np = y.to_numpy()
#
#for train_index, test_index in cv.split(X_np, y_np):
#    print("TRAIN:", train_index, "TEST:", test_index)
#    X_train, X_test = X_np[train_index], X_np[test_index]
#    y_train, y_test = y_np[train_index], y_np[test_index]
#
#n_scores = cross_val_score(est_native, X, y, scoring="accuracy",
#                           cv=cv, n_jobs=-1, error_score="raise")


# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native.score(X_test, y_test):.2f}")
r2 = est_native.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs = est_native.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs = hgb_probs[:, 1]
# Calculate scores
ns_auc = roc_auc_score(y_test, ns_probs)
hgb_auc = roc_auc_score(y_test, hgb_probs) 
# Summarise scores
print("No-skill classifier: ROC AUC=%.3f" % (ns_auc))
print("Full-model classifier: ROC AUC=%.3f" % (hgb_auc))
# Calculate ROC curves
ns_fpr, ns_tpr, _ = roc_curve(y_test, ns_probs)
hgb_fpr, hgb_tpr, _ = roc_curve(y_test, hgb_probs)


# ==== Repeat for full model less kappa

# Append categorical_columns (left) and continuous_columns (right)
vif_continuous_pass_no_kappa = [vif_continuous_pass[0]] + vif_continuous_pass[2:]
merged_DF_features_no_kappa = merged_DF_features[vif_categorical_pass + vif_continuous_pass_no_kappa]

# Set categorical_column type to "category" 
merged_DF_features_no_kappa[vif_categorical_pass] = merged_DF_features_no_kappa[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features_no_kappa.select_dtypes(include="category").shape[1]
n_continuous_features_no_kappa = merged_DF_features_no_kappa.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_kappa, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_kappa = [True] * n_categorical_features + [False] * n_continuous_features_no_kappa

# Make the HistGradientBoostingClassifier estimator
est_native_no_kappa = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_kappa
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_kappa.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_kappa.score(X_test, y_test):.2f}")
r2_no_kappa = est_native_no_kappa.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_kappa = est_native_no_kappa.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_kappa = hgb_probs_no_kappa[:, 1]
# Calculate scores
hgb_auc_no_kappa = roc_auc_score(y_test, hgb_probs_no_kappa) 
# Summarise scores
print("Full-model-less-kappa classifier: ROC AUC=%.3f" % (hgb_auc_no_kappa))
# Calculate ROC curves
hgb_fpr_no_kappa, hgb_tpr_no_kappa, _ = roc_curve(y_test, hgb_probs_no_kappa)


# ==== Repeat for full model less coexpress_mod_size

# Append categorical_columns (left) and continuous_columns (right)
vif_continuous_pass_no_coexpress_mod_size = vif_continuous_pass[0:10] + vif_continuous_pass[11:]
merged_DF_features_no_coexpress_mod_size = merged_DF_features[vif_categorical_pass + vif_continuous_pass_no_coexpress_mod_size]

# Set categorical_column type to "category" 
merged_DF_features_no_coexpress_mod_size[vif_categorical_pass] = merged_DF_features_no_coexpress_mod_size[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features_no_coexpress_mod_size.select_dtypes(include="category").shape[1]
n_continuous_features_no_coexpress_mod_size = merged_DF_features_no_coexpress_mod_size.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_coexpress_mod_size, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_coexpress_mod_size = [True] * n_categorical_features + [False] * n_continuous_features_no_coexpress_mod_size

# Make the HistGradientBoostingClassifier estimator
est_native_no_coexpress_mod_size = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_coexpress_mod_size
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_coexpress_mod_size.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_coexpress_mod_size.score(X_test, y_test):.2f}")
r2_no_coexpress_mod_size = est_native_no_coexpress_mod_size.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_coexpress_mod_size = est_native_no_coexpress_mod_size.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_coexpress_mod_size = hgb_probs_no_coexpress_mod_size[:, 1]
# Calculate scores
hgb_auc_no_coexpress_mod_size = roc_auc_score(y_test, hgb_probs_no_coexpress_mod_size) 
# Summarise scores
print("Full-model-less-coexpress_mod_size classifier: ROC AUC=%.3f" % (hgb_auc_no_coexpress_mod_size))
# Calculate ROC curves
hgb_fpr_no_coexpress_mod_size, hgb_tpr_no_coexpress_mod_size, _ = roc_curve(y_test, hgb_probs_no_coexpress_mod_size)


# ==== Repeat for full model less paralog_ID

# Append categorical_columns (left) and continuous_columns (right)
vif_continuous_pass_no_paralog_ID = vif_continuous_pass[0:11] + vif_continuous_pass[12:]
merged_DF_features_no_paralog_ID = merged_DF_features[vif_categorical_pass + vif_continuous_pass_no_paralog_ID]

# Set categorical_column type to "category" 
merged_DF_features_no_paralog_ID[vif_categorical_pass] = merged_DF_features_no_paralog_ID[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features_no_paralog_ID.select_dtypes(include="category").shape[1]
n_continuous_features_no_paralog_ID = merged_DF_features_no_paralog_ID.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_paralog_ID, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_paralog_ID = [True] * n_categorical_features + [False] * n_continuous_features_no_paralog_ID

# Make the HistGradientBoostingClassifier estimator
est_native_no_paralog_ID = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_paralog_ID
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_paralog_ID.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_paralog_ID.score(X_test, y_test):.2f}")
r2_no_paralog_ID = est_native_no_paralog_ID.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_paralog_ID = est_native_no_paralog_ID.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_paralog_ID = hgb_probs_no_paralog_ID[:, 1]
# Calculate scores
hgb_auc_no_paralog_ID = roc_auc_score(y_test, hgb_probs_no_paralog_ID) 
# Summarise scores
print("Full-model-less-paralog_ID classifier: ROC AUC=%.3f" % (hgb_auc_no_paralog_ID))
# Calculate ROC curves
hgb_fpr_no_paralog_ID, hgb_tpr_no_paralog_ID, _ = roc_curve(y_test, hgb_probs_no_paralog_ID)


# ==== Repeat for full model less paralog_clust_size

# Append categorical_columns (left) and continuous_columns (right)
vif_continuous_pass_no_paralog_clust_size = vif_continuous_pass[0:14]
merged_DF_features_no_paralog_clust_size = merged_DF_features[vif_categorical_pass + vif_continuous_pass_no_paralog_clust_size]

# Set categorical_column type to "category" 
merged_DF_features_no_paralog_clust_size[vif_categorical_pass] = merged_DF_features_no_paralog_clust_size[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features_no_paralog_clust_size.select_dtypes(include="category").shape[1]
n_continuous_features_no_paralog_clust_size = merged_DF_features_no_paralog_clust_size.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_paralog_clust_size, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_paralog_clust_size = [True] * n_categorical_features + [False] * n_continuous_features_no_paralog_clust_size

# Make the HistGradientBoostingClassifier estimator
est_native_no_paralog_clust_size = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_paralog_clust_size
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_paralog_clust_size.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_paralog_clust_size.score(X_test, y_test):.2f}")
r2_no_paralog_clust_size = est_native_no_paralog_clust_size.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_paralog_clust_size = est_native_no_paralog_clust_size.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_paralog_clust_size = hgb_probs_no_paralog_clust_size[:, 1]
# Calculate scores
hgb_auc_no_paralog_clust_size = roc_auc_score(y_test, hgb_probs_no_paralog_clust_size) 
# Summarise scores
print("Full-model-less-paralog_clust_size classifier: ROC AUC=%.3f" % (hgb_auc_no_paralog_clust_size))
# Calculate ROC curves
hgb_fpr_no_paralog_clust_size, hgb_tpr_no_paralog_clust_size, _ = roc_curve(y_test, hgb_probs_no_paralog_clust_size)


# ==== Repeat for full model less express_breadth

# Append categorical_columns (left) and continuous_columns (right)
vif_continuous_pass_no_express_breadth = vif_continuous_pass[0:13] + [vif_continuous_pass[14]]
merged_DF_features_no_express_breadth = merged_DF_features[vif_categorical_pass + vif_continuous_pass_no_express_breadth]

# Set categorical_column type to "category" 
merged_DF_features_no_express_breadth[vif_categorical_pass] = merged_DF_features_no_express_breadth[vif_categorical_pass].astype("category")

n_categorical_features = merged_DF_features_no_express_breadth.select_dtypes(include="category").shape[1]
n_continuous_features_no_express_breadth = merged_DF_features_no_express_breadth.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_express_breadth, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_express_breadth = [True] * n_categorical_features + [False] * n_continuous_features_no_express_breadth

# Make the HistGradientBoostingClassifier estimator
est_native_no_express_breadth = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_express_breadth
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_express_breadth.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_express_breadth.score(X_test, y_test):.2f}")
r2_no_express_breadth = est_native_no_express_breadth.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_express_breadth = est_native_no_express_breadth.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_express_breadth = hgb_probs_no_express_breadth[:, 1]
# Calculate scores
hgb_auc_no_express_breadth = roc_auc_score(y_test, hgb_probs_no_express_breadth) 
# Summarise scores
print("Full-model-less-express_breadth classifier: ROC AUC=%.3f" % (hgb_auc_no_express_breadth))
# Calculate ROC curves
hgb_fpr_no_express_breadth, hgb_tpr_no_express_breadth, _ = roc_curve(y_test, hgb_probs_no_express_breadth)


# ==== Repeat for full model less alpha_WGD

# Append categorical_columns (left) and continuous_columns (right)
vif_categorical_pass_no_alpha_WGD = vif_categorical_pass[0:3] + vif_categorical_pass[4:]
merged_DF_features_no_alpha_WGD = merged_DF_features[vif_categorical_pass_no_alpha_WGD + vif_continuous_pass]

# Set categorical_column type to "category" 
merged_DF_features_no_alpha_WGD[vif_categorical_pass_no_alpha_WGD] = merged_DF_features_no_alpha_WGD[vif_categorical_pass_no_alpha_WGD].astype("category")

n_categorical_features_no_alpha_WGD = merged_DF_features_no_alpha_WGD.select_dtypes(include="category").shape[1]
n_continuous_features = merged_DF_features_no_alpha_WGD.select_dtypes(include="number").shape[1]

seed=42

X, y = merged_DF_features_no_alpha_WGD, merged_DF_target
print(X.shape, y.shape)

# Define training (90% of rows) and test subsets (the other 10% of rows)
# Consulted for use of HistGradientBoostingClassifier() :
# https://scikit-learn.org/stable/auto_examples/inspection/plot_partial_dependence.html#sphx-glr-auto-examples-inspection-plot-partial-dependence-py
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=seed)

categorical_mask_no_alpha_WGD = [True] * n_categorical_features_no_alpha_WGD + [False] * n_continuous_features

# Make the HistGradientBoostingClassifier estimator
est_native_no_alpha_WGD = make_pipeline(
    ordinal_encoder,
    HistGradientBoostingClassifier(
         random_state=seed,
         categorical_features=categorical_mask_no_alpha_WGD
    ),
)

# Train the estimator on the training set (90% of loci)
print("Training HistGradientBoostingClassifier...")
tic = time()
est_native_no_alpha_WGD.fit(X_train, y_train)
print(f"Training done in {time() - tic:.3f}s")
# Evaluate performance on the test set (the other 10% of loci)
print(f"Test R2 score: {est_native_no_alpha_WGD.score(X_test, y_test):.2f}")
r2_no_alpha_WGD = est_native_no_alpha_WGD.score(X_test, y_test)

# For below see https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
# Predict probabilities
hgb_probs_no_alpha_WGD = est_native_no_alpha_WGD.predict_proba(X_test)
# Keep probabilities for the positive outcome only
hgb_probs_no_alpha_WGD = hgb_probs_no_alpha_WGD[:, 1]
# Calculate scores
hgb_auc_no_alpha_WGD = roc_auc_score(y_test, hgb_probs_no_alpha_WGD) 
# Summarise scores
print("Full-model-less-alpha_WGD classifier: ROC AUC=%.3f" % (hgb_auc_no_alpha_WGD))
# Calculate ROC curves
hgb_fpr_no_alpha_WGD, hgb_tpr_no_alpha_WGD, _ = roc_curve(y_test, hgb_probs_no_alpha_WGD)


# Plot the ROC curves
plt.plot(hgb_fpr, hgb_tpr, label=str("Full model: ROC AUC=%.3f" % (hgb_auc)))
plt.plot(hgb_fpr_no_kappa, hgb_tpr_no_kappa, label=str("Excl. kappa: ROC AUC=%.3f" % (hgb_auc_no_kappa)))
plt.plot(hgb_fpr_no_coexpress_mod_size, hgb_tpr_no_coexpress_mod_size, label=str("Excl. coexpress. mod. size: ROC AUC=%.3f" % (hgb_auc_no_coexpress_mod_size)))
plt.plot(hgb_fpr_no_paralog_ID, hgb_tpr_no_paralog_ID, label=str("Excl. paralog identity: ROC AUC=%.3f" % (hgb_auc_no_paralog_ID)))
plt.plot(hgb_fpr_no_paralog_clust_size, hgb_tpr_no_paralog_clust_size, label=str("Excl. paralog cluster size: ROC AUC=%.3f" % (hgb_auc_no_paralog_clust_size)))
plt.plot(hgb_fpr_no_express_breadth, hgb_tpr_no_express_breadth, label=str("Excl. express. breadth: ROC AUC=%.3f" % (hgb_auc_no_express_breadth)))
plt.plot(hgb_fpr_no_alpha_WGD, hgb_tpr_no_alpha_WGD, label=str("Excl. alpha WGD: ROC AUC=%.3f" % (hgb_auc_no_alpha_WGD)))
plt.plot(ns_fpr, ns_tpr, linestyle="--", label=str("No skill: ROC AUC=%.3f" % (ns_auc)))
# Axis labels
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.legend(loc=4)
plt.savefig(outDir +
            "HistGradientBoostingClassifier_lethal.pdf",
            bbox_inches="tight")
plt.close()




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



