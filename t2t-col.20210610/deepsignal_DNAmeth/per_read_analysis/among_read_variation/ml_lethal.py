#!/usr/bin/env python

# Build supervised learning model to predict lethal-phenotype genes

from sklearn.neighbors import KNeighborsClassifier
import pandas as pd


# Read in Lloyd et al. (2015) Plant Cell gene features
ds1 = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/TPC2015-00051-LSBR3_Supplemental_Data_set_1_Sheet1.csv")
ds1.columns = ["gene", "phenotype", "predicted_lethal"]

ds3 = pd.read_csv("Lloyd_2015_Plant_Cell_SupplData/plcell_v27_8_2133_s1/plcell_v27_8_2133_s1/TPC2015-00051-LSBR3_Supplemental_Data_set_3_Sheet1.csv")
ds3.rename(columns = {"Locus number":"gene"}, inplace = True)

df = pd.merge(ds1, ds3, how = "inner", on = "gene")
df.replace(to_replace = "?", value = "NaN", inplace = True)
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

