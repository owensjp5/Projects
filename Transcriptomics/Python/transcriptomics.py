import csv
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import seaborn as sns

# read in gene expression data

df = pd.read_csv('/Users/jackowens/desktop/Projects/Transcriptomics/ref/GSE150910_gene-level_count_file.csv')

genes = []
data = []

for row in df.itertuples(index=False):
    genes.append(row[0])
    data.append(row[1:])

genes = np.array(genes)
samples = np.array(df.columns[1:].to_numpy())
data = np.array(data).astype(float)

labels = []

for i in range(len(samples)):
    tmp = samples[i].split("_")
    labels.append(tmp[0])

labels = np.array(labels)

# remove chp samples

data = data[:, labels != "chp"]
samples = samples[labels != "chp"]
labels = labels[labels != "chp"]

# normalize data (counts per million)

for j in range(data.shape[1]):
    column_sum = sum(data[:,j])
    data[:,j] = data[:,j] / column_sum * 10**6

# filter out low-expression genes

mean_CPM_control = data[:,labels == "control"].mean(axis=1)
mean_CPM_ipf = data[:,labels == "ipf"].mean(axis=1)

to_keep = (mean_CPM_ipf >= 5) | (mean_CPM_control >= 5)

data = data[to_keep, :]
genes = genes[to_keep]