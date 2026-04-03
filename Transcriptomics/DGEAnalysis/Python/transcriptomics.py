import csv
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import seaborn as sns

# read in gene expression data
df = pd.read_csv('/Users/jackowens/desktop/Projects/Transcriptomics/DGEAnalysis/ref/GSE150910_gene-level_count_file.csv')

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

# K-S Tests evaluate likelihood that two sample distributions were drawn from the same distribution
p_values = []
log2_FCs = []

for i in range(data.shape[0]):
    control = data[i, labels == "control"]
    ipf = data[i, labels == "ipf"]

    ks_statistic, p_value = ks_2samp(control, ipf)

    p_values.append(p_value)

    epsilon = 1
    fold_change = (np.mean(ipf) + epsilon) / (np.mean(control) + epsilon)
    log2_FC = np.log2(fold_change)

    log2_FCs.append(log2_FC)

p_values = np.array(p_values)
log2_FCs = np.array(log2_FCs)

# Bonferroni, multiple test correction
p_values_bonf = p_values * len(genes) #multiplying p_values by # of genes helps account for random chance in large dataset
bonf_threshold = 0.05
# threshold = 0.05 / len(genes)
to_keep = (p_values_bonf <= bonf_threshold) & (np.abs(log2_FCs) >= 1)
# to_keep = (p_values <= threshold) & (np.abs(log2_FCs) >= 1) #Alternate way to perform bonf

sig_genes = genes[to_keep]
sig_log2_FCs = log2_FCs[to_keep]
sig_p_values = p_values_bonf[to_keep]
sig_data = data[to_keep, :]

print("Num DEGs:", len(sig_genes))
print("Total Genes:", len(genes))
print("DEG Percentage:", len(sig_genes) / len(genes) * 100)
print()
print("Number of upregulated DEGs:", np.sum(sig_log2_FCs>0))
print("Number of downregulated DEGs:", np.sum(sig_log2_FCs<0))

# Save results to CSV      
output_file = "output/upregulated_DEGs.csv"
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    writer.writerow(["gene", "p_value", "log2_FC"])

    for i in range(len(sig_genes)):
        if sig_log2_FCs[i] > 0:
            writer.writerow([sig_genes[i], sig_p_values[i], sig_log2_FCs[i]])
     
output_file = "output/downregulated_DEGs.csv"
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    writer.writerow(["gene", "p_value", "log2_FC"])

    for i in range(len(sig_genes)):
        if sig_log2_FCs[i] < 0:
            writer.writerow([sig_genes[i], sig_p_values[i], sig_log2_FCs[i]])


# Generate heatmap plot
sig_data_log2FCs = np.zeros_like(sig_data)

