# Python program to read gene expression data from .txt file(s) and load into CSV
# Jack Owens, Oct 24 2024

# Load Libraries
import sys
import pandas as pd

# Check Command Line Arguments
assert len(sys.argv) > 1, "No directory given"
assert len(sys.argv) > 2, "No files given"

# Initialize Dataframe
genes=sys.argv[2:]
df = pd.DataFrame(columns=genes)
df.index.name = "REPORTER_ID"

# Read Text Data
for gene in genes:
    filePath = "E-GEOD-32367-small/" + gene + "_sample_table.txt"
    with open(filePath) as dataFile:
        next(dataFile)
        for line in dataFile:
            values = line.split()
            label = values[0]
            df.loc[label, gene] = values[1]

# Format Data
df = df.sort_index()

# Output Data into CSV
outputDirectory = sys.argv[1]
outputFilePath = outputDirectory + "/out.csv"
df.to_csv(outputFilePath)