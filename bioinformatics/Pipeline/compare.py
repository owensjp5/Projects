# Python program for comparing expression between genes
# Jack Owens, Oct 25 2024

# Load Libraries
import sys
import pandas as pd
import numpy as np

# Check Command Line Arguments
assert len(sys.argv) > 1, "No directory given"

# Import Loaded Data
outputDirectory = sys.argv[1]
csvFilePath = outputDirectory + "/out.csv"
df = pd.read_csv(csvFilePath)

# 