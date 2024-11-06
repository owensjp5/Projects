#!/bin/bash

# Create Output Directory
outdir=$1
mkdir -p $outdir

# Parse Command Line Arguments
gene1=$2
gene2=$3

# Pass Arguments to Each Stage of the Pipeline
python load.py $outdir $gene1 $gene2
python compare.py $outdir
