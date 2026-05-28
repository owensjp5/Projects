#!/bin/bash

inputFilePath=$1
outdir=$2

mkdir -p pipeline_outputs/$outdir
touch pipeline_outputs/$outdir/output.txt

python3 central_dogma.py $inputFilePath $outdir