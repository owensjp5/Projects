#!/bin/bash

inputFilePath=$1
outdir=$2

mkdir -p outputs/$outdir

python3 central_dogma.py $inputFilePath $outdir