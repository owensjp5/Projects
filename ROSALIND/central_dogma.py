#!/usr/bin/env python3

import rosalind
import sys

numArgs = len(sys.argv)
if numArgs < 2:
    sys.exit("Error: Missing first command line argument: input file")
if numArgs < 3:
    sys.exit("Error: Missing second command line argument: output directory")

genomeFile = sys.argv[1]
outdir = sys.argv[2]
outputFilePath = "pipeline_outputs/" + outdir + "/output.txt"

DNA_seqs = rosalind.readSequencesFromFasta(genomeFile)

RNA_seqs = {key: rosalind.transcribe(DNA_seq) for key, DNA_seq in DNA_seqs.items()}

PROT_seqs = {}
for key, RNA_seq in RNA_seqs.items():
    ORFs = rosalind.allORFs(RNA_seq, "RNA")
    if ORFs:
        PROT_seqs[key] = ORFs

with open(outputFilePath, 'w') as outputFile:
    prevStdout = sys.stdout
    sys.stdout = outputFile
    for key in DNA_seqs:
        print(key)
        print(DNA_seqs[key])
        print(RNA_seqs[key])
        try:
            print(PROT_seqs[key])
        except:
            print("'No ORFs'")
        print()
    sys.stdout = prevStdout