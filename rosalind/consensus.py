#########################
# Consensus and Profile #
#########################
from .read_fasta import readSequencesFromFasta
import numpy as np

def profileSequences(fastaFile):
    seqs = list(readSequencesFromFasta(fastaFile).values())
    n = len(seqs[1])
    totalSeqs = len(seqs)

    profile = np.zeros((4,n))
    consensusString = ""

    for i in range(n):
        for j in range(totalSeqs):
            if seqs[j][i] == "A":
                profile[0][i] += 1
            if seqs[j][i] == "C":
                profile[1][i] += 1
            if seqs[j][i] == "G":
                profile[2][i] += 1
            if seqs[j][i] == "T":
                profile[3][i] += 1
        consensus = ""
        max = 0
        if profile[0][i] > max:
            max = profile[0][i]
            consensus = "A"
        if profile[1][i] > max:
            max = profile[1][i]
            consensus = "C"
        if profile[2][i] > max:
            max = profile[2][i]
            consensus = "G"
        if profile[3][i] > max:
            max = profile[3][i]
            consensus = "T"
        consensusString += consensus
    print(consensusString)
    return profile
