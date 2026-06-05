##############################
# Creating a Distance Matrix #
##############################
from .hamming_distance import getHammingDistance
from .read_fasta import readSequencesFromFasta

def distanceMatrixFromFasta(fastaFile):
    seqs = readSequencesFromFasta(fastaFile)
    seqs = list(seqs.values())
    distanceMatrix(seqs)

def distanceMatrix(seqs):
    n = len(seqs)
    m = len(seqs[0])
    dp = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            dp[i][j] = round(getHammingDistance(seqs[i], seqs[j]) / m, 5)
    for i in range(n):
        for j in range(n):
            print(dp[i][j], end=" ")
        print()

def main():
    distanceMatrixFromFasta("rosalind/reads.fasta")

if __name__ == "__main__":
    main()