#####################
# k-Mer Composition #
#####################
from .read_fasta import readSequencesFromFasta

def kMerComposition(s, k):
    n = len(s)
    total = 4 ** k #4 possible bases for all k indices
    A = [0] * total

    ranks = {'A':0, 'C':1, 'G':2, 'T':3}

    hash = 0
    for base in s[:k]:
        hash = hash * 4 + ranks[base] #Multiply previous hash by 4 so each base has unique significance (associated power of 4)

    A[hash] = 1
    high = 4 ** (k-1)

    for i in range(k, n):
        hash = (hash - ranks[s[i-k]] * high) * 4 + ranks[s[i]] #Subtract hash of previous leftmost digit, increment significance(power of 4), and add hash of new rightmost digit
        A[hash] += 1

    return A

def main():
    read = readSequencesFromFasta("read.txt")
    print(read[">Rosalind_0523"])
    A = kMerComposition(read[">Rosalind_0523"], 4)
    for i in range(len(A)):
        print(A[i], end=" ")
    print()

if __name__ == "__main__":
    main()