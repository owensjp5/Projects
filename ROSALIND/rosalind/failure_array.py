#############################
# Speeding Up Motif Finding #
#############################
from .read_fasta import readSequencesFromFasta

def failureArray(seq):
    n = len(seq)
    dp = [0] * n
    for i in range(1,n):
        for j in range(i,n):
            if seq[j] == seq[j-i]:
                dp[j] = max(j-i+1, dp[j])
            else:
                break
    return dp

def main():
    seqs = readSequencesFromFasta("rosalind_kmp.txt")
    seq = seqs[">Rosalind_4708"]

    failure_array = failureArray("CAGCATGGTATCACAGCAGAG")
    for i in range(len(failure_array)):
        print(failure_array[i], end=" ")

if __name__ == "__main__":
    main()