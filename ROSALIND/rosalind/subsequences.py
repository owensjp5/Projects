##################################
# Longest Increasing Subsequence #
##################################
import numpy as np
import math

def longestIncreasingSubsequence(n, sequence):
    dp = np.ones(n+1)
    sequence.append(math.inf)
    traceback = [0]*(n+1)
    for i in range(n+1):
        for j in range(i):
            if sequence[i] > sequence[j] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                traceback[i] = j
    answer = []
    cur = traceback[n]
    for _ in range(int(dp[cur])):
        answer.append(sequence[cur])
        cur = traceback[cur]
    return answer[::-1]

def longestDecreasingSubsequence(n, sequence):
    dp = np.ones(n+1)
    sequence.append(-math.inf)
    traceback = [0]*(n+1)
    for i in range(n+1):
        for j in range(i):
            if sequence[i] < sequence[j] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                traceback[i] = j
    answer = []
    cur = traceback[n]
    for _ in range(int(dp[cur])):
        answer.append(sequence[cur])
        cur = traceback[cur]
    return answer[::-1]

###########################
# Finding a Spliced Motif #
###########################

def subsequence(seq, subseq):
    ans = []
    subseqLength = len(subseq)
    subseqPos = 0
    for seqPos in range(len(seq)):
        if seq[seqPos] == subseq[subseqPos]:
            ans.append(seqPos+1)
            subseqPos += 1
            if subseqPos == subseqLength:
                return ans

def main():
    print(subsequence("ACGTACGTGACG", "GTA"))

if __name__ == "__main__":
    main()