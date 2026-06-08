################################################
# Catalan Numbers and RNA Secondary Structures #
################################################
import numpy as np

def catalanNumber(seq):
    complementBase = {"A":"U","U":"A","C":"G","G":"C"}
    n = len(seq)
    dp = np.zeros(shape=(n, n), dtype=object)

    for i in range(n-1):
        if complementBase[seq[i]] == seq[i+1]:
            dp[i, i+1] = 1

    for j in range(3, n, 2):    # window size
        for i in range(n-j):    # start index
            for k in range(i+1, i+j+1, 2): # step through window
                if complementBase[seq[i]] == seq[k]:
                    if k > i+1:
                        left  = dp[i+1, k-1]
                    else:
                        left = 1
                    if k < i+j:
                        right = dp[k+1, i+j]
                    else:
                        right = 1
                    dp[i, i+j] = (dp[i, i+j] + left * right) % 1000000

    return int(dp[0, n-1])

def main():
    print(catalanNumber("CGAUCG"))
    print(catalanNumber("AUAUAU"))
    print(catalanNumber("GUUUUACACGCGUGAUAUAAUCGUGUACAUGUCGUUAAACAUAAAUCGGCAUGGCCAAUAUCGGCUCAAUUGCCAAUUGGUAGGCGCUUAAUACCGACAUGUGUGAUGCGCCAAAGCUGCUACGGCUCGAUGCAUCGUAUAUAUGCGCAGCCGAUAUCCGCGUGCGCAUAACGGUCCGCUAGGAUAGCGCGCUACCAUACGUAUAGGCCAGCUUAGCUCGCGGAUCUAUGCAAUUAUAUAUCGCAUGCUAUAUAGUAUAGCAGGCCUAGCCAUG"))

if __name__ == "__main__":
    main()