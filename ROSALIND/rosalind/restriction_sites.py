##############################
# Locating Restriction Sites #
##############################
from .reverse_complement import reverseComplement
from .gc_content import getGCContent
import math

def reversePalindromes(seq):
    complementBase = {"A":"T","T":"A","C":"G","G":"C"}
    reversePalindromes = []
    for i in range(len(seq)-3):
        if reverseComplement(seq[i:i+4]) == seq[i:i+4]:
            reversePalindrome = (i,i+3)
            reversePalindromes.append(reversePalindrome)
            if i > 0:
                for j in range(1,i+1):
                    if i+3+j >= len(seq):
                        break
                    if seq[i-j] == complementBase[seq[i+3+j]]:
                        reversePalindromes.append((i-j, i+3+j))
                    else:
                        break
    return reversePalindromes

########################################
# Expected Number of Restriction Sites #
########################################
def expectedRestrictionSites(n, s, A):
    t = len(s)
    m = len(A)
    B = [1] * m
    for i in range(m):
        gc = A[i]
        for base in s:
            if base == "A" or base == "T":
                B[i] = B[i] * (0.5 * (1 - gc))
            else:
                B[i] = B[i] * (0.5 * gc)
        B[i] = math.ceil(B[i] * (n-(t-1)) * 1000) / 1000
    for b in B:
        print(b, end=" ")
    print()