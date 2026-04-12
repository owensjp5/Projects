##############################
# Locating Restriction Sites #
##############################
from .reverse_complement import reverseComplement

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