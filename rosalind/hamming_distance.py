############################
# Counting Point Mutations #
############################
def getHammingDistance(seq1, seq2):
    length1 = len(seq1)
    length2 = len(seq2)
    if length1 > 1000 or length2 > 1000:
        print("Error: improper input for getHammingDistance fn, input contains sequence longer than 1 kbp")
        return
    if length1 != length2:
        print("Error: improper input for getHammingDistance fn, sequences are different length")
        return
    hammingDistance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hammingDistance += 1
    return hammingDistance