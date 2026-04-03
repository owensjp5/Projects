##################################################
# Perfect Matchings and RNA Secondary Structures #
##################################################
import math

def perfectMatchings(seq):
    occ = {"A":0,"C":0,"G":0,"U":0}
    length = len(seq)
    for i in range(length):
        occ[seq[i]] += 1
    GC_pairs = min(occ["G"],occ["C"])
    AU_pairs = min(occ["A"],occ["U"])
    return math.factorial(GC_pairs) * math.factorial(AU_pairs)