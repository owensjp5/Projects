##################################################
# Perfect Matchings and RNA Secondary Structures #
##################################################
import math

def perfectMatchings(seq):
    occ = {"A":0,"C":0,"G":0,"U":0}
    length = len(seq)
    for i in range(length):
        occ[seq[i]] += 1
    return math.factorial(occ["C"]) * math.factorial(occ["A"])

##################################################
# Maximum Matchings and RNA Secondary Structures #
##################################################
def maximumMatchings(seq):
    occ = {"A":0,"C":0,"G":0,"U":0}
    length = len(seq)
    for i in range(length):
        occ[seq[i]] += 1
    max_GC_pairs = max(occ["G"], occ["C"])
    min_GC_pairs = min(occ["G"],occ["C"])
    max_AU_pairs = max(occ["A"], occ["U"])
    min_AU_pairs = min(occ["A"],occ["U"])
    possible_GC_pairs = math.perm(max_GC_pairs, min_GC_pairs)
    possible_AU_pairs = math.perm(max_AU_pairs, min_AU_pairs)
    return possible_GC_pairs * possible_AU_pairs