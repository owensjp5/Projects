########################
# Partial Permutations #
########################
import math

def partialPermutations(n, k):
    numerator = math.factorial(n)
    denominator = math.factorial(n-k)
    return int(numerator // denominator) % 1000000