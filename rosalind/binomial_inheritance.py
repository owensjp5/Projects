#######################
# Independent Alleles #
#######################
import numpy as np
import math

def simulateIndependentInheritance(k, N):
    x = np.array([0,1,0])
    M = np.array([[0.5,0.5,0],[0.25,0.5,0.25],[0,0.5,0.5]])
    for _ in range(k):
        x = x @ M
    return x

def binomialIndependentInheritance(k, N):
    n = 2**k
    x = 1/4
    y = 3/4
    if k > 7:
        return "k must be ≤ 7"
    if N > n:
        return "N must be ≤ 2^k"
    sum = 0
    for i in range(N,n+1):
        sum += math.comb(n,i) * x**i * y**(n-i)
    return sum