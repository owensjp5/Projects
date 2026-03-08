import numpy as np
import math

####################################
# Rabbits and Recurrence Relations #
####################################
def recurrenceRelationIterative(n, k):
    cur_reproductive_age = 1
    cur_young = 0
    for i in range(n-2):
        prev_young = cur_young
        prev_reproductive_age = cur_reproductive_age
        cur_young = prev_reproductive_age * k
        cur_reproductive_age = cur_reproductive_age + prev_young
    return cur_reproductive_age + cur_young

def recurrenceRelationDP(n, k):
    buffer = np.zeros(shape=n)
    for i in range(2,n):
        buffer[i] = buffer[i-1] + (buffer[i-2] * k)
    print(buffer)
    return buffer[n-1]

######################
# Mendel's First Law #
######################
def mendelianInheritance(k, m, n):
    total = k + m + n
    probDom1 = k/total + (m/total * 0.5)
    probRec1Dom2 = ((n/total) * (k/(total-1) + m/(total-1)*0.5)) + ((m/total*0.5) * (k/(total-1) + ((m-1)/(total-1)*0.5)))
    probDom = probDom1 + probRec1Dom2
    # probRec = (n/total * ((n-1)/(total-1) + (m/(total-1)*0.5))) + ((m/total*0.5) * ((n/(total-1) + (m-1)/(total-1)*0.5)))
    # print("Probabilities should sum to 1, Sum:", probDom+probRec)
    return probDom

def mendelianInheritanceSimulation(k, m, n, epochs=7, numOffspring=1):
    x = np.array([k, m, n])
    for _ in range(epochs):
        AA = np.array([1/2, 1/2, 0 ]) * x[0]
        Aa = np.array([1/4, 1/2, 1/4]) * x[1]
        aa = np.array([0, 1/2, 1/2]) * x[2]
        x = numOffspring * (AA + Aa + aa)
    print(x)

############################
# Mortal Fibonacci Rabbits #
############################
def mortalRecurrenceRelation(n, m, k=1):
    buffer1 = [0] * (n+1)
    buffer2 = [0] * (n+1)
    buffer1[1] = k
    buffer2[0] = k
    for i in range(2,n+1):
        buffer2[i-1] = buffer1[i-1] - buffer2[i-2]
        buffer1[i] = buffer1[i-1] + (buffer1[i-1] - buffer2[i-2]) - buffer2[i-m-1]
    print("Total Rabbits: ", buffer1)
    return buffer1[n]

#######################
# Independent Alleles #
#######################
def independentInheritance(k, N):
    pass

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

def main():
    # print(recurrenceRelationIterative(35, 5))
    # print(mendelianInheritance(29,24,21))
    # print(mortalRecurrenceRelation(93, 20))
    print(independentInheritance(2,1))
    print(binomialIndependentInheritance(5,8))

if __name__ == "__main__":
    main()