import numpy as np

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
    buffer[0] = 1
    buffer[1] = 1
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

############################
# Mortal Fibonacci Rabbits #
############################
def mortalRecurrenceRelation(n, m, k=1):
    buffer = np.zeros(shape=(2,n+1))
    buffer[0][1] = k
    buffer[1][0] = k
    for i in range(2,n+1):
        buffer[1][i-1] = buffer[0][i-1] - buffer[1][i-2]
        buffer[0][i] = buffer[0][i-1] + (buffer[0][i-1] - buffer[1][i-2]) - buffer[1][i-m-1]
    print("Total Rabbits: ", buffer[0])
    return buffer[0][n]

def main():
    # print(recurrenceRelationDP(5, 3))
    # print(mendelianInheritance(2,2,2))
    print(mortalRecurrenceRelation(6, 3))

if __name__ == "__main__":
    main()