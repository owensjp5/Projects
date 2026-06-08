####################################
# Rabbits and Recurrence Relations #
####################################
import numpy as np

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