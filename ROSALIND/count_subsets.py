####################
# Counting Subsets #
####################
def totalSubsets(n):
    return 2**n % 1000000

def powerSet(n):
    powerSet = [[]]
    for i in range(1,n+1):
        powerSet.append([i])
        if len(powerSet) > 1:
            for j in range(1,len(powerSet)-1): #Start at 1 to ignore empty set
                powerSet.append(powerSet[j]+ [i])
    return powerSet