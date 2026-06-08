#####################
# Reversal Distance #
#####################
def reverseSequence(sequence):
    newSeq = []
    for i in range(len(sequence)):
        newSeq.append(sequence[-i-1])
    return newSeq

def reverseSubSequence(sequence, start, end):
    if start == end:
        return sequence
    return sequence[:start] + reverseSequence(sequence[start:end+1]) + sequence[end+1:]

def reversalDistance(permutation1, permutation2):
    n = len(permutation1)
    if n != len(permutation2):
        print("permutations must be same length")
        return -1

    if permutation1 == permutation2:
        return 0
    
    reversals = [(permutation1, 0)]
    visited = {tuple(permutation1)}

    while True:
        newReversals = []
        for k in range(len(reversals)):
            for i in range(n):
                for j in range(i+1,n):
                    newPerm = reverseSubSequence(reversals[k][0],i,j)
                    if newPerm == permutation2:
                        return reversals[k][1]+1
                    tupleNewPerm = tuple(newPerm)
                    if tupleNewPerm not in visited:
                        visited.add(tupleNewPerm)
                        newReversals.append((newPerm, reversals[k][1]+1))
        reversals = newReversals

def main():
    print(reversalDistance([1,2,3,4,5,6,7,8,9,10],[3,1,5,2,7,4,9,6,10,8]))
    # print(reverseSubSequence([1,2,3,4],0,1))

if __name__ == "__main__":
    main()