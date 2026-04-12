###########################
# Enumerating Gene Orders #
###########################
import math

def permutations(N):
    total = math.factorial(N)
    base = list(range(1,N+1))
    permutationList = []
    permute(base, [], permutationList)
    return total, permutationList
    
def permute(arr, rec, permutationList):
    if len(arr) != 1:
        for i in range(len(arr)):
            base = arr.copy()
            tmp = base[0]
            base[0] = base[i]
            base[i] = tmp
            permute(base[1:], rec=rec+[base[0]], permutationList=permutationList)
    else:
        permutationList.append(rec + arr)

#######################################
# Enumerating Oriented Gene Orderings #
#######################################
def signedPermutations(N):
    signedPermutationList = []
    total, permutatationList = permutations(N)
    total = total * 2**N
    for permutation in permutatationList:
        signedPermutation = [permutation]
        signedPermute(signedPermutation, N, N)
        signedPermutationList += signedPermutation
    return total, signedPermutationList

def signedPermute(permutationList, N, k):
    for i in range(len(permutationList)):
        negative = permutationList[i].copy()
        negative[N-k] = -negative[N-k]
        permutationList.append(negative)
    if k > 1:
        signedPermute(permutationList, N, k-1)

def main():
    # total, permutationList = permutations(3)
    # print(total)
    # for permutation in permutationList:
    #     print(*(permutation))

    total, signedPermutationList = signedPermutations(2)
    print(total)
    for signedPermutation in signedPermutationList:
        print(*(signedPermutation))

if __name__ == "__main__":
    main()