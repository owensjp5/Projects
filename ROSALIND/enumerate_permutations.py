###########################
# Enumerating Gene Orders #
###########################
import math

def permutatations(N):
    total = math.factorial(N)
    print(total)
    base = list(range(1,N+1))
    permute(base, [])
    
def permute(arr, rec):
    if len(arr) != 1:
        for i in range(len(arr)):
            base = arr.copy()
            tmp = base[0]
            base[0] = base[i]
            base[i] = tmp
            permute(base[1:], rec=rec+[base[0]])
    else:
        print(*(rec + arr))