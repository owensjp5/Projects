###############################################
# Wobble Bonding and RNA Secondary Structures #
###############################################
def wobbleBonding(seq):
    complementBase = {"A":"U","U":"A","C":"G","G":"C"}
    wobbleBase = {"U":"G","G":"U","A":'-',"C":'-'}
    n = len(seq)
    dp = [[1] * (n+2) for _ in range(n+2)]

    for j in range(1, n):
        for i in range(n-j):
            total = dp[i+2][i+j+1]
            for k in range(i+4,i+j+1):
                if complementBase[seq[i]] == seq[k] or wobbleBase[seq[i]] == seq[k]:
                    total += dp[i+2][k] * dp[k+2][i+j+1]
            dp[i+1][i+j+1] = total

    return dp[1][n]

def main():
    print(wobbleBonding("AGGGCCACCAAAGGGGCGGUACGCCACACGAAAAUUCUGAACCACUGAAACAGCAAAGGCCCUCAAAUGAGAAACUAGGGUGCAUACAGAGUCUGCGACAGCAAUUAAUGUGGGUGAUGUCCGAAGGGGUACUAGGCGCGAGCUUCACUAACUUAAUGCAGC"))
    
if __name__ == "__main__":
    main()