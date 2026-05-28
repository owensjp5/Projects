##########################
# Finding a Shared Motif #
##########################
from .read_fasta import readSequencesFromFasta

def sharedMotif(fastaFile):
    seqs, max_len, longest_seq = readSequencesFromFasta(fastaFile, True)
    
    motifs = {}

    for i in range(max_len+1):
        for j in range(i):
            if seqs[longest_seq][j:i] not in motifs:
                motifs[seqs[longest_seq][j:i]] = False
            shared = True
            for seq in seqs:
                if seq != longest_seq:
                    if seqs[longest_seq][j:i] not in seqs[seq]:
                        shared = False
                        break
            if shared:
                motifs[seqs[longest_seq][j:i]] = True
    
    motifs = sorted(list(motifs.items()), key = lambda item: len(item[0]), reverse=True)

    for motif in motifs:
        if motif[1]:
            return motif[0]
        
    return "No shared motifs"

##################################
# Finding a Shared Spliced Motif #
##################################
def longestCommonSubsequence(seq1, seq2):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n+1) for _ in range(m+1)]
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                dp[i+1][j+1] = dp[i][j]+1
            else:
                dp[i+1][j+1] = max(dp[i][j+1],dp[i+1][j])

    # Display Memoization Table            
    # print("", end="      ")
    # for char in seq2:
    #     print(char, end="  ")
    # print()
    # print(" ", dp[0])
    # for i in range(1, m+1):
    #     print(seq1[i-1], dp[i])

    xpos = m
    ypos = n
    traceback = ""
    while xpos > 0 and ypos > 0:
        if seq1[xpos-1] == seq2[ypos-1]:
            traceback = seq1[xpos-1] + traceback
            xpos -= 1
            ypos -= 1
        elif dp[xpos-1][ypos] > dp[xpos][ypos-1]:
            xpos -= 1
        else:
            ypos -= 1
    return traceback

def sharedSplicedMotif(seq1,seq2):
    return longestCommonSubsequence(seq1,seq2)

def main():
    # seqs = readSequencesFromFasta("rosalind_lcsq.txt")
    # print(longestCommonSubsequence(seqs[">Rosalind_5877"],seqs[">Rosalind_5630"]))
    print(sharedSplicedMotif("AACCTTGG", "ACACTGTGA"))

if __name__ == "__main__":
    main()