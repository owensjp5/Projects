##################
# Overlap Graphs #
##################
from .read_fasta import readSequencesFromFasta
import numpy as np

def overlap(head, tail):
    len1 = len(head)
    len2 = len(tail)
    if len1 == len2:
        size = len1
        diff = 0
    else:
        diff = len1-len2
        if diff > 0:
            size = len2
        else:
            size = len1
            diff = 0
    
    buffer = np.zeros(shape=(size+1,size+1))

    for i in range(1,size+1):
        for j in range(1,size+1):
            if buffer[i-1][j-1]+1 == j:
                if head[diff+i-1] == tail[j-1]:
                    buffer[i][j] = j

    if buffer[size][size] == size:
        # Strings Match or One String is Prefix/Suffix of the Other
        return 0

    return max(buffer[size])

def overlapZ(head, tail):
    # Z algorithm to speed up overlap graph construction + assembly
    s = tail + "$" + head
    n = len(s)
    z = [0] * n
    l, r = 0, 0
    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]

    # The overlap is the z-value at the position where head starts,
    # but only if it extends to the end of head
    head_start = len(tail) + 1
    for i in range(head_start, n):
        if z[i] == n - i:  # z-value reaches end of string
            return n - i   # length of the matching suffix/prefix

    return 0

def constructOverlapGraph(fastaFile, k):
    sequences = readSequencesFromFasta(fastaFile)
    num_sequences = len(sequences)
    sequence_names = list(sequences.keys())

    adjacency_list = []

    for i in range(num_sequences):
        for j in range(num_sequences):
            if i!=j:
                overlap_value = overlapZ(sequences[sequence_names[i]],sequences[sequence_names[j]])
                if overlap_value == k:
                    adjacency = sequence_names[i] + " " + sequence_names[j] + " " + str(overlap_value)
                    adjacency_list.append(adjacency)
                    print(sequence_names[i][1:] + " " + sequence_names[j][1:])
    
    return adjacency_list