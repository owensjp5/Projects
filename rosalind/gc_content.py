########################
# Computing GC Content #
########################
from .read_fasta import readSequencesFromFasta

def getGCContent(sequence):
    GCQuantity = 0
    for i in range(len(sequence)):
        if sequence[i] == "G" or sequence[i] == "C":
            GCQuantity += 1
    return GCQuantity/len(sequence) * 100

def computeGCContent(fastaFile):
    sequences = readSequencesFromFasta(fastaFile)
    highestGCContent = 0
    for sequence in sequences:
        GCContent = getGCContent(sequences[sequence])
        if GCContent > highestGCContent:
            highestGCContentRead = sequence
            highestGCContent = GCContent
    return highestGCContentRead, round(highestGCContent, 6)