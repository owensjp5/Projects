import numpy as np

##########################
# Counting DNA Nucleotides
##########################
def countNucleotides(DNA):
    nucleotide_counts = {"A":0,"C":0,"G":0,"T":0}
    for nucleotide in DNA:
        nucleotide_counts[nucleotide] += 1
    return nucleotide_counts["A"], nucleotide_counts["C"], nucleotide_counts["G"], nucleotide_counts["T"]

###########################
# Transcribing DNA into RNA
###########################
def transcribeDNA(DNA):
    RNA = ""
    for nucleotide in DNA:
        if nucleotide == "T":
            RNA += "U"
        else:
            RNA += nucleotide
    return RNA

##############################
#Complementing a Strand of DNA
##############################
def reverseComplement(DNA):
    complementBase = {"A":"T","T":"A","C":"G","G":"C"}
    reverseComplement = ""
    for nucletide in DNA:
        reverseComplement = complementBase[nucletide] + reverseComplement
    return reverseComplement

#################################
#Rabbits and Recurrence Relations
#################################
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

#####################
#Computing GC Content
#####################
def getGCContent(DNA):
    GCQuantity = 0
    for i in range(len(DNA)):
        if DNA[i] == "G" or DNA[i] == "C":
            GCQuantity += 1
    return GCQuantity/len(DNA) * 100

def readSequencesFromFasta(fastaFilePath):
    sequences = {}
    with open(fastaFilePath) as fastaFile:
        for line in fastaFile:
            if line[0] == ">":
                currentSequence = line
                sequences[currentSequence] = ""
            else:
                sequences[currentSequence] += line.rstrip()
    return sequences

def computeGCContent(fastaFile):
    sequences = readSequencesFromFasta(fastaFile)
    highestGCContent = 0
    for sequence in sequences:
        GCContent = getGCContent(sequences[sequence])
        if GCContent > highestGCContent:
            highestGCContentRead = sequence
            highestGCContent = GCContent
    return highestGCContentRead, highestGCContent


def main():
    print(countNucleotides("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
    print(transcribeDNA("GATGGAACTTGACTACGTAAATT"))
    print(reverseComplement("AAAACCCGGT"))
    print(recurrenceRelationDP(5, 3))
    print(computeGCContent("reads.fasta"))

if __name__ == "__main__":
    main()