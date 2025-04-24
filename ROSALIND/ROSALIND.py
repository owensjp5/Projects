import numpy as np

############################
# Counting DNA Nucleotides #
############################
def countNucleotides(sequence):
    nucleotide_counts = {"A":0,"C":0,"G":0,"T":0}
    for nucleotide in sequence:
        nucleotide_counts[nucleotide] += 1
    return nucleotide_counts["A"], nucleotide_counts["C"], nucleotide_counts["G"], nucleotide_counts["T"]

#############################
# Transcribing DNA into RNA #
#############################
def transcribe(sequence):
    RNA = ""
    for nucleotide in sequence:
        if nucleotide == "T":
            RNA += "U"
        else:
            RNA += nucleotide
    return RNA

#################################
# Complementing a Strand of DNA #
#################################
def reverseComplement(sequence):
    complementBase = {"A":"T","T":"A","C":"G","G":"C"}
    reverseComplement = ""
    for nucletide in sequence:
        reverseComplement = complementBase[nucletide] + reverseComplement
    return reverseComplement

####################################
# Rabbits and Recurrence Relations #
####################################
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

########################
# Computing GC Content #
########################
def getGCContent(sequence):
    GCQuantity = 0
    for i in range(len(sequence)):
        if sequence[i] == "G" or sequence[i] == "C":
            GCQuantity += 1
    return GCQuantity/len(sequence) * 100

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

############################
# Counting Point Mutations #
############################
def getHammingDistance(seq1, seq2):
    length1 = len(seq1)
    length2 = len(seq2)
    if length1 > 1000 or length2 > 1000:
        print("Error: improper input for getHammingDistance fn, input contains sequence longer than 1 kbp")
    if length1 != length2:
        print("Error: improper input for getHammingDistance fn, sequences are different length")
        return
    hammingDistance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hammingDistance += 1
    return hammingDistance

######################
# Mendel's First Law #
######################
def mendelianInheritance(k, m, n):
    total = k + m + n
    probDom1 = k/total + (m/total * 0.5)
    probRec1Dom2 = ((n/total) * (k/(total-1) + m/(total-1)*0.5)) + ((m/total*0.5) * (k/(total-1) + ((m-1)/(total-1)*0.5)))
    probDom = probDom1 + probRec1Dom2
    # probRec = (n/total * ((n-1)/(total-1) + (m/(total-1)*0.5))) + ((m/total*0.5) * ((n/(total-1) + (m-1)/(total-1)*0.5)))
    # print("Probabilities should sum to 1, Sum:", probDom+probRec)
    return probDom

################################
# Translating RNA into Protein #
################################
def translate(sequence):
    codons = {
    "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
    "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
    "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
    "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
    "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
    "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
    "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
    "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
    "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
    "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
    "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
    "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
    if sequence[:3] != "AUG":
        print("Error: improper input for translate fn, mRNA sequence does not begin with start codon (AUG)")
        return
    try:
        if codons[sequence [-3:]] != "Stop":
            print("Error: improper input for translate fn, mRNA sequence does not end with stop codon (UAA, UAG, UGA)")
            return
    except:
        print("Error: improper input for translate fn, mRNA contains invalid nucleotide(can only be As, Cs, Us, or Gs)")
        return
    protein = ""
    for i in range(0, len(sequence)-3, 3):
        try:
            protein += codons[sequence[i:i+3]]
        except:
            print("Error: improper input for translate fn, mRNA contains invalid nucleotide(can only be As, Cs, Us, or Gs)")
            return
    return protein

def main():
    # print(countNucleotides("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
    # print(transcribeDNA("GATGGAACTTGACTACGTAAATT"))
    # print(reverseComplement("AAAACCCGGT"))
    # print(recurrenceRelationDP(5, 3))
    # print(computeGCContent("reads.fasta"))
    # print(getHammingDistance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
    # print(mendelianInheritance(2,2,2))
    # print(translate("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
    print(substring("GATATATGCATATACTT", "ATAT"))

if __name__ == "__main__":
    main()