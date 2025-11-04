# Solutions to Genomics Problems from ROSALIND
import numpy as np
import math

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

########################
# Computing GC Content #
########################
def getGCContent(sequence):
    GCQuantity = 0
    for i in range(len(sequence)):
        if sequence[i] == "G" or sequence[i] == "C":
            GCQuantity += 1
    return GCQuantity/len(sequence) * 100

def readSequencesFromFasta(fastaFilePath, provideMaxLen=False):
    sequences = {}
    with open(fastaFilePath) as fastaFile:
        for line in fastaFile:
            if line[0] == ">":
                currentSequence = line.rstrip()
                sequences[currentSequence] = ""
            else:
                sequences[currentSequence] += line.rstrip()
    if provideMaxLen:
        max_len = 0
        longest_seq = ""
        for seq in sequences:
            cur_len = len(sequences[seq])
            if cur_len > max_len:
                longest_seq = seq
                max_len = cur_len
        return sequences, max_len, longest_seq
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
        return
    if length1 != length2:
        print("Error: improper input for getHammingDistance fn, sequences are different length")
        return
    hammingDistance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hammingDistance += 1
    return hammingDistance

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
        "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }
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

##########################
# Finding a Motif in DNA #
##########################
def initializeSkipTable(substring, substringLength, skips, increment=0):
    if increment == 0:
        for i in range(substringLength):
            skips["A"][i] = substringLength
            skips["C"][i] = substringLength
            skips["G"][i] = substringLength
            skips["T"][i] = substringLength
            skips["U"][i] = substringLength
    for i in range(substringLength):
        skips[substring[i]][increment] = substringLength-(i+1)
    if substringLength > 1:
        return initializeSkipTable(substring[:-1], substringLength-1, skips, increment+1)
    else:
        return skips

def alignSubstring(string, substring):
    length = len(string)
    substringLength = len(substring)

    # Preprocess Substring to use Bad Character Heuristic from Boyer-Moore
    skipTable = {"A":{},"C":{},"G":{},"T":{},"U":{}}
    skipTable = initializeSkipTable(substring, substringLength, skipTable)

    indices = []

    pos = substringLength-1

    while pos < length:
        for i in range(substringLength):
            curNucleotide = string[pos-i]
            skip = skipTable[curNucleotide][i]
            pos += skip
            if skip > 0:
                break
            if i == substringLength-1:
                indices.append(pos-substringLength+2)
                pos+=1
    return indices

#########################
# Consensus and Profile #
#########################
def profileSequences(seqs):
    n = len(seqs[1])
    totalSeqs = len(seqs)

    profile = np.zeros((4,n))
    consensusString = ""

    for i in range(n):
        for j in range(totalSeqs):
            if seqs[j][i] == "A":
                profile[0][i] += 1
            if seqs[j][i] == "C":
                profile[1][i] += 1
            if seqs[j][i] == "G":
                profile[2][i] += 1
            if seqs[j][i] == "T":
                profile[3][i] += 1
        consensus = ""
        max = 0
        if profile[0][i] > max:
            max = profile[0][i]
            consensus = "A"
        if profile[1][i] > max:
            max = profile[1][i]
            consensus = "C"
        if profile[2][i] > max:
            max = profile[2][i]
            consensus = "G"
        if profile[3][i] > max:
            max = profile[3][i]
            consensus = "T"
        consensusString += consensus
    print(consensusString)
    return profile

##################
# Overlap Graphs #
##################
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

def constructOverlapGraph(fastaFile):
    sequences = readSequencesFromFasta(fastaFile)
    num_sequences = len(sequences)
    sequence_names = list(sequences.keys())

    adjacency_list = []

    for i in range(num_sequences):
        for j in range(num_sequences):
            if j != i:
                overlap_value = overlap(sequences[sequence_names[i]],sequences[sequence_names[j]])
                if overlap_value > 0:
                    adjacency = sequence_names[i] + " " + sequence_names[j] + " " + str(overlap_value)
                    adjacency_list.append(adjacency)
    
    return adjacency_list

################################
# Calculate Expected Offspring #
################################
def expectedOffspring(population):
    total = 0
    for i in range(len(population)):
        if i < 3:
            total += 2*population[i]
        elif i==4:
            total += population[i]
        elif i == 3:
            total += 1.5*population[i]
    return total

##########################
# Finding a Shared Motif #
##########################
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
                    if alignSubstring(seqs[seq], seqs[longest_seq][j:i]) == []:
                        shared = False
                        break
            if shared:
                motifs[seqs[longest_seq][j:i]] = True
    
    motifs = sorted(list(motifs.items()), key = lambda item: len(item[0]), reverse=True)

    for motif in motifs:
        if motif[1]:
            return motif[0]
        
    return "No shared motifs"

###############################
# Inferring mRNA from Protein #
###############################
def reverseTranslate(protein):
    codonCounts = {
        "A": 4,  # Alanine
        "R": 6,  # Arginine
        "N": 2,  # Asparagine
        "D": 2,  # Aspartic acid
        "C": 2,  # Cysteine
        "Q": 2,  # Glutamine
        "E": 2,  # Glutamic acid
        "G": 4,  # Glycine
        "H": 2,  # Histidine
        "I": 3,  # Isoleucine
        "L": 6,  # Leucine
        "K": 2,  # Lysine
        "M": 1,  # Methionine (also start)
        "F": 2,  # Phenylalanine
        "P": 4,  # Proline
        "S": 6,  # Serine
        "T": 4,  # Threonine
        "W": 1,  # Tryptophan
        "Y": 2,  # Tyrosine
        "V": 4,  # Valine
    }
    potentialStrings = math.log10(3) #Initialize as number of stop codons, sum logs to avoid computational constraints of large numbers
    for AA in protein:
        potentialStrings += math.log10(codonCounts[AA])
    return int(10**potentialStrings % 1000000)

def main():
    # print(countNucleotides("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
    # print(transcribeDNA("GATGGAACTTGACTACGTAAATT"))
    # print(reverseComplement("AAAACCCGGT"))
    # print(computeGCContent("reads.fasta"))
    # print(getHammingDistance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
    # print(translate("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
    # print(alignSubstring("GATATATGCATATACTT", "ATAT"))
    # print(profileSequences(["ATCCAGCT","GGGCAACT","ATGGATCT","AAGCAACC","TTGGAACT","ATGCCATT","ATGGCACT"]))
    # print(constructOverlapGraph("reads2.fasta"))
    # print(expectedOffspring([1, 0, 0, 1, 0, 1]))
    # print(sharedMotif("reads3.fasta"))
    print(reverseTranslate("CYIQNCPLG")) #Oxytocin Peptide

if __name__ == "__main__":
    main()