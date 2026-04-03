#######################
# Open Reading Frames #
#######################
from .reverse_complement import reverseComplement

def ORF(sequence):
    codons = {
        "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
        "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
        "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
        "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
        "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
        "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
        "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "CAA": "Q", "AAA": "K", "GAA": "E",
        "CAG": "Q", "AAG": "K", "GAG": "E",
        "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
        "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "CGA": "R", "AGA": "R", "GGA": "G",
        "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",'':''
    }
    answers = []
    ORF1, ORF2, ORF3 = "", "", ""
    isORF1, isORF2, isORF3 = False, False, False
    seqLen = len(sequence)
    for i in range(3,seqLen+2,3):
        RF1 = sequence[i-3:i]
        RF2 = sequence[:seqLen-2][i-2:i+1]
        RF3 = sequence[:seqLen-1][i-1:i+2]
        if isORF1:
            if RF1 in ["TGA","TAG","TAA"]:
                answers.append(ORF1)
                ORF1 = ""
                isORF1 = False
            else:
                ORF1 += codons[RF1]
        elif RF1 == "ATG":
            isORF1 = True
            ORF1 += "M"
        if isORF2:
            if RF2 in ["TGA","TAG","TAA"]:
                answers.append(ORF2)
                ORF2 = ""
                isORF2 = False
            else:
                ORF2 += codons[RF2]
        elif RF2 == "ATG":
            isORF2 = True
            ORF2 += "M"
        if isORF3:
            if RF3 in ["TGA","TAG","TAA"]:
                answers.append(ORF3)
                ORF3 = ""
                isORF3 = False
            else:
                ORF3 += codons[RF3]
        elif RF3 == "ATG":
            isORF3 = True
            ORF3 += "M"
    return answers

def allORFs(sequence):
    answers = ORF(sequence) + ORF(reverseComplement(sequence))
    uniqueAnswers = []
    for answer in answers:
        for i in range(len(answer)):
            if answer[i] == 'M':
                if answer[i:] not in uniqueAnswers:
                    uniqueAnswers.append(answer[i:])
    return uniqueAnswers