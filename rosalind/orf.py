#######################
# Open Reading Frames #
#######################
from .reverse_complement import reverseComplement

DNA_codons = {
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
RNA_codons = {
        "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
        "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
        "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
        "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
        "UCU": "S", "CCU": "P", "ACU": "U", "GCU": "A",
        "UCC": "S", "CCC": "P", "ACC": "U", "GCC": "A",
        "UCA": "S", "CCA": "P", "ACA": "U", "GCA": "A",
        "UCG": "S", "CCG": "P", "ACG": "U", "GCG": "A",
        "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
        "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "CAA": "Q", "AAA": "K", "GAA": "E",
        "CAG": "Q", "AAG": "K", "GAG": "E",
        "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
        "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "CGA": "R", "AGA": "R", "GGA": "G",
        "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",'':''
    }

def ORF(sequence, mode="DNA"):
    if mode == "DNA":
        codons = DNA_codons
        start_codon = "ATG"
        stop_codons = ["TGA","TAG","TAA"]
    elif mode == "RNA":
        codons = RNA_codons
        start_codon = "AUG"
        stop_codons = ["UGA","UAG","UAA"]
    answers = []
    ORF1, ORF2, ORF3 = "", "", ""
    isORF1, isORF2, isORF3 = False, False, False
    seqLen = len(sequence)
    for i in range(3,seqLen+2,3):
        RF1 = sequence[i-3:i]
        RF2 = sequence[:seqLen-2][i-2:i+1]
        RF3 = sequence[:seqLen-1][i-1:i+2]
        if isORF1:
            if len(RF1) != 3: #If codon is < 3 chars, seq is cut off before a stop codon
                ORF1 = ""
                isORF1 = False
            elif RF1 in stop_codons:
                answers.append(ORF1)
                ORF1 = ""
                isORF1 = False
            else:
                ORF1 += codons[RF1]
        elif RF1 == start_codon:
            isORF1 = True
            ORF1 += "M"
        if isORF2:
            if len(RF2) != 3:
                ORF2 = ""
                isORF2 = False
            elif RF2 in stop_codons:
                answers.append(ORF2)
                ORF2 = ""
                isORF2 = False
            else:
                ORF2 += codons[RF2]
        elif RF2 == start_codon:
            isORF2 = True
            ORF2 += "M"
        if isORF3:
            if len(RF3) != 3:
                ORF3 = ""
                isORF3 = False
            elif RF3 in stop_codons:
                answers.append(ORF3)
                ORF3 = ""
                isORF3 = False
            else:
                ORF3 += codons[RF3]
        elif RF3 == start_codon:
            isORF3 = True
            ORF3 += "M"
    return answers

def allORFs(sequence, mode="DNA"):
    answers = ORF(sequence, mode) + ORF(reverseComplement(sequence, mode), mode)
    uniqueAnswers = []
    for answer in answers:
        for i in range(len(answer)):
            if answer[i] == 'M':
                if answer[i:] not in uniqueAnswers:
                    uniqueAnswers.append(answer[i:])
    return uniqueAnswers