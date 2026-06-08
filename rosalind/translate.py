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