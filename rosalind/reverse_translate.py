###############################
# Inferring mRNA from Protein #
###############################
def potentialmRNAs(protein):
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
    potentialStrings = 3 #Initialize as 3 to account for all 3 stop codons
    for AA in protein:
        potentialStrings *= codonCounts[AA]
        potentialStrings = potentialStrings % 1000000
    return potentialStrings