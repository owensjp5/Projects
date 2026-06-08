############################
# Counting DNA Nucleotides #
############################
def countNucleotides(sequence):
    nucleotide_counts = {"A":0,"C":0,"G":0,"T":0}
    for nucleotide in sequence:
        nucleotide_counts[nucleotide] += 1
    return nucleotide_counts["A"], nucleotide_counts["C"], nucleotide_counts["G"], nucleotide_counts["T"]
