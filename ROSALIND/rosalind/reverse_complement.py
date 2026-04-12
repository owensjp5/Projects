#################################
# Complementing a Strand of DNA #
#################################
def reverseComplement(sequence, mode="DNA"):
    if mode == "DNA":
        complementBase = {"A":"T","T":"A","C":"G","G":"C"}
    elif mode == "RNA":
        complementBase = {"A":"U","U":"A","C":"G","G":"C"}
    reverseComplement = ""
    for nucleotide in sequence:
        reverseComplement = complementBase[nucleotide] + reverseComplement
    return reverseComplement