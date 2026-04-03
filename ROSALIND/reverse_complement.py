#################################
# Complementing a Strand of DNA #
#################################
def reverseComplement(sequence):
    complementBase = {"A":"T","T":"A","C":"G","G":"C"}
    reverseComplement = ""
    for nucletide in sequence:
        reverseComplement = complementBase[nucletide] + reverseComplement
    return reverseComplement