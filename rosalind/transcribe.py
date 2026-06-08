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