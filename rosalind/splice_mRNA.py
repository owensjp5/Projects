################
# RNA Splicing #
################
from .sequence_alignment import alignSubstring

def spliceExons(seq, introns):
    for intron in introns:
        start = alignSubstring(seq, intron)[0]
        end = start+len(intron)
        seq = seq[:start-1] + seq[end-1:]
    return seq