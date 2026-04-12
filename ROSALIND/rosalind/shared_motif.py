##########################
# Finding a Shared Motif #
##########################
from .read_fasta import readSequencesFromFasta

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
                    if seqs[longest_seq][j:i] not in seqs[seq]:
                        shared = False
                        break
            if shared:
                motifs[seqs[longest_seq][j:i]] = True
    
    motifs = sorted(list(motifs.items()), key = lambda item: len(item[0]), reverse=True)

    for motif in motifs:
        if motif[1]:
            return motif[0]
        
    return "No shared motifs"