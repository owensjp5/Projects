###########################
# Finding a Protein Motif #
###########################
import subprocess
import re

def getSeqsFromUnitProt(accessIds):
    seqs = {}
    for id in accessIds:
        accession = id.split("_")[0]
        path = "https://rest.uniprot.org/uniprotkb/"  + accession + ".fasta"
        cmd = ["curl", "-L", path]
        rawSeq = subprocess.run(cmd, capture_output=True, text=True)
        seq = "".join(rawSeq.stdout.splitlines()[1:])
        seqs[id] = seq
    return seqs

def findMotif(accessIds, regex=r"N(?=[^P][ST][^P])"):
    seqs = getSeqsFromUnitProt(accessIds)
    for seqName in seqs:
        loci = re.finditer(regex, seqs[seqName])
        print(seqName)
        for loc in loci:
            print(loc.start()+1, end=" ")
        print()