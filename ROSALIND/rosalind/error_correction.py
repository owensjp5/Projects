#############################
# Error Correction in Reads #
#############################
from .read_fasta import readSequencesFromFasta
from .reverse_complement import reverseComplement
from .hamming_distance import getHammingDistance

def correctPointMutations(fastaFile):
    erroneousReads = []
    correctReads = []
    rawReads = readSequencesFromFasta(fastaFile)
    reads = []
    revComps = []
    for read in rawReads:
        reads.append(rawReads[read])
        revComps.append(reverseComplement(rawReads[read]))
    n = len(reads)
    for i in range(n):
        found = False
        for j in range(i):
            if reads[j] == reads[i] or revComps[j] == reads[i] or reads[j] == revComps[i]:
                # print(f"\033[1m{reads[j]}\033[0m", end=" ")
                found = True
                if reads[i] not in [r[0] for r in correctReads]: correctReads.append((reads[i], i))
                break
            # else:
            #     print(reads[j], end=" ")
        # print(f"\033[1m{reads[i]}\033[0m", end=" ")
        if found:
            # print()
            found = False
            continue
        for j in range(i+1,n):
            if reads[j] == reads[i] or revComps[j] == reads[i] or reads[j] == revComps[i]:
                # print(f"\033[1m{reads[j]}\033[0m", end=" ")
                found = True
                if reads[i] not in [r[0] for r in correctReads]: correctReads.append((reads[i], i))
                break
            # else:
            #     print(reads[j], end=" ")
        
        if not found:
            erroneousReads.append((reads[i], i))
    
    for erroneousReadTuple in erroneousReads:
        erroneousRead = erroneousReadTuple[0]
        i = erroneousReadTuple[1]
        for readTuple in correctReads:
            correctRead = readTuple[0]
            j = readTuple[1]
            if i != j:
                if getHammingDistance(erroneousRead, correctRead) == 1:
                    print(erroneousRead + "->" + correctRead)
                    break
                elif getHammingDistance(erroneousRead, revComps[j]) == 1:
                    print(erroneousRead + "->" + revComps[j])
                    break

def main():
    correctPointMutations("reads.fasta")

if __name__ == "__main__":
    main()