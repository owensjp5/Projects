def readSequencesFromFasta(fastaFilePath, provideMaxLen=False):
    sequences = {}
    with open(fastaFilePath) as fastaFile:
        for line in fastaFile:
            if line[0] == ">":
                currentSequence = line.rstrip()
                sequences[currentSequence] = ""
            else:
                sequences[currentSequence] += line.rstrip()
    if provideMaxLen:
        max_len = 0
        longest_seq = ""
        for seq in sequences:
            cur_len = len(sequences[seq])
            if cur_len > max_len:
                longest_seq = seq
                max_len = cur_len
        return sequences, max_len, longest_seq
    return sequences