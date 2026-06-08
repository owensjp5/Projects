##########################
# Finding a Motif in DNA #
##########################
def initializeSkipTable(substring, substringLength, skips, increment=0):
    if increment == 0:
        for i in range(substringLength):
            skips["A"][i] = substringLength
            skips["C"][i] = substringLength
            skips["G"][i] = substringLength
            skips["T"][i] = substringLength
            skips["U"][i] = substringLength
    for i in range(substringLength):
        skips[substring[i]][increment] = substringLength-(i+1)
    if substringLength > 1:
        return initializeSkipTable(substring[:-1], substringLength-1, skips, increment+1)
    else:
        return skips

def alignSubstring(string, substring, Display=False):
    length = len(string)
    substringLength = len(substring)

    # Preprocess Substring to use Bad Character Heuristic from Boyer-Moore
    skipTable = {"A":{},"C":{},"G":{},"T":{},"U":{}}
    skipTable = initializeSkipTable(substring, substringLength, skipTable)

    indices = []

    pos = substringLength-1

    while pos < length:
        for i in range(substringLength):
            curNucleotide = string[pos-i]
            skip = skipTable[curNucleotide][i]
            pos += skip
            if skip > 0:
                break
            if i == substringLength-1:
                indices.append(pos-substringLength+2)
                pos+=1
    if Display:
        string = ""
        for index in indices:
            string += str(index) + ' '
        print(string)
    else:
        return indices