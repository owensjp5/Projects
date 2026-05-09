################################
# Transitions and Translations #
################################

def mutationRatio(seq1, seq2):
    purines = ["A","G"]
    pyrimidines = ["C","T"]
    transitions = 0
    translations = 0

    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            continue
        if seq1[i] in purines:
            if seq2[i] in purines:
                transitions += 1
            else:
                translations += 1
        else:
            if seq2[i] in purines:
                translations += 1
            else:
                transitions += 1

    return transitions / translations

def main():
    print(mutationRatio("GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT", "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"))

if __name__ == "__main__":
    main()