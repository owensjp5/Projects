########################################
# Enumerating k-mers Lexicographically #
########################################
def lexicographicalKmers(alphabet, N, current=""):
    if N == 0:
        print(current)
    else:
        for char in alphabet:
            cur = current + char
            lexicographicalKmers(alphabet, N-1, cur)