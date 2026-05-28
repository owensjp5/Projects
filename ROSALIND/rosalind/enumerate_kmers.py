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

########################################################
# Ordering Strings of Varying Length Lexicographically #
########################################################
def enumerateLexicographicalStrings(alphabet, N, current=""):
    if current != "":
        print(current) #Only difference w/ lexicographicalKmers is printing regardless of len
    if N != 0:
        for char in alphabet:
            cur = current + char
            enumerateLexicographicalStrings(alphabet, N-1, cur)

def main():
    enumerateLexicographicalStrings(['E', 'M', 'B', 'Q', 'G', 'U', 'O', 'D', 'K', 'F', 'A'], 3)

if __name__ == "__main__":
    main()