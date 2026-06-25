########################################
# Introduction to Alternative Splicing #
########################################
def alternateSplices(n, m):
    if m > n:
        print("Error: m must be less than or equal to n")
        return
    elif n > 2000:
        print("Error: n must be less than 2000")
        return
    total = n + 1
    prevTerm = n
    for i in range(1,(n-m)):
        newTerm = prevTerm * (n-i)
        newTerm = newTerm // (i+1)
        prevTerm = newTerm
        total = (total + newTerm) % 1000000
    return total

def main():
    print(alternateSplices(6, 3))

if __name__ == "__main__":
    main()