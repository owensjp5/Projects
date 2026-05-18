###################################
# Counting Phylogenetic Ancestors #
###################################
import math

def unrootedBinaryTreeRecursive(n):
    interNodes = 0
    if math.log2(n) % 1 == 0:
        temp = n
        while temp != 2:
            temp = int(temp / 2)
            interNodes += temp
        return interNodes
    else:
        m = (2 ** (math.log2(n) // 1))
        if n - m < 3:
            return unrootedBinaryTreeRecursive(m) + (n-m)
        else:
            return unrootedBinaryTreeRecursive(m) + unrootedBinaryTreeRecursive(n-m) + 2

def unrootedBinaryTree(n):
    return n - 2

def main():
    print(unrootedBinaryTreeRecursive(4))

if __name__ == "__main__":
    main()