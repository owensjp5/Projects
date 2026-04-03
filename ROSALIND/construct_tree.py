#####################
# Completing a Tree #
#####################
def getAdjacencyList(inputFilePath):
    adjacencyList = []
    with open(inputFilePath) as inputfile:
        n = int(inputfile.readline())
        for line in inputfile:
            adjacency = [int(x) for x in line.rstrip().split()] #Split as ints instead of strs
            adjacencyList.append(adjacency)
    return n, adjacencyList

def connectEdges(inputFile):
    n, adjacencyList = getAdjacencyList(inputFile)
    branches = []
    visited = set()
    for adjacency in adjacencyList:
        commonGroup = -1
        first = adjacency[0]
        second = adjacency[1]
        visited.add(first)
        visited.add(second)
        for i in range(len(branches)):
            if first in branches[i]:
                # Add to branch
                branches[i].append(second)
                commonGroup = i
                break
        for i in range(len(branches)):
            if second in branches[i]:
                if commonGroup != -1:
                    if i != commonGroup:
                        # Merge branches
                        branches[commonGroup].pop()
                        branches[i] += branches[commonGroup]
                        branches.pop(commonGroup)
                        break
                else:
                    # Add to branch
                    branches[i].append(first)
                    commonGroup = True
                    break
        if commonGroup == -1:
            # Create new branch
            branches.append([first, second])
        commonGroup = -1
    unvisited = [[x] for x in (set(range(1,n+1)) - visited)]
    return branches + unvisited