def readListInput(filePath):
    with open(filePath) as file:
        listStr = ""
        for line in file:
            listStr += line.rstrip()
    list = [float(x) for x in listStr.split()[1:]]
    return list