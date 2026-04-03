##################################
# Introduction to Random Strings #
##################################
import math

def randomString(seq, gcContent):
        logProbability = 0
        for base in seq:
                if base == "C" or base == "G":
                        logProbability += math.log10(gcContent/2)
                else:
                        logProbability += math.log10((1-gcContent)/2)
        return logProbability
