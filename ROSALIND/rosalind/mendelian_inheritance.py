######################
# Mendel's First Law #
######################
def dominantOffspringProbability(k, m, n):
    total = k + m + n
    probDom1 = k/total + (m/total * 0.5)
    probRec1Dom2 = ((n/total) * (k/(total-1) + m/(total-1)*0.5)) + ((m/total*0.5) * (k/(total-1) + ((m-1)/(total-1)*0.5)))
    probDom = probDom1 + probRec1Dom2
    # probRec = (n/total * ((n-1)/(total-1) + (m/(total-1)*0.5))) + ((m/total*0.5) * ((n/(total-1) + (m-1)/(total-1)*0.5)))
    # print("Probabilities should sum to 1, Sum:", probDom+probRec)
    return probDom