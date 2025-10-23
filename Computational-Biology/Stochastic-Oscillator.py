# Stochastic 3-Gene Goodwin Oscillator Model
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random

def SSA():
    G = [[0],[0],[0]]
    t = [0]

    t_end = 1000

    k = [2,2,2]
    gamma = [0.1,0.1,0.1]
    n = 9
    c = 1

    while t[-1] < t_end:
        current_G1 = G[0][-1]
        current_G2 = G[1][-1]
        current_G3 = G[2][-1]

        r = [(c**n         /(c**n + current_G3**n))*k[0], gamma[0]*current_G1, \
             (current_G1**n/(c**n + current_G1**n))*k[1], gamma[1]*current_G2, \
             (current_G2**n/(c**n + current_G2**n))*k[2], gamma[2]*current_G3]

        r_sum = sum(r)

        tau = np.random.exponential(scale=1/r_sum)
        t.append(t[-1] + tau)

        rand = random.uniform(0,1)

        if rand * r_sum < r[0]:
            G[0].append(current_G1 + 1)
            G[1].append(current_G2)
            G[2].append(current_G3)
        elif rand * r_sum > r[0] and rand * r_sum < sum(r[:2]):
            G[0].append(current_G1 - 1)
            G[1].append(current_G2)
            G[2].append(current_G3)
        elif rand * r_sum > sum(r[:2]) and rand * r_sum < sum(r[:3]):
            G[0].append(current_G1)
            G[1].append(current_G2 + 1)
            G[2].append(current_G3)
        elif rand * r_sum > sum(r[:3]) and rand * r_sum < sum(r[:4]):
            G[0].append(current_G1)
            G[1].append(current_G2 - 1)
            G[2].append(current_G3)
        elif rand * r_sum > sum(r[:4]) and rand * r_sum < sum(r[:5]):
            G[0].append(current_G1)
            G[1].append(current_G2)
            G[2].append(current_G3 + 1)
        elif rand * r_sum > sum(r[:5]) and rand * r_sum < r_sum:
            G[0].append(current_G1)
            G[1].append(current_G2)
            G[2].append(current_G3 - 1)

    return t, G

def main():
    t, G = SSA()
    
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
    line1, = ax1.plot(t, G[0], color="b", label="Gene 1")
    line2, = ax2.plot(t, G[1], color="r", label="Gene 2")
    line3, = ax3.plot(t, G[2], color="g", label="Gene 3")
    ax1.set_ylabel("Quantity")
    ax1.set_xlabel("Time")
    ax1.legend(handles=[line1, line2, line3])
    plt.show()

if __name__ == "__main__":
    main()