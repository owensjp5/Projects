# Stochastic Simulation Involving
# Reaction Network w/ Chemicals A,B,C,D
# Rxn 1: A + B --> 2C
# Rxn 2: 2B + C --> D
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random

def SSA():
    A = [100]
    B = [100]
    C = [100]
    D = [100]

    t = [0]
    t_end = 10

    k1 = 0.3
    k2 = 0.01

    while t[-1] < t_end:
        propensities = [k1 * A[-1] * B[-1], k2 * B[-1]**2 * C[-1]]
        
        prop_sum = sum(propensities)

        if prop_sum == 0:
            break

        tau = np.random.exponential(scale=1/prop_sum)

        t.append(t[-1] + tau)

        rand = random.uniform(0,1)

        if rand * prop_sum <= propensities[0]:
            A.append(A[-1] - 1)
            B.append(B[-1] - 1)
            C.append(C[-1] + 2)
            D.append(D[-1])
        elif rand * prop_sum > propensities[0] and rand * prop_sum <= propensities[0] + propensities[1]:
            A.append(A[-1])
            B.append(B[-1] - 2)
            C.append(C[-1] - 1)
            D.append(D[-1] + 1)

    return t, A, B, C, D

def main():
    t, A, B, C, D = SSA()

    line1, = plt.plot(t, A, label="A")
    line2, = plt.plot(t, B, label="B")
    line3, = plt.plot(t, C, label="C")
    line4, = plt.plot(t, D, label="D")

    plt.legend(handles=[line1,line2,line3,line4])

    plt.xlabel("Time(s)")
    plt.ylabel("Abundance")
    plt.show()

if __name__ == "__main__":
    main()