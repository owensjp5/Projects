# Gillespie Stochastic Model
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random

def SSA():
    X = [0]
    t = [0]
    t_end = 1000
    k = 2
    gamma_X = 0.1
    print("Steady State:", k/gamma_X)

    while t[-1] < t_end:
        cur_X = X[-1]

        r = [k, gamma_X*cur_X]
        r_sum = sum(r)

        # Randomly select time increment, tau
        tau = np.random.exponential(scale=1/r_sum)
        t.append(t[-1] + tau)

        # Randomly select event
        rand = random.uniform(0,1)

        # Production event
        if rand * r_sum > 0 and rand * r_sum < r[0]:
            X.append(cur_X + 1)

        # Degredation event
        elif rand * r_sum > r[0] and rand * r_sum < r_sum:
            X.append(cur_X - 1)

    return t,X

def main():
    t,X = SSA()
    plt.plot(t,X)
    plt.xlabel("Time")
    plt.ylabel("Quantity mRNA")
    plt.show()

if __name__ == "__main__":
    main()