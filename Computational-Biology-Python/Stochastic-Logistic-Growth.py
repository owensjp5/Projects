# Gillespie's Algorithm (SSA) Logistic Growth Model
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random

def SSA():
    X = [1]
    t = [0]

    t_end = 1000

    r = 0.05
    K = 1000

    while t[-1] < t_end:
        current_X = X[-1]

        rates = [r*current_X, (r*current_X**2)/K]

        r_sum = sum(rates)

        tau = np.random.exponential(scale=1/r_sum)
        t.append(t[-1] + tau)

        rand = random.uniform(0,1)

        if rand * r_sum <= rates[0]:
            X.append(current_X+1)
        else:
            X.append(current_X-1)

    return t, X

def main():
    t, x = SSA()

    plt.plot(t,x)
    plt.xlabel("Time(s)")
    plt.ylabel("Cells")
    plt.show()


if __name__ == "__main__":
    main()