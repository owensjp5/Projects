# 3 Gene Oscillating Network ODE Model
# Gene 1 activates Gene 2
# Gene 2 activates Gene 3
# Gene 3 represses Gene 1
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def simulate(variables, t, params):
    G1, G2, G3 = variables[0], variables[1], variables[2]
    k_1, gamma_1, k_2, gamma_2, k_3, gamma_3, n, c = params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7]
    
    dG1dt = (c**n / (c**n+G3**n)) * k_1 - (gamma_1 * G1)
    dG2dt = (G1**n / (c**n+G1**n)) * k_2 - (gamma_2 * G2)
    dG3dt = (G2**n / (c**n+G2**n)) * k_3 - (gamma_3 * G3)
    
    return ([dG1dt,dG2dt,dG3dt])

def main():
    y0 = [0,0,0]
    t = np.linspace(0,200,num=100)
    k_1 = 0.5
    gamma_1 = 0.1
    k_2 = 0.5
    gamma_2 = 0.1
    k_3 = 0.5
    gamma_3 = 0.1
    n = 9
    c = 1

    params = [k_1, gamma_1, k_2, gamma_2, k_3, gamma_3, n, c]

    y = odeint(simulate, y0, t, args=(params,))

    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
    line1, = ax1.plot(t,y[:,0], color='b',label="Gene 1")
    line2, = ax2.plot(t,y[:,1], color='r',label="Gene 2")
    line3, = ax3.plot(t,y[:,2], color='g',label="Gene 3")

    ax1.set_ylabel("Quantity")
    ax1.set_xlabel("Time")

    ax1.legend(handles=[line1,line2,line3])

    plt.show()

if __name__ == "__main__":
    main()