# Logistic Growth ODE for Modeling Cell Culture Population Growth
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

def simulate(variables, t, params):    
    X = variables[0]
    
    r = params[0]
    K = params[1]

    dxdt = r*X*(1-(X/K))

    return([dxdt])

def main():
    x0 = 1
    
    t = np.linspace(0,500,num=1000)

    r = 0.05
    K = 1000

    params = [r,K]

    x = odeint(simulate, [x0], t, args=(params,))

    plt.plot(t,x[:,0])
    plt.xlabel("Time")
    plt.xlabel("Cells")
    plt.show()


if __name__ == "__main__":
    main()