# Lotka-Voltera Predator-Prey Model
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def simulate(variables, t, params):
    x = variables[0]
    y = variables[1]

    alpha = params[0]
    beta = params[1]
    delta = params[2]
    gamma = params[3]

    dxdt = alpha*x - beta*x*y
    dydt = delta*x*y - gamma*y

    return([dxdt, dydt])

def main():
    y0 = [10,1] #[fish, bears]
    
    t = np.linspace(0,50,num=1000)

    alpha = 1.1
    beta = 0.4
    delta = 0.1
    gamma = 0.4

    params = [alpha, beta, delta, gamma]

    y = odeint(simulate, y0, t, args=(params,))

    f,(ax1,ax2) = plt.subplots(2)

    line1, = ax1.plot(t,y[:,0], color="b")
    line2, = ax2.plot(t,y[:,1], color="r")

    ax1.set_ylabel("Fish (hundreds)")
    ax2.set_ylabel("Bears (hundreds)")
    ax2.set_xlabel("Time")

    plt.show()

if __name__ == "__main__":
    main()