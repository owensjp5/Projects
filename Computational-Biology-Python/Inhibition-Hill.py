# Two Gene ODE Inhibition Model
# Gene 1 encodes for transcription factor that represses expression of Gene 2
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def simulate(variables, t, params):
    G1 = variables[0]
    G2 = variables[1]

    k_1 = params[0]
    gamma_1 = params[1]
    k_2 = params[2]
    gamma_2 = params[3]
    n = params[4]
    c = params[5]

    dG1dt = k_1 - (gamma_1 * G1)                            # Change in Gene 1 protein - function of Gene 1 expression rate(k_1), Gene 1 protein decay rate(gamma_1), and Gene 1 protein quantity(G1)
    dG2dt = (c**n/(c**n+G1**n)) * k_2 - (gamma_2 * G2)     # Change in Gene 2 protein - function of Hill function, Gene 2 expression rate(k_2), Gene 2 protein decay rate(gamma_2), and Gene 2 protein quantity(G2)

    return ([dG1dt, dG2dt])

def getSteadyState(params):
    k_1 = params[0]
    gamma_1 = params[1]
    k_2 = params[2]
    gamma_2 = params[3]
    n = params[4]
    c = params[5]

    G1_SS = k_1/gamma_1
    G2_SS = ((c**n/(c**n + G1_SS**n)) * k_2) / gamma_2
    
    return G1_SS, G2_SS

def main():
    #Define Parameters
    y0 = [0,0]
    t = np.linspace(0,200,num=100)

    k_1 = 0.4
    gamma_1 = 0.1
    k_2 = 0.9
    gamma_2 = 0.1
    n = 4
    c = 1

    params = [k_1,gamma_1,k_2,gamma_2,n,c]
    steadyStates = getSteadyState(params)
    print("Steady States:", steadyStates)

    #Run Simulation
    y = odeint(simulate, y0, t, args=(params,))

    #Plot Results
    f, ax = plt.subplots(1)
    line1, = ax.plot(t,y[:,0], color='b',label="Gene 1")
    line2, = ax.plot(t,y[:,1], color='r',label="Gene 2")

    ax.set_ylabel("Quantity")
    ax.set_xlabel("Time(s)")
    ax.legend(handles=[line1,line2])

    plt.show()

if __name__ == "__main__":
    main()