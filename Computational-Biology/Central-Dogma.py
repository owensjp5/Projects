# Central Dogma ODE Model
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def simulate(variables, t, params):
    m = variables[0]
    p = variables[1]
    
    k_m = params[0]
    gamma_m = params[1]
    k_p = params[2]
    gamma_p = params[3]

    dmdt = k_m - (gamma_m * m)          # Change in mRNA - function of transcription rate(k_m), mRNA decay rate(gamma_m), and mRNA quantity(m)
    dpdt = (k_p * m) - (gamma_p * p)    # Change in protein - function of translation rate(k_p), mRNA quantity(m), protein decay rate(gamma_p), and protein quantity(p)

    return ([dmdt, dpdt])
    
def getSteadyState(params):
    k_m = params[0]
    gamma_m = params[1]
    k_p = params[2]
    gamma_p = params[3]

    m_SS = k_m/gamma_m
    p_SS = (k_p * m_SS) / gamma_p
    
    return m_SS, p_SS

def main():
   #Initialize Variables
    y0 = [0,0]
    t = np.linspace(0,200,100)

    k_m = 0.2
    gamma_m = 0.05
    k_p = 0.4
    gamma_p = 0.2

    rates = [k_m, gamma_m, k_p, gamma_p]

    print("Steady States:", getSteadyState(rates))

    #Run Simulation
    y = odeint(simulate, y0, t, args=(rates,))
    
    #Plot Results
    f, ax = plt.subplots(1)
    line1, = ax.plot(t,y[:,0], color='b',label="mRNA")
    line2, = ax.plot(t,y[:,1], color='r',label="Protein")

    ax.set_ylabel("Quantity")
    ax.set_xlabel("Time")
    ax.legend(handles=[line1,line2])

    plt.show()

if __name__ == "__main__":
    main()