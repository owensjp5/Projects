# Central Dogma ODE Model w/ Euler's Method instead of 
import matplotlib.pyplot as plt

def eulerOneVar():
    m_init = 0

    t_end = 200

    k = 0.2
    gamma = 0.05

    M = [m_init]
    t = [0]

    delta_t = 0.5

    while t[-1] < t_end:
        delta_M = (k-gamma*M[-1]) * delta_t
        next_M = M[-1] + delta_M

        M.append(next_M)

        next_t = t[-1] + delta_t

        t.append(next_t)

    return t, M

def eulerTwoVars():
    m_init = 0
    p_init = 0

    t_end = 200

    k_m = 0.4
    gamma_m = 0.1
    k_p = 0.2
    gamma_p = 0.05

    M = [m_init]
    P = [p_init]
    t = [0]

    delta_t = 0.5

    while t[-1] < t_end:
        delta_M = (k_m - gamma_m*M[-1]) * delta_t    
        delta_P = (k_p*M[-1] - gamma_p*P[-1]) * delta_t
        
        next_M = M[-1] + delta_M
        next_P = P[-1] + delta_P
        next_t = t[-1] + delta_t

        M.append(next_M)
        P.append(next_P)
        t.append(next_t)

    return t, M, P

def main():
    t, M = eulerOneVar()

    f,ax = plt.subplots(1)
    line1, = ax.plot(t,M, color='b', label="Euler")
    ax.set_ylabel("Abundance")
    ax.set_xlabel("Time")
    ax.legend(handles=[line1])
    plt.show()

    t, M, P = eulerTwoVars()

    f,ax = plt.subplots(1)
    line1, = ax.plot(t,M, color="b",label="Euler M")
    line2, = ax.plot(t,P, color="r",label="Euler P")
    
    ax.set_ylabel("Quantity")
    ax.set_xlabel("Time(s)")
    ax.legend(handles=[line1,line2])
    plt.show()

if __name__ == "__main__":
    main()