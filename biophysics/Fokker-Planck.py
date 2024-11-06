# Python program to simulate the stochastic movement of a molecular motor with the Fokker-Planck Equation
# Jack Owens, Nov 5 2024

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve1d

# Define Parameters
v = 3                               #Drift velocity
D = 0.2                             #Diffusion constant
dx = 1.0                            #Space interval
dt = 0.1                            #Time interval
minimum_x, maximum_x = -50, 300    #Spacial domain
total_time = 150                    #Total simulation runtime

# Initialize Arrays
x_domain = np.arange(minimum_x, maximum_x + dx, dx) #Initialize Domain
prob_dist = np.zeros_like(x_domain)                 #Initialize probability distribution
dist_center = -1 * minimum_x                        #Initialize starting position probability as 1
prob_dist[dist_center] = 1.0 / dx 

# Calculate Finite Difference Constants
drift_constant = v * dt / (2 * dx)  #Asymmetric finite difference approximation for first derivative
diffusion_constant = D * dt / dx**2 #Finite difference approximation for second derivative

# Kernels for Diffusion and Drift
diffusion_kernel = np.array([diffusion_constant, 1-2*diffusion_constant, diffusion_constant]) #Centered distribution
drift_kernel = np.array([0, 1-drift_constant, drift_constant]) #Asymmetric distribution
# drift_kernel = np.array([-drift_constant, 1, drift_constant])

plt.figure(figsize=(10,6))

# Run Simulation
epochs = int(total_time / dt)
for epoch in range(epochs):
    prob_dist = convolve1d(prob_dist,diffusion_kernel,mode='constant')
    prob_dist = convolve1d(prob_dist,drift_kernel,mode='constant')
    prob_dist = prob_dist / np.sum(prob_dist) #Normalize
    if (epoch+1) % 200 == 0:
        plt.plot(x_domain, prob_dist, label=f'{epoch*dt:.1f} ms')

# Plot Results
plt.plot(x_domain, prob_dist, label=f'{total_time} ms')
plt.xlabel("Position (nm)")
plt.ylabel("Probability")
plt.title("Probability Distributions of Molecular Motor Position")
plt.legend()
plt.show()