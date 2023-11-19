#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:37:42 2023

@author: peterm
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Constants
g_big_cgs = 6.67430e-8  # dyn·cm^2/g^2
m_sol_cgs = 1.98847e33  # g
c_cgs = 3e10  # cm/s
sigma_cgs = 5.670374419e-5  # erg·cm^(-2)·s^(-1)·K^(-4)
h = 6.62607015e-27  # erg·s
k = 1.380649e-16  # erg/K

# Accretion rate
m_dot = 1e18

# Define the Schwarzschild radius
m = 10*m_sol_cgs  # Replace with the desired mass
r_g = g_big_cgs * m / c_cgs**2

# Define the integration limits
r_in = 6
r_out = 1e5

# Define the range of nu values
nu_min = 1e14  
nu_max = 1e19  

# Define the temperature function
def T(rad):
    numerator = 3 * g_big_cgs * m_dot * 10 * m_sol_cgs
    denominator = 8 * np.pi * (rad)**3 * r_g**3 * sigma_cgs
    bracket = 1 - np.sqrt(r_in/(rad))
    return ((numerator/denominator)*bracket)**0.25


# Define the number of points
r_num_points = 100

# Create a radial grid with logarithmic spacing, excluding zero
r_grid = np.logspace(np.log10(r_in), np.log10(r_out), num=r_num_points, endpoint=False)

#print(r_grid)


# Calculate midpoints of the bins

r_midpoints = []
for i in range(0,len (r_grid)-1):
    m = (r_grid[i] + r_grid[i+1]) / 2.0
    r_midpoints.append(m)

#print(r_midpoints)

# Calculate temperatures for each midpoint using the T(rad) function
temperatures = [T(rad) for rad in r_midpoints]

log_r_midpoints = np.log10(r_midpoints)

plt.figure(dpi=600)
plt.scatter(log_r_midpoints, temperatures, s=1)
plt.xlabel('r')
plt.ylabel('Temperature (K)')
plt.grid(False)
plt.show()
 

#print(temperatures)

# Plot temperature against midpoints
plt.figure(dpi=600)
plt.semilogx(r_midpoints, temperatures, 0.1)
plt.xlabel('r')
plt.ylabel('Temperature (K)')
plt.grid(False)
plt.show()

# Define the flux function
def F(nu, temperature):
        numerator = (2 * np.pi * h * nu**3) / c_cgs**2
        denominator = np.exp((h * nu)/(k * temperature))-1
        return numerator/denominator
    
# Plotting a temporary nu value for the flux against r

temp_nu = 10e15

fluxes = [F(temp_nu,t) for t in temperatures]

plt.figure(dpi=600)
plt.semilogx(r_midpoints, fluxes)
plt.xlabel('r')
plt.ylabel('Flux (10e15 Hz) (Mx)')
plt.grid(False)
plt.show()


# Create log spaced frequency values

nu_grid = np.logspace(np.log10(nu_min), np.log10(nu_max), num=r_num_points, endpoint=False)

# Define the luminosity integral

L_values = []
nu_L_nu_values = []

for nu_value in nu_grid:
    def integrand(r):
        return F(nu_value, T(r)) * 4 * np.pi * r_g**2 * r

    result, error = quad(integrand, r_in, r_out)
    L_values.append(result)
    nu_L_nu_values.append(result*nu_value)
    
log_nu_L_nu = np.log10(nu_L_nu_values)

# Plotting luminosity against frequency


plt.figure(dpi=600)
plt.semilogx(nu_grid, nu_L_nu_values)
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\nu L_\nu$')
plt.grid(False)
plt.show()



plt.figure(dpi=600)
plt.semilogx(nu_grid, log_nu_L_nu)
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\nu L_\nu$')
plt.grid(False)
#plt.ylim(30, 40)
plt.show()

print(sum(L_values))