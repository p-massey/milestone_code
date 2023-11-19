#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:08:28 2023

@author: peterm
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:09:29 2023

@author: peterm
"""

import numpy as np
from scipy.integrate import simpson
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


# Define the number of points for radial and frequency

r_num_points = 2000
nu_num_points = 2000


# Create a radial grid with logarithmic spacing

r_grid = np.linspace(np.log10(r_in), np.log10(r_out), num = r_num_points)


# Calculate midpoints of the bins

r_midpoints = []
for i in range(0,len (r_grid)-1):
    m = 10**(r_grid[i] + ((r_grid[i+1] - r_grid[i]) / 2.0))
    r_midpoints.append(m)


# Calculate temperatures for each midpoint using the T(rad) function

temperatures = [T(rad) for rad in r_midpoints]


# Make a logarithmic set of r_midpoints

log_r_midpoints = np.log10(r_midpoints)

# Plot the temperature against the radius

plt.figure(dpi=1000)
plt.scatter(log_r_midpoints, temperatures, s=1)
plt.xlabel(r'$log_{10}(r)$')
plt.ylabel('Temperature (K)')
plt.grid(False)
plt.show()


# Define the flux function

def F(nu, r):
        numerator = (2 * np.pi * h * nu**3) / (c_cgs**2)
        denominator = np.exp((h * nu)/(k * T(r)))-1
        return numerator/denominator
    
    
print(F(1e15, 10))
    

# Plotting a temporary nu value for the flux against r

temp_nu = 1e15

fluxes = [F(temp_nu,r) for r in r_midpoints]

plt.figure(dpi=600)
plt.scatter(log_r_midpoints, fluxes, s=1)
plt.xlabel(r'$log_10(r)$')
plt.ylabel('Flux (1e15 Hz) (Mx)')
plt.grid(False)
plt.show()

# Plotting frequency against temperature

plt.figure(dpi=600)
plt.scatter(temperatures, fluxes, s=1)
plt.xlabel('temperatures')
plt.ylabel('fluxes')
plt.grid(False)
plt.show()


# Create a frequency grid with logarithmic spacing

nu_grid = np.linspace(np.log10(nu_min), np.log10(nu_max), num = nu_num_points)


# define the integrand for luminosity function

def integrand(nu,r):
    return F(nu, r) * 4 * np.pi * r_g**2 * r

# Plotting a temporary nu value for the integrand against radius

integrand_values = [integrand(temp_nu,i) for i in r_midpoints]

plt.figure(dpi=600)
plt.scatter(log_r_midpoints, integrand_values, s=1)
plt.xlabel('temperatures')
plt.ylabel('fluxes')
plt.grid(False)
plt.show()


# Integral for luminosity for a given frequency

def integral(nu):
    area =[]
    for i in range(0,len(r_midpoints)):
        rect = integrand(nu, r_midpoints[i])*(10**r_grid[i+1]-10**r_grid[i])
        area.append(rect)
    return sum(area)


# plot of luminosity against frequency

nu_midpoints = []
for i in range(0,len (nu_grid)-1):
    m = 10**(nu_grid[i] + ((nu_grid[i+1] - nu_grid[i]) / 2.0))
    nu_midpoints.append(m)
    
integral_values = [integral(i) for i in nu_midpoints]

# Make a logarithmic set of nu_midpoints

log_nu_midpoints = np.log10(nu_midpoints)

plt.figure(dpi=600)
plt.scatter(log_nu_midpoints, integral_values, s=1)
plt.xlabel('nu')
plt.ylabel('L')
plt.grid(False)
plt.show()

print(integral(temp_nu))

# Plot of nuL_nu against frequency

nu_l_nu = []
for i in range(0,len(nu_midpoints)):
    x = nu_midpoints[i]*integral_values[i]
    nu_l_nu.append(x)
    
log_nu_l_nu = np.log10(nu_l_nu)




plt.figure(dpi=1000)
plt.scatter(log_nu_midpoints, log_nu_l_nu, s=1)
plt.xlabel(r'log[$\nu$ (Hz)]')
plt.ylabel(r'log[$\nu L_{\nu}$ (erg $s^{-1}$)]')
plt.grid(False)
plt.show()


# Calculation for total luminosity

def total_luminosity():
    area =[]
    for i in range(0,len(nu_midpoints)):
        rect = integral_values[i]*(10**nu_grid[i+1]-10**nu_grid[i])
        area.append(rect)
    return sum(area)

print(total_luminosity())