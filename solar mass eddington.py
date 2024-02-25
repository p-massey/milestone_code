#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:54:06 2024

@author: peterm
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants

g_big_cgs = 6.67430e-8  # dyn·cm^2/g^2
m_sol_cgs = 1.98847e33  # g
c_cgs = 3e10  # cm/s
sigma_cgs = 5.670374419e-5  # erg·cm^(-2)·s^(-1)·K^(-4)
h = 6.62607015e-27  # erg·s
k = 1.380649e-16  # erg/K
m_p = 1.6726219e-24  # Mass of a proton in grams
epsilon_r = 0.1  # Radiative efficiency
sigma_T = 6.6524587158e-25  # Thomson scattering cross-section in cm^2


# input black hole mass in solar masses
def eddington_accretion_rate(bh_mass):
    return (4 * np.pi * g_big_cgs * bh_mass*m_sol_cgs*m_p) / (sigma_T*epsilon_r * c_cgs)

print(eddington_accretion_rate(1e8))


# Define the integration limits

r_in = 6
r_out = 1e5

# Define the range of nu values

nu_min = 1e5
nu_max = 1e20

def nuLnu(m):
    
    r_g = g_big_cgs * (m*m_sol_cgs) / c_cgs**2
    
    def T(rad):
        numerator = 3 * g_big_cgs * eddington_accretion_rate(m) * m * m_sol_cgs
        denominator = 8 * np.pi * (rad)**3 * r_g**3 * sigma_cgs
        bracket = 1 - np.sqrt(r_in/(rad))
        return ((numerator/denominator)*bracket)**0.25
    
    # Create a radial grid with logarithmic spacing

    r_grid = np.linspace(np.log10(r_in), np.log10(r_out), num = 2000)


    # Calculate midpoints of the bins

    r_midpoints = []
    for i in range(0,len (r_grid)-1):
        n = 10**(r_grid[i] + ((r_grid[i+1] - r_grid[i]) / 2.0))
        r_midpoints.append(n)

    # Define the flux function

    def F(nu, r):
            numerator = (2 * np.pi * h * nu**3) / (c_cgs**2)
            denominator = np.exp((h * nu)/(k * T(r)))-1
            return numerator/denominator


    # Create a frequency grid with logarithmic spacing

    nu_grid = np.linspace(np.log10(nu_min), np.log10(nu_max), num = 2000)


    # define the integrand for luminosity function

    def integrand(nu,r):
        return F(nu, r) * 4 * np.pi * r_g**2 * r

    # Integral for luminosity for a given frequency

    def integral(nu):
        area =[]
        for i in range(0,len(r_midpoints)):
            rect = integrand(nu, r_midpoints[i])*(10**r_grid[i+1]-10**r_grid[i])
            area.append(rect)
        return sum(area)
    
    # Find the midpoints for nu

    nu_midpoints = []
    for i in range(0,len (nu_grid)-1):
        n = 10**(nu_grid[i] + ((nu_grid[i+1] - nu_grid[i]) / 2.0))
        nu_midpoints.append(n)
        
    integral_values = [integral(i) for i in nu_midpoints]

    # find set of nuLnu values
    
    nu_l_nu = []
    for i in range(0,len(nu_midpoints)):
        x = nu_midpoints[i]*integral_values[i]
        nu_l_nu.append(x)
    
    return integral_values


nu_grid = np.linspace(np.log10(nu_min), np.log10(nu_max), num = 2000)
nu_midpoints = []
for i in range(0,len (nu_grid)-1):
    n = 10**(nu_grid[i] + ((nu_grid[i+1] - nu_grid[i]) / 2.0))
    nu_midpoints.append(n)

def total_luminosity(mass):
    area =[]
    y=nuLnu(mass)
    for i in range(0,len(nu_midpoints)):
        
        rect = y[i]*(10**nu_grid[i+1]-10**nu_grid[i])
        area.append(rect)
    return sum(area)

xs = [3, 6, 11, 20, 37, 78]
    
mass_grid = [i for i in xs]



nu_grid = np.linspace(np.log10(nu_min), np.log10(nu_max), num = 2000)

nu_midpoints = []

for i in range(0,len (nu_grid)-1):
    n = 10**(nu_grid[i] + ((nu_grid[i+1] - nu_grid[i]) / 2.0))
    nu_midpoints.append(n)


luminosities = [np.log10(total_luminosity(i)) for i in mass_grid]

plt.figure(dpi=1000)
plt.scatter(np.log10(mass_grid), luminosities, s=10)
plt.xlabel(r'$log_{10}$[Mass of black hole ($M_\odot$)]')
plt.ylabel(r'$log_{10}$[Total luminosity (erg)]')
plt.grid(False)
plt.show()

print(total_luminosity(20))