#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:00:42 2023

@author: peterm
"""

import numpy
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Constants
g_big_cgs = 6.67430e-8  # dyn·cm^2/g^2
m_sol_cgs = 1.98847e33  # g
c_cgs = 3e10  # cm/s
sigma_cgs = 5.670374419e-5  # erg·cm^(-2)·s^(-1)·K^(-4)
h = 6.62607015e-27  # erg·s
k = 1.380649e-16  # erg/K
m_dot = 

# Define the Schwarzschild radius
m = 1e5  # Replace with the desired mass
r_g = g_big_cgs * m / c_cgs**2

# Define the integration limits
r_in = 6
r_out = 1e6

# Define the range of nu values
nu_min = 1e14  # Replace with the desired nu_min
nu_max = 1e19  # Replace with the desired nu_max

# Create an array to store the results
nu_values = numpy.logspace(numpy.log10(nu_min), numpy.log10(nu_max), num=100)  # Adjust the 'num' parameter as needed
nu_L_nu_values = []

# Loop over the range of nu values and calculate nu * L_nu for each nu
for nu_value in nu_values:
    def integrand(r):
        temperature = (((3 * g_big_cgs * m * m_dot) / (8 * numpy.pi * r**3 * r_g**3 * sigma_cgs)) * (1 - numpy.sqrt(r_in / r)))**0.25
        numerator = (2 * numpy.pi * h * nu_value**3) / c_cgs**2
        denominator = numpy.expm1(h * nu_value / (k * temperature))-1
        return (numerator / denominator) * 4 * numpy.pi * r_g**2 * r

    result, error = quad(integrand, r_in, r_out)
    nu_L_nu = nu_value * result
    nu_L_nu_values.append(nu_L_nu)

# Create the plot
plt.figure(figsize=(8, 6))
plt.loglog(nu_values, nu_L_nu_values)
plt.xlabel(r'$\nu$ (Hz)')
plt.ylabel(r'$\nu \cdot L_{\nu}$ (erg/s)')
plt.title(r'$\nu \cdot L_{\nu}$ vs. $\nu$')
plt.grid(True)
plt.show()



