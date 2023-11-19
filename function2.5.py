# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy 
import scipy
import matplotlib.pyplot as plt
from scipy.special import logsumexp
import scipy.integrate as integrate
import scipy.special as special

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants in cgs units
g_big_cgs = 6.67430e-8  # dyn·cm^2/g^2
m_sol_cgs = 1.98847e33  # g
c_cgs = 3e10  # cm/s
sigma_cgs = 5.670374419e-5  # erg·cm^(-2)·s^(-1)·K^(-4)
h = 6.62607015e-27  # erg·s
k = 1.380649e-16  # erg/K

# Define a range of nu values on a logarithmic scale in cgs units
nu_min_cgs = 1e14  # Hz
nu_max_cgs = 1e15  # Hz
nu_values = np.logspace(np.log10(nu_min_cgs), np.log10(nu_max_cgs), 500)

# Initialize an empty list to store L_nu * nu values
L_nu_times_nu_values = []

# Define the function to be integrated
def integrand(r, nu):
    temperature = T(r)
    numerator = (2 * np.pi * h * nu**3) / c_cgs**2
    denominator = np.expm1(h * nu / (k * temperature))
    return numerator / denominator

# Define the temperature function in cgs units
def T(rad):
    return (((3 * g_big_cgs * m_sol_cgs * 10 * m_sol_cgs) / (8 * np.pi * rad**3 * (g_big_cgs * m_sol_cgs / c_cgs**2)**3 * sigma_cgs)) * (1 - np.sqrt(6e5 / rad))**0.25)
# Specify the inner and outer radii in cm
r_in = 6e5  # 6 km in cm
r_out = 1e10  # 10,000 km in cm

# Calculate L_nu * nu for each nu value
for nu in nu_values:
    result, _ = quad(integrand, r_in, r_out, args=(nu,))
    L_nu_times_nu_values.append(result * nu)

# Convert L_nu_times_nu_values to a NumPy array
L_nu_times_nu_values = np.array(L_nu_times_nu_values)

# Create the logarithmic plot
plt.figure()
plt.loglog(nu_values, L_nu_times_nu_values)
plt.xlabel('log(nu) (log Hz)')
plt.ylabel('log(L_nu * nu) (log erg/s)')
plt.title('Logarithmic plot of L_nu * nu as a function of log(nu)')
plt.grid(True)
plt.show()
