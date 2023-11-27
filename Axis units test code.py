#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 17:10:48 2023

@author: peterm
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 15:46:27 2023

@author: peterm
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 18:20:18 2023

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

# Accretion rate

m_dot = 1e18

# Define the Schwarzschild radius

m = 10*m_sol_cgs  # Replace with the desired mass
r_g = g_big_cgs * m / c_cgs**2

# Define the integration limits

r_in = 6
r_out = 1e5

# Define the range of nu values

nu_min = 1e10  
nu_max = 1e30 

x_grid = np.linspace(1, 100, 10)
z_grid = np.linspace(2, 20, 10)

def y(i,j):
    return i**j


plt.figure(dpi=600)
for a in z_grid:
    plt.scatter(x_grid, y(x_grid, a), s=1, label= f'$\dot M$ = {a} ' r'$g^{-1}$')
plt.xlabel(r'$log_{10}$[$\nu$ (Hz)]')
plt.ylabel(r'$log_{10}$[$\nu L_{\nu}$ (erg $s^{-1}$)]')
plt.grid(False)
plt.legend()
plt.show()



