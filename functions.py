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

g_big = 6.67430E-11 #N m^2 kg^-2
m_sol = 1.98847E30 #kg
c = 3E8 #m s^-1
sigma = 5.670374419E-8 #W m^-2 K^-4
h = 6.62607015E-34 
k = 1.380649E-23 

m     = 10*m_sol
m_dot = 10.0E15  #kg s^-1
r_g   = g_big*m/c**2
r_in = 6 
r_out = 10.0E5
nu_min = 10.0E14
nu_max = 10.0E19



    