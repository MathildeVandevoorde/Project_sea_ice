# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 10:16:28 2021

@author: Mathilde

This code calculates the surface temperature.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

#Physical constants : 
Lfus= 3.35e5 #Latent  heat  of  fusion  for  water  [J/kg] 
rhoi = 917 #Sea ice density [kg/m3]
ki = 2.2 #Sea ice thermal conductivity [W/m/K]
ks = 0.31 #Snow thermal conductivity [W/m/K]
sec_per_day = 86400 #Seconds in one day [s/day]
epsilon = 0.99 #surface emissivity 
sigma = 5.67e-8 #Stefan-Boltzmann constant
Kelvin = 273.15 #Conversion from Celsius to Kelvin
alb = 0.8 #Surface albedo

#Boundary condition : 
T_bo = -1.8 + Kelvin #Bottom temperature[K]
h_i = 0.5 #Sea ice thickness [m]
days_number = 365 #number of days we're working with
days = np.arange(days_number) #Number of days we're working with
T_su = np.zeros(days_number) #Surface temperature[K]

    
def getQ_sol(day):
    return 314 *math.exp(-(day-164)**2/4608)
def getQ_nsol(day):
    return 118 * math.exp(-0.5 * (day - 206)**2 / (53**2)) + 179


for i in days:
    resolv = np.roots([-epsilon*sigma, 0, 0,-ki/h_i, ki*T_bo/h_i + getQ_sol(i+1)*(1-alb) + getQ_nsol(i+1)])
    sol = np.max(resolv)
    T_su[i] = np.minimum(sol,273.15)
 

plt.plot(days,T_su)