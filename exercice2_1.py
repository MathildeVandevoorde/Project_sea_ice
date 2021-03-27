# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 09:57:27 2021

@author: Mathilde Vandevoorde

This code gives us the solar and non solar heat fluxes depending on the days of the year
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

def getQ_sol(day):
    return 314 *math.exp(-(day-164)**2/4608)
def getQ_nsol(day):
    return 118 * math.exp(-0.5 * (day - 206)**2 / (53**2)) + 179

day = np.arange(365)
Q_sol = np.zeros(365) #solar heat flux [W/m2]
Q_nsol = np.zeros(365) #non-solar heat flux [W/m2]

for i in day:
    Q_sol[i] = getQ_sol(i)
    Q_nsol[i] = getQ_nsol(i)

plt.plot(day, Q_sol)
plt.plot(day,Q_nsol)