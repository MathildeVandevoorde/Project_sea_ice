# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:23:57 2021

@author: Mathilde
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

#Bottom boundary condition : 
T_bo = -1.8 + Kelvin #Bottom temperature[K]

#Set of variables for the problem 
years_number = 1 #number of year we're working with
days_number = 365*years_number #number of days
h_i = np.zeros(days_number) #ice thickness
h_i[0] = 0.1 #initial conditions on the ice thickness[m]
Q_w = 2 #Ocean heat flux [W/m2]
h_s = 0 #Snow thickness [m]
days = np.arange(days_number)
T_su = np.zeros(days_number) #Surface temperature[K]

    
def getQ_sol(day):
    return 314 *math.exp(-(day-164)**2/4608)
def getQ_nsol(day):
    return 118 * math.exp(-0.5 * (day - 206)**2 / (53**2)) + 179



#Resolution of the problem : 
for i in np.arange(days_number-1):
    doy = i%365 + 1
    T_w = T_bo #mixed layer temperature
    Q_sol = getQ_sol(doy) #solar heat flux [W/m2]
    Q_nsol = getQ_nsol(doy) #non solar heat flux [W/m2]
    resolv = np.roots([-epsilon*sigma,0,0,-ki/h_i[i],ki*T_bo/h_i[i] + Q_sol*(1-alb) + Q_nsol])
    T_su[i] = np.max(resolv) #Temperature[K]
    flux_surf = 0 #surface heat flux [W/m2]
    
    #Melting of the ice if the temperature is above 0 degree celsius :
    if(T_su[i]>=273.15):
        flux_surf = (1-alb)*Q_sol + Q_nsol - epsilon*sigma*T_su[i]**4
        T_su[i] = 273.15 
    
    keff = ki*ks*(h_i[i]+h_s)/(ki*h_s+ks*h_i[i])
    Q_c = keff * (T_su[i] - T_bo)/(h_i[i]+h_s)#sensible heat flux[W/m2]
    h_i[i+1] = h_i[i] - (Q_c + Q_w + flux_surf)*sec_per_day/(rhoi*Lfus)

#We calculate the last temperature, because we only do it until days_number-1 before
resolv = np.max(np.roots([-epsilon*sigma,0,0,-ki/h_i[-1],ki*T_bo/h_i[-1] + getQ_sol(365)*(1-alb) + getQ_nsol(365)]))
T_su[-1] = np.minimum(resolv,273.15)

plt.plot(days,h_i)