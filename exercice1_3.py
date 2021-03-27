# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 09:16:23 2021

@author: Mathilde

This code adds a constant amount of snow on top of the ice
"""
import numpy as np
import matplotlib.pyplot as plt

#Physical constants : 
Lfus= 3.35e5 #Latent  heat  of  fusion  for  water  [J/kg] 
rhoi = 917 #Sea ice density [kg/m3]
ki = 2.2 #Sea ice thermal conductivity [W/m/K]
ks = 0.31 #Snow thermal conductivity [W/m/K]
sec_per_day = 86400 #Seconds in one day [s/day]

#Bottom boundary condition : 
T_bo = -1.8 #Bottom temperature[C]

#Set of variables for the problem : 
T_air = -10 #Temperature of the air
days_number = 30 #number of days
h_i = np.zeros(30) #ice thickness
h_i[0] = 0.1 #initial conditions on the ice thickness[m]
Q_w = 5 #Ocean heat flux [W/m2]
h_s = 0.05 #Snow thickness [m]
days = np.arange(30)

#Resolution of the problem : 
for i in np.arange(days_number-1):
    keff = ki*ks*(h_i[i]+h_s)/(ki*h_s+ks*h_i[i])
    Q_c = keff * (T_air - T_bo)/(h_i[i]+h_s) #sensible heat flux
    dhdt = - (Q_c+Q_w)*sec_per_day /(rhoi*Lfus)
    h_i[i+1] = h_i[i] + dhdt

fig, axs = plt.subplots(2,2)
axs[0,0].plot(days,h_i)
axs[0,0].set_title("Ice thickness")
axs[0,0].set_xlabel('time(day)')
axs[0,0].set_ylabel('thickness(m)')
fig.tight_layout()