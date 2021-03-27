# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 09:51:13 2021

@author: Mathilde
"""


import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import optimize

#Physical constants : 
Lfus= 3.35e5 #Latent  heat  of  fusion  for  water  [J/kg] 
rhoi = 917 #Sea ice density [kg/m3]
ki = 2.2 #Sea ice thermal conductivity [W/m/K]
ks = 0.31 #Snow thermal conductivity [W/m/K]
sec_per_day = 86400 #Seconds in one day [s/day]
epsilon = 0.99 #surface emissivity 
sigma = 5.67e-8 #Stefan-Boltzmann constant
Kelvin = 273.15 #Conversion from Celsius to Kelvin
alb = 0.6 #Surface albedo
alb_w = 0.1 #Albedo of water
h_w = 50 #depth of the mixed layer [m]
rhow = 1025 #water density [kg/m3]
c_w = 4000 #Heat capacity of water [J/(kg K)]


#Bottom boundary condition : 
T_bo = -1.8 + Kelvin #Bottom temperature[K]

#Set of variables for the problem 
years_number = 50 #number of years we're working with
days_number = 365*years_number #number of days
h_i = np.zeros(days_number) #ice thickness
h_i[0] = 0.1 #initial conditions on the ice thickness[m]
Q_w = 5 #Ocean heat flux [W/m2]
h_s = np.zeros(days_number) #Snow thickness [m]
days = np.arange(days_number)
T_su = np.zeros(days_number) #Surface temperature[K]
T_w = np.zeros(days_number)
T_w[0] = T_bo #mixed layer temperature


    
def getQ_sol(day):
    return 314 *math.exp(-(day-164)**2/4608)
def getQ_nsol(day):
    return 118 * math.exp(-0.5 * (day - 206)**2 / (53**2)) + 179

x0 = 200.15
def f(T_su,h,day):
    return (-(h*epsilon*sigma)/ki*T_su**4 - T_su + h/ki*((getQ_sol(day))*(1-alb)+getQ_nsol(day))+T_bo)
 
def df(T_su,h,day):
    return (-(4*h*epsilon*sigma)/ki*T_su**3-1)


def get_root(h,day):
    root = optimize.newton(f,x0,fprime = df, args= (h,day))
    return(root)

#Resolution of the problem : 
for i in np.arange(days_number-1):
    if (h_i[i]>0): #if there is ice
        doy = i%365 + 1 #day of the year
        Q_sol = getQ_sol(doy) #solar heat flux [W/m2]
        Q_nsol = getQ_nsol(doy) #non-solar heat flux [W/m2]
        T_su[i] = get_root(h_i[i],doy) #surface temperature[K]
        keff = ki*ks*(h_i[i]+h_s[i])/(ki*h_s[i]+ks*h_i[i])
        flux_surf = 0 #surface heat flux [W/m2]
        T_w[i+1] = T_bo
        
        #Melting of the ice if the temperature is above 0 degree celsius :
        if(T_su[i]>=273.15):
            T_su[i] = 273.15
            flux_surf = (1-alb)*Q_sol + Q_nsol - epsilon*sigma*T_su[i]**4
            
        Q_c = keff * (T_su[i] - T_w[i])/(h_i[i]+h_s[i]) #sensible heat flux[W/m2]
        dhdt = - (Q_c + Q_w + flux_surf)*sec_per_day/(rhoi*Lfus)
        h_i[i+1] = h_i[i] + dhdt
        
    if (h_i[i]<=0): #if there is no ice
        doy = i%365 + 1
        Q_sol = getQ_sol(doy)
        Q_nsol = getQ_nsol(doy)
        flux_surf = Q_sol*(1 - alb_w) + Q_nsol - epsilon*sigma*Kelvin**4
        T_w[i+1] = T_w[i] + (flux_surf)*sec_per_day/(c_w*rhow*h_w) #Water temperature[K], it changes if there is no ice
        h_i[i] = 0 #we can't have negative ice
        T_su[i] = Kelvin
        #Regrowth of ice if the temperature of the ocean is below the freezing temperature :
        if (T_w[i+1] <= T_bo):
            dhdt = -(flux_surf)*sec_per_day/(Lfus*rhoi)
            h_i[i+1] = h_i[i] + dhdt
            T_w[i+1] = T_bo
        else :
            h_i[i+1] = 0

        
#We calculate the last temperature, because we only do it until days_number-1 before
resolv = np.max(np.roots([-epsilon*sigma,0,0,-ki/h_i[-1],ki*T_bo/h_i[-1] + getQ_sol(365)*(1-alb) + getQ_nsol(365)]))
T_su[-1] = np.minimum(resolv,273.15)

#Plots :
fig, axs = plt.subplots(2,2)
axs[0,0].plot(days,h_i)
axs[0,0].set_title("Ice thickness")
axs[0,0].set_xlabel('time(day)')
axs[0,0].set_ylabel('thickness(m)')
axs[0,1].plot(days,h_s)
axs[0,1].set_title("Snow thickness")
axs[0,1].set_xlabel('time(day)')
axs[0,1].set_ylabel('thickness(m)')
axs[1,0].plot(days,T_w)
axs[1,0].set_title('Water temperature')
axs[1,0].set_xlabel('time (day)')
axs[1,0].set_ylabel('temperature(K)')
axs[1,1].plot(days,T_su)
axs[1,1].set_title('Surface temperature')
axs[1,1].set_xlabel('time (day)')
axs[1,1].set_ylabel('temperature(K)')
fig.tight_layout()
