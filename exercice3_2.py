# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:29:13 2021

@author: Mathilde Vandevoorde

Code final, avec formation de snow-ice ou non
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
alb = 0.8 #Surface albedo
alb_w = 0.1 #Albedo of water
h_w = 50 #depth of the mixed layer [m]
rhow = 1025 #water density [kg/m3]
c_w = 4000 #Heat capacity of water [J/(kg K)]
rhos = 330#density of snow[kg/m3]

#Bottom boundary condition : 
T_bo = -1.8 + Kelvin #Bottom temperature[K]

#Set of variables for the problem 
years_number = 50 #number of years we're working with
days_number = 365*years_number #number of days
h_i = np.zeros(days_number) #ice thickness
h_i[0] = 0.1 #initial conditions on the ice thickness[m]
Q_w = 5 #Ocean heat flux [W/m2]
h_s = np.zeros(days_number)#Snow thickness [m]
days = np.arange(days_number)
T_su = np.zeros(days_number) #Surface temperature[K]
T_w = np.zeros(days_number) #Water Temperature [K]
T_w[0] = T_bo #mixed layer temperature
snow_ice = 0 #0 if you don't want snow-ice formation, 1 if you want snow ice formation

    
def getQ_sol(day):
    return 314 *math.exp(-(day-164)**2/4608)
def getQ_nsol(day):
    return 118 * math.exp(-0.5 * (day - 206)**2 / (53**2)) + 179

def snowfall(day):
    if (doy<120):
        return 0.05/181
    if (doy>304):
        return 0.05/181
    if (doy>232):
        return 0.3/71
    if(doy<151):
        return 0.05/31
    else:
        return 0

#Function to resolve the quadratic equation : 
x0 = 205.15
def f(T_su,h,day):
    return (-(h*epsilon*sigma)/ki*T_su**4 - T_su + h/ki*((getQ_sol(day))*(1-alb)+getQ_nsol(day))+T_bo)
 
def df(T_su,h,day):
    return (-(4*h*epsilon*sigma)/ki*T_su**3-1)

def get_root(h,day):
    root = optimize.newton(f,x0,fprime = df, args= (h,day))
    return(root)

#Resolution of the problem : 
for i in np.arange(days_number-1):
    doy = i%365 + 1
    if (h_i[i]>0): #if there is ice
        doy = i%365 + 1 #day of the year
        Q_sol = getQ_sol(doy) #solar heat flux [W/m2]
        Q_nsol = getQ_nsol(doy) #non-solar heat flux [W/m2]
        T_su[i] = get_root(h_i[i],doy) #surface temperature[K]
        flux_surf = 0 #surface heat flux [W/m2]
        T_w[i+1] = T_bo
        
        #Melting of the ice/snow if the temperature is above 0 degree celsius :
        if(T_su[i]>=273.15):
            T_su[i] = 273.15
            flux_surf = (1-alb)*Q_sol + Q_nsol - epsilon*sigma*T_su[i]**4
        
        keff = ki*ks*(h_i[i]+h_s[i])/(ki*h_s[i]+ks*h_i[i])
        
        #If there is snow and it can melt : 
        if h_s[i] > 0 and T_su[i]==273.15:
            h_s[i+1] = h_s[i]+snowfall(doy) #We add the snowfall to the snow
            h_s[i+1] = h_s[i+1] - flux_surf/(rhos*Lfus)*sec_per_day #We melt the snow with the surface heat flux
            Q_c = keff * (T_su[i] - T_w[i+1])/(h_i[i]+h_s[i]) #sensible heat flux[W/m2]
            h_i[i+1] = h_i[i] - (Q_c + Q_w)*sec_per_day/(rhoi*Lfus)#Ice still have the 2 other fluxes
            if h_s[i+1] <= 0:
                h_s[i+1] = 0  #Snow can't go beneath 0
        #If the snow can't melt/ if there is no snow:
        else :
            #To avoid a bug where the ice can't refreeze :
            if h_i[i] <= 0.1: 
                Q_c = keff * (T_su[i] - T_w[i+1])/(h_i[i]+h_s[i]) #sensible heat flux[W/m2]
                h_i[i+1] = h_i[i] - (Q_c + Q_w + flux_surf)*sec_per_day/(rhoi*Lfus)
            else:
                Q_c = keff*(T_su[i] - T_w[i+1])/(h_i[i]+h_s[i]) #sensible heat flux[W/m2]
                h_s[i+1] = h_s[i] +snowfall(doy)
                h_i[i+1] = h_i[i] - (Q_c + Q_w+ flux_surf)/(rhoi*Lfus)*sec_per_day  

    if (h_i[i]<=0): #if there is no ice
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
            
    #Snow-ice formation
    h_d = (rhos*h_s[i] + rhoi*h_i[i])/rhow #snow-ice draft
    #We check if there is snow beneath the sea level :
    if h_d>h_i[i] and snow_ice == 1: 
        deltah = (rhos*h_s[i] + rhoi*h_i[i] - rhow*h_i[i])/(rhos+rhow-rhoi) #Snow that is beneath the sea level
        h_s[i+1] = h_s[i] - deltah #We remove the snow beneath the sea level
        h_i[i+1] = h_i[i] + deltah #We add the snow beneath the sea level to the ice

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


