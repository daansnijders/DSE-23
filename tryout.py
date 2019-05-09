# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:16:32 2019

@author: Lisa
"""

import numpy as np
import matplotlib.pyplot as plt
from math import * 
#Input in CL from Roskam 
CL_clean_min=1.2
CL_clean_max=1.8
CL_clean_av=1.5

CL_TO_min=1.6
CL_TO_max=2.2
CL_TO_av=1.9

CL_land_min=1.8
CL_land_max=2.8
CL_land_av=2.3

#flight conditions
rho_0=1.225
rho=0.3482320451
V=0.75*(1.4*216.65*287.05)**0.5

#take off and landing
Sland=1500/0.3048           #field length in feet
Stakeoff=2000/0.3048        # field length in feet
TOP=175*4.8824*9.81    #N/M^2
sigma=1

#statistical input
f=0.84
Swet_S=6
Cfe=0.0045



d_CD0_to=0.015
d_CD0_landing=0.065
d_CD0_gear=0.02

CD0=Cfe*Swet_S
CD0_TO=CD0+d_CD0_to
CD0_land=CD0+d_CD0_landing

e=0.85
c=15


Alist=[8,10,12,14,16]


#calculate the  stall speed at landing
V_a=(Sland/0.3)**0.5*0.5144444
V_s=V_a/1.3

#using the stall speed to calculate the max W/S
WS_landing=0.5*rho_0*V_s**2*CL_land_max

WS_TO=(CL_land_av*rho_0*Sland*0.3048/0.5915)/(2*f)



WS=np.arange(100,7100,100)


#take off parameter
TW=WS*1/CL_TO_av*1/sigma*1/TOP

#cruise performance and climb performance 
for i in range(len(Alist)):
    print (i)
    A=Alist[i]
    CD=CD0+CL_clean_av**2/(pi*e*A)
    
    CD_climb=4*CD0
    CL_climb=(3*CD0*pi*A*e)**0.5
    
    TW_TO=(rho_0/rho)**0.75*(CD0*0.5*rho*V**2/(WS)+(WS)*1/(pi*A*e*0.5*rho*V**2))
    TW_climb=c/(WS*2/rho_0*1/CL_climb)**0.5 + CD_climb/CL_climb
    TW_cV=c/V+2*(CD0/(pi*A*e))**0.5
    
    
    plt.figure(1)
    plt.plot(WS,TW_TO,label='cruise condition at A= %2.3f' %A )      
    plt.plot(WS,TW_climb,label='climb condition at A= %2.3f' %A ) 
    plt.axhline(TW_cV,label='CV condition at A= %2.3f' %A ) 
 
    
plt.figure(1)
plt.axvline(WS_TO, label='take off W/S')
plt.axvline(WS_landing, label='landing W/s')
plt.plot(WS,TW,label='min CL takeoff TW in function of WS')
plt.legend()
plt.show()



#weights of the concepts and configurations

W_baseline=506305.22
W_maxconf=687128.52