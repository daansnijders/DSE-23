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
TOP=175*4.8824*9.81         #N/M^2
sigma=1

#statistical input
f=0.84
Swet_S=6
Cfe=0.0045

d_CD0_to=0.015
d_CD0_landing=0.065
d_CD0_gear=0.02



e=0.85
c=15
A=11

def dragcoefficient(Cfe, Swet_S):
    CD0=Cfe*Swet_S
    CD0_TO=CD0+d_CD0_to
    CD0_land=CD0+d_CD0_landing
    
    return CD0, CD0_TO, CD_land


class wingthrustloadingdiagram:
    def stallspeedlanding(Fieldlenghtlanding):
        V_a=(fieldlenghtlanding/0.3)**0.5*0.5144444
        V_s=V_a/1.3
        return V_s
    
    def wingloading(rho,CL,V):
        WS=0.5*rho*V**2*CL
        return WS
    
    def wingloading_takeoff(rho,fieldlenghtlanding, CL,f):
        WS_TO=(CL*rho*fieldlenghtlanding*0.3048/0.5915)/(2*f)
        
        return WS_TO
    def thurstloading_takeoff(TOP,sigma,CL,WSrange):
        TW_TO=WSrange/(sigma*TOP*CL)    
        return TW_TO
    
    def thrustloading_cruise(rho_0,rho,CD0,V,A,e,WSrange):
        TW_cruise=(rho_0/rho)**0.75*(CD0*0.5*rho*V**2/(WS)+(WS)*1/(pi*A*e*0.5*rho*V**2))
        return TW_cruise
    
    def thrustloading_climbrate(c,rho,CL,CD,WSrange):
        TW_climb=c/(WSrange*2/rho*1/CL)**0.5 + CD/CL
        return TW_climb
    def thrustloading_climbgradient(c,V,CD0,A,e):
        TW_cV=c/(V)+2*(CD0/(pi*A*e))**0.5
        return TW_cV
    
    def WSrange(WSstart,WSend,D_WS):
        WSrange=np.arange(WSstart,WSend,D_WS)
    
    
    
    def plotwingloadingdiagram(Sland,rho,rho_0,Cl_TO,Cl_clean,Cl_land,c,Vcruise,Vclimb,f,sigma, TOP, A, e, CD0):
        CD_climb=4*CD0
        CL_climb=(3*CD0*pi*A*e)**0.5
        
        WSrange=WSrange(WSstart,Wsend,D_Ws)
        
        V_s=stallspeedlanding(Sland)
        WS_L=wingloading(rho_0,CL_land, V_s)
        WS_TO=wingloading(rho_0,Sland,Cl_TO,f)
        TW_TO=thurstloading_takeoff(TOP,sigma,Cl_TO,WSrange)
        TW_cruise=thrustloading_cruise(rho_0,rho,CD0,Vcruise,A,e,WSrange)
        TW_climb=thrustloading_climbrate(c,rho_0,CL_climb,CD_climb,WSrange)
        TW_cV= thrustloading_climbgradient(c,Vclimb,CD0,A,e)
        
        plt.figure(1)
        plt.plot(WSrange,TW_TO, label='Cruise condition ')      
        plt.plot(WSrange,TW_climb,label='Climb rate condition' ) 
        plt.axhline(TW_cV,label='Climb gradient condition' ) 
        plt.axvline(WS_TO, label='Take off W/S',color='r')
        plt.axvline(WS_L, label='Landing W/S')
        plt.plot(WSrange,TW,label='Take-off condition ')
        plt.title('Wing/Thrust Loading diagram')
        plt.xlabel('W/S[-]',size=16)
        plt.ylabel('T/W[-]',size=16)
        plt.legend()
        plt.show()
                           
   plt.plot(rangeWS,)
    plt.plot




WS=np.arange(100,7100,100)

WS_landing=0.5*rho_0*V_s**2*CL_land_max

WS_TO=(CL_land_av*rho_0*Sland*0.3048/0.5915)/(2*f)
#take off parameter
TW=WS*1/CL_TO_av*1/sigma*1/TOP

#cruise performance and climb performance 

CD=CD0+CL_clean_av**2/(pi*e*A)

CD_climb=4*CD0
CL_climb=(3*CD0*pi*A*e)**0.5

TW_cruise=(rho_0/rho)**0.75*(CD0*0.5*rho*V**2/(WS)+(WS)*1/(pi*A*e*0.5*rho*V**2))
TW_climb=c/(WS*2/rho_0*1/CL_climb)**0.5 + CD_climb/CL_climb
TW_cV=c/(300/3.6)+2*(CD0/(pi*A*e))**0.5


plt.figure(1)
plt.plot(WS,TW_TO, label='Cruise condition ')      
plt.plot(WS,TW_climb,label='Climb rate condition' ) 
plt.axhline(TW_cV,label='Climb gradient condition' ) 
plt.axvline(WS_TO, label='Take off W/S',color='r')
plt.axvline(WS_landing, label='Landing W/S')
plt.plot(WS,TW,label='Take-off condition ')
plt.title('Wing/Thrust Loading diagram')
plt.xlabel('W/S[-]',size=16)
plt.ylabel('T/W[-]',size=16)
plt.legend()
plt.show()


