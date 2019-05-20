# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:16:32 2019

@author: Lisa
"""

import numpy as np
import matplotlib.pyplot as plt
from math import * 
from inputs.constants import *
from inputs.performance_inputs import *

def dragcoefficient(Cfe, Swet_S):
    CD0=Cfe*Swet_S
    CD0_TO=CD0+d_CD0_to
    CD0_land=CD0+d_CD0_landing
    
    return CD0, CD0_TO, CD0_land


def stallspeedlanding(fieldlenghtlanding):
        V_a=(fieldlenghtlanding/0.3)**0.5*0.5144444
        V_s=V_a/1.3
        return V_s
    
def wingloading_landing(CL,V):
    WS=0.5*rho_0*V**2*CL
    return WS

def wingloading_takeoff(fieldlenghtlanding, CL,f):
    WS_TO=(CL*rho_0*fieldlenghtlanding*0.3048/0.5915)/(2*f)
    
    return WS_TO

def thurstloading_takeoff(TOP,sigma,CL,WSrange):
    TW_TO= (WSrange) / (sigma*TOP*CL)    
    return TW_TO

def thrustloading_cruise(CD0,WSrange):
    TW_cruise=(rho_0/rho)**0.75*(CD0*0.5*rho*V_cruise**2/(WSrange)+(WSrange)*1/(pi*A*e*0.5*rho*V_cruise**2))
    return TW_cruise

def thrustloading_climbrate(c,CL,CD,WSrange):
    TW_climb=c/(WSrange*2/(rho_0*CL))**0.5 + CD/CL
    return TW_climb
def thrustloading_climbgradient(c,V,CD0):
    TW_cV=c/(V)+2*(CD0/(pi*A*e))**0.5
    return TW_cV

def get_WSrange(WSstart,WSend,D_WS):
    return np.linspace(WSstart,WSend,D_WS)


    
def plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,c,f,sigma, TOP, CD0,WSstart,WSend,D_Ws):
    
    CD_climb=4*CD0
    Cl_climb=(3*CD0*pi*A*e)**0.5
    
    WSrange=get_WSrange(WSstart,WSend,D_Ws)
    
    V_s=stallspeedlanding(Sland)

    V_climb=300/3.6
    WS_L=wingloading_landing(Cl_land, V_s)
    WS_TO=wingloading_takeoff(Sland,Cl_land,f)
    
    TW_TO=thurstloading_takeoff(TOP,sigma,Cl_TO,WSrange)
    TW_cruise=thrustloading_cruise(CD0,WSrange)
    
    TW_climb=thrustloading_climbrate(c,Cl_climb,CD_climb,WSrange)
    TW_cV= thrustloading_climbgradient(c,V_climb,CD0)
    
    fig = plt.figure(figsize = (12,5))
    ax = fig.add_subplot(111)
    ax.plot(WSrange,TW_cruise, label='Cruise condition ')      
    ax.plot(WSrange,TW_climb,label='Climb rate condition' ) 
    ax.axhline(TW_cV,label='Climb gradient condition' ) 
    ax.axvline(WS_TO, label='Take off W/S',color='r')
    ax.axvline(WS_L, label='Landing W/S')
    ax.plot(WSrange,TW_TO,label='Take-off condition ')
    ax.fill_between(WSrange, TW_climb, 2, where = (TW_cruise <= TW_climb) * (WSrange < WS_L), facecolor='green', interpolate=False)
    ax.fill_between(WSrange, TW_cruise, 2, where = (TW_cruise >= TW_climb) * (WSrange < WS_L), facecolor='green', interpolate=False)
    ax.scatter(4359,0.33, color = 'r',linewidths=5, label = 'Design point')
    ax.set(title = 'Wing/Thrust Loading Diagram', xlabel = 'W/S[-]', ylabel = 'T/W[-]', ylim = (0,max(TW_climb)), xlim = (0,5500))
    ax.legend(loc='upper right')
    plt.savefig('./Output/loading_diagram.pdf')
    return TW_cruise,TW_climb,WSrange,WS_L


CD0, CD0_TO, CD0_land=dragcoefficient(Cfe,Swet_S)
#
#loadingdiagram=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,c,f,sigma, TOP, CD0,100,7100,100)
#WS=np.arange(100,7100,100)
#
#WS_landing=0.5*rho_0*V_s**2*CL_land_max
#
#WS_TO=(CL_land_av*rho_0*Sland*0.3048/0.5915)/(2*f)
##take off parameter
#TW=WS*1/CL_TO_av*1/sigma*1/TOP
#
##cruise performance and climb performance 
#
#CD=CD0+CL_clean_av**2/(pi*e*A)
#
#CD_climb=4*CD0
#CL_climb=(3*CD0*pi*A*e)**0.5
#
#TW_cruise=(rho_0/rho)**0.75*(CD0*0.5*rho*V**2/(WS)+(WS)*1/(pi*A*e*0.5*rho*V**2))
#TW_climb=c/(WS*2/rho_0*1/CL_climb)**0.5 + CD_climb/CL_climb
#TW_cV=c/(300/3.6)+2*(CD0/(pi*A*e))**0.5
#
#
#plt.figure(1)
#plt.plot(WS,TW_TO, label='Cruise condition ')      
#plt.plot(WS,TW_climb,label='Climb rate condition' ) 
#plt.axhline(TW_cV,label='Climb gradient condition' ) 
#plt.axvline(WS_TO, label='Take off W/S',color='r')
#plt.axvline(WS_landing, label='Landing W/S')
#plt.plot(WS,TW,label='Take-off condition ')
#plt.title('Wing/Thrust Loading diagram')
#plt.xlabel('W/S[-]',size=16)
#plt.ylabel('T/W[-]',size=16)
#plt.legend()
#plt.show()
#
#
