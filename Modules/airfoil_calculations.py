# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:01:38 2019

@author: Anique
"""
from math import *
from inputs.constants import *
from inputs.performance_inputs import *

#### Inputs needed for the definition, they are from requirements or from previous sizing ####
#### This is an example for concept 1, configuration 1 ####



Cr = 4.91489505      # Root chord [m]
Ct = 1.520444671     # Tip chord [m]
#MTOW = 484363.0809   # Maximum Take-off Weight [N]




#### Definition calculates different Reynold numbers for different locations 
#### and calculates Cl design for the airfoils
def airfoil( Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, sweep_c2, b, Taper):
    avgC = [(Cr[i] + Ct[i])/2   for i in range(3)]   # Average chord [m]
    Re1 =[ (rho*V_cruise*Cr[i])/mu   for i in range(3)]   # Reynolds number at root chord [-]
    Re2 = [(rho*V_cruise*avgC[i])/mu for i in range(3)]   # Reynolds number at avg chord [-]
    Re3 = [(rho*V_cruise*Ct[i])/mu   for i in range(3)]   # Reynolds number at tip chord [-]
    q = 0.5*rho*V_cruise**2        # Dynamic pressure [Pa]
    WSbegin = [(MTOW[i]*FF1*FF2*FF3*FF4)/S[i]   for i in range(3)]        # Wing Loading begin cruise  
    WSend = [(MTOW[i]*FF1*FF2*FF3*FF4*FF5)/S[i]  for i in range(3)]        # Wing loading end cruise
    CLdes =[ 1.1*(1/q)*(0.5*(WSbegin[i] + WSend[i]))   for i in range(3)] # Wing CL design [-]
                                             
    Veff =[ V_cruise*cos(sweep_c2[i])  for i in range(3)]              # Effective velocity [m/s]
    qeff =[ 0.5*rho*Veff[i]**2   for i in range(3)]                    # Effective Dynamic Pressure [Pa]
    Cl_des = [(q*CLdes[i])/qeff  [i] for i in range(3)]                   # Airfoil Cl design [-]
    
    return(Re1, Re2, Re3, Cl_des)
    

