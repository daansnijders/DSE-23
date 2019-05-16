# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:01:38 2019

@author: Anique
"""
from math import *
from inputs.constants import *
from inputs.performance_inputs import *

#### airfoil calculates different Reynold numbers for different locations 
#### and calculates Cl design for the airfoils and CL max for the wing
def airfoil( Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, sweep_le, sweep_c2, b, Taper, A, Cl_max):
    avgC = [(Cr[i] + Ct[i])/2   for i in range(3)]   # Average chord [m]
    Re1 = [(rho*V_cruise*Cr[i])/mu   for i in range(3)]   # Reynolds number at root chord [-]
    Re2 = [(rho*V_cruise*avgC[i])/mu for i in range(3)]   # Reynolds number at avg chord [-]
    Re3 = [(rho*V_cruise*Ct[i])/mu   for i in range(3)]   # Reynolds number at tip chord [-]
    q = 0.5*rho*V_cruise**2        # Dynamic pressure [Pa]
    
    WSbegin = [(MTOW[i]*g*FF1*FF2*FF3*FF4)/S[i]   for i in range(3)]        # Wing Loading begin cruise  
    WSend = [(MTOW[i]*g*FF1*FF2*FF3*FF4*FF5)/S[i]  for i in range(3)]        # Wing loading end cruise
    CLdes =[ 1.1*(1/q)*(0.5*(WSbegin[i] + WSend[i]))   for i in range(3)] # Wing CL design [-]
                                             
    Veff =[ V_cruise*cos(sweep_le[i])  for i in range(3)]              # Effective velocity [m/s]
    qeff =[ 0.5*rho*Veff[i]**2   for i in range(3)]                    # Effective Dynamic Pressure [Pa]
    Cl_des = [(q*CLdes[i])/qeff  [i] for i in range(3)]                   # Airfoil Cl design [-]
    
    beta = sqrt(1-M_cruise**2)      # Prandtl-Glauert compressibility correction factor
    CL_alpha = [(2*pi*A[i])/(2 + sqrt(4+(A[i]*beta/0.95)**2*(1+(tan(sweep_c2[i])**2)/beta**2)))  for i in range(3)]
    
    #Wing CLmax for different Re numbers
    CLmax = [0.8*Cl_max[i]-0.24 for i in range(3)]
    
    return(Re1, Re2, Re3, CLdes, Cl_des, CL_alpha, CLmax)
       

#### drag1 calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 1   
def drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = []
    
    for i in range(3):
        if i == 0:
            Wing.append(1.07*2*S[0]*0.003)
        if i == 1:
            Wing.append(1.07*2*S[2]*0.003)
        if i == 2:
            Wing.append(1.07*2*S[2]*0.003 + (((45*inch_to_m)/2)**2*pi*2 + (45*inch_to_m*pi*308.4*inch_to_m))*0.006)
            
    l2 = [l_fuselage[i] - l_nose[i] - l_tailcone[i] for i in range(3)]
    Fuselage = [0.0024*(pi*D[i]/4)*(1/(3*l_nose[i]**2)*((4*l_nose[i]**2+D[i]**2/4)**(1.5)-D[i]**3/8)-D[i]+4*l2[i]+2*sqrt(l_tailcone[i]**2+D[i]**2/4)) for i in range(3)]
    Nacelle = (2*(Dnacel/2)**2*pi + Dnacel*pi*Lnacel)*0.0060
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/S[0]*(Wing[i] + Fuselage[i] + Nacelle + Tailplane[i])*1.1  for i in range(3)]
    
    e = [4.61*(1-0.045*A[i]**0.68)*(cos(sweep_le[i]))**0.15 -3.1 for i in range(3)]
    
    CDcruise = [CD0[i] + 1/(pi*A[i]*e[i])*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD)
    

#### drag calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 2  
def drag2(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = [1.07*2*S[i]*0.003 for i in range(3)]  
    l2 = [l_fuselage[i] - l_nose[i] - l_tailcone[i] for i in range(3)]
    Fuselage = [0.0024*(pi*D[i]/4)*(1/(3*l_nose[i]**2)*((4*l_nose[i]**2+D[i]**2/4)**(1.5)-D[i]**3/8)-D[i]+4*l2[i]+2*sqrt(l_tailcone[i]**2+D[i]**2/4)) for i in range(3)]
    Nacelle = (2*(Dnacel/2)**2*pi + Dnacel*pi*Lnacel)*0.0060
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/S[i]*(Wing[i] + Fuselage[i] + Nacelle + Tailplane[i])*1.1  for i in range(3)]
    
    e = [4.61*(1-0.045*A[i]**0.68)*(cos(sweep_le[i]))**0.15 -3.1 for i in range(3)]
    
    CDcruise = [CD0[i] + 1/(pi*A[i]*e[i])*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD)
    
#### drag calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 3  
def drag3(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = [1.07*2*S[i]*0.003 for i in range(3)]  
    l2 = [l_fuselage[i] - l_nose[i] - l_tailcone[i] for i in range(3)]
    Fuselage = [0.0024*(pi*D[i]/4)*(1/(3*l_nose[i]**2)*((4*l_nose[i]**2+D[i]**2/4)**(1.5)-D[i]**3/8)-D[i]+4*l2[i]+2*sqrt(l_tailcone[i]**2+D[i]**2/4)) for i in range(3)]
    Nacelle = []
    for i in range(3):
        if i == 2:
            Nacelle.append((2*(Dnacel/2)**2*pi + Dnacel*pi*Lnacel)*0.006 + (((45*inch_to_m)/2)**2*pi*2 + (45*inch_to_m*pi*308.4*inch_to_m))*0.006)
        else:
            Nacelle.append((2*(Dnacel/2)**2*pi + Dnacel*pi*Lnacel)*0.0060)
            
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/S[0]*(Wing[i] + Fuselage[i] + Nacelle[i] + Tailplane[i])*1.1  for i in range(3)]
    
    e = [4.61*(1-0.045*A[i]**0.68)*(cos(sweep_le[i]))**0.15 -3.1 for i in range(3)]
    
    CDcruise = [CD0[i] + 1/(pi*A[i]*e[i])*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD, Nacelle)