# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:01:38 2019

@author: Anique
"""
from math import *

#### Inputs needed for the definition, they are from requirements or from previous sizing ####
#### This is an example for concept 1, configuration 1 ####

h = 37000            # Altitude [ft]
rho = 0.348331       # Density [kg/m^3]
T = 216.65           # Temperature [K]
a = 295.07           # Speed of Sound [m/s]
mu = 0.0000143226    # Dynamic Viscosity [Pa*s]
M = 0.75             # Mach number [-]
Cr = 4.91489505      # Root chord [m]
Ct = 1.520444671     # Tip chord [m]
MTOW = 484363.0809   # Maximum Take-off Weight [N]
FF1 = 0.99           # Fuel Fraction phase 1 [-]
FF2 = 0.99           # Fuel Fraction phase 2 [-]
FF3 = 0.995          # Fuel Fraction phase 3 [-]
FF4 = 0.98           # Fuel Fraction phase 4 [-]
FF5 = 0.8366         # Fuel Fraction phase 5 [-]
S = 113.8873926      # Surface Area [m^2]
Sw25 = 25.96803627   # Quarterchord sweep [deg]
b = 35.39436847      # Wing span [m]
Taper = 0.3093544534 # Taper ratio [-]


#### Definition calculates different Reynold numbers for different locations 
#### and calculates Cl design for the airfoils
def airfoil(h, rho, T, a, mu, M, Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, Sw25, b, Taper):
    V = M*a                 # Cruise velocity [m/s]
    avgC = (Cr + Ct)/2      # Average chord [m]
    Re1 = (rho*V*Cr)/mu     # Reynolds number at root chord [-]
    Re2 = (rho*V*avgC)/mu   # Reynolds number at avg chord [-]
    Re3 = (rho*V*Ct)/mu     # Reynolds number at tip chord [-]
    q = 0.5*rho*V**2        # Dynamic pressure [Pa]
    WSbegin = (MTOW*FF1*FF2*FF3*FF4)/S          # Wing Loading begin cruise  
    WSend = (MTOW*FF1*FF2*FF3*FF4*FF5)/S        # Wing loading end cruise
    CLdes = 1.1*(1/q)*(0.5*(WSbegin + WSend))   # Wing CL design [-]
    Sweep = atan(tan((pi/180)*Sw25) - ((Cr/(2*b))*(Taper - 1)))*(180/pi)
                                                # Sweep angle [deg]
    Veff = V*cos((pi/180)*Sweep)                # Effective velocity [m/s]
    qeff = 0.5*rho*Veff**2                      # Effective Dynamic Pressure [Pa]
    Cl_des = (q*CLdes)/qeff                     # Airfoil Cl design [-]
    
    return(Re1, Re2, Re3, Cl_des)
    
print(airfoil(h, rho, T, a, mu, M, Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, Sw25, b, Taper))
