# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:02:42 2019

@author: Lisa
"""
from inputs.constants import *
#lift coefficient
#Input in CL from Roskam 
#CL_clean_min=1.2
#CL_clean_max=1.8
Cl_clean=1.5

#CL_TO_min=1.6
#CL_TO_max=2.2
Cl_TO=1.9

#CL_land_min=1.8
Cl_land=2.8
#CL_land_av=2.8


#take off and landing performance
Sland=1500/ft_to_m           #field length in feet
Stakeoff=2000/ft_to_m        # field length in feet
TOP=175*4.8824*9.81         #N/M^2
sigma=1
h_screen = 50 * ft_to_m         # [m] Screen height

Vto1 = 135  #Initial guess Vtake-off                                            # [kts]
VL1 = 172   #Initial guess Vtake-off                                            # [kts]

#statistical input for cd0
f=0.84
Swet_S=6
Cfe=0.0045

#statistical input from roskam
d_CD0_to=0.015
d_CD0_landing=0.065
d_CD0_gear=0.02

#assumed/chosen values
c=15

#fuel fractions
FF1 = 0.99           # Fuel Fraction phase 1 [-]
FF2 = 0.99           # Fuel Fraction phase 2 [-]
FF3 = 0.995          # Fuel Fraction phase 3 [-]
FF4 = 0.98           # Fuel Fraction phase 4 [-]
FF5 = 0.8366         # Fuel Fraction phase 5 [-]

#engine characteristics
thrust_max = 108.5E3        # [N] maximum engine thrust (1 engine) a.k.a. rated output
