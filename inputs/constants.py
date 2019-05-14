# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:27:30 2019

@author: Lisa
"""

from math import * 

#conversion
inch_to_m=0.0254
ft_to_m=0.3048


#FLIGHT PARAMETERS parameters
H_ft=37000   #[feet]
H_m=H_ft*ft_to_m      #[m]
M_cruise=0.75
M_x=0.935

#air constants
gamma=1.4
R=287.05
T_0=288.15

p_0=101325
g=9.80665
rho_0=1.225

if H_m<11000:

    T=T_0-0.0065*(H_m)

    p=p_0*(T/T_0)**(-g/(R*-0.0065))

    rho= rho_0*(T/T_0)**(-g/(R*-0.0065)-1)

elif H_m<20000:
    T1=216.65
    a1=0
    p1=22632


    T=T1+a1*(H_m)
    p=p1*e**(-g/(R*T)*(H_m-11000))
    rho= rho_0*(T1/T_0)**(-g/(R*-0.0065)-1)*p/p1
    


#get cruise velocities and speed of sound
a=(gamma*R*T)**0.5
V_cruise=M_cruise*a

#fuselage interior design 
#data on the interior
N_sa = 5                        # Number of seats abreast
N_aisle = 1                     # Number of aisles
aisle_width = 0.51              #[m]
s_clearance = 0.02              #[m] Clearance between seat and fuselage
seat_pitch = 32*inch_to_m       #[m]
seat_width = 20*inch_to_m       #[m]
armrest = 2*inch_to_m           #[m] Armrest width


