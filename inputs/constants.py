# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:27:30 2019

@author: Lisa
"""

from math import * 

#conversion
inch_to_m=0.0254
ft_to_m=0.3048


#flight parameters
H_ft=37000   #[feet]
H_m=H_ft*ft_to_m      #[m]
M_cr=0.75
M_x = 0.935

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
    
    

