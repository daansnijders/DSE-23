# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:02:42 2019

@author: Lisa
"""

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

#take off and landing
Sland=1500/ft_to_m           #field length in feet
Stakeoff=2000/ft_to_m        # field length in feet
TOP=175*4.8824*9.81         #N/M^2
sigma=1

#statistical input
f=0.84
Swet_S=6
Cfe=0.0045

d_CD0_to=0.015
d_CD0_landing=0.065
d_CD0_gear=0.02

#assumed/chosen values
A=11
e=0.85
c=15