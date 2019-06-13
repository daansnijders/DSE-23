# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:41:25 2019

@author: Lisa
"""
from modules.sustainability.noise_defs import *
r_observer=120
theta_observer=radians(90)
phi_observer=radians(1)
flap_deflection=np.radians(40)


b_flap=b/2*0.4
c_flap=0.35*Cr
S_flap=b_flap*c_flap

V_approach=V_TO[2]
r1,r2,r3,theta_1,theta_2,theta_3=simulate_flight_path(V_approach)

OSPL_dBA_flap, OSPL_dBA_slat,OSPL_dBA_wing,OSPL_dBA_mlg , OSPL_dBA_nlg ,OSPL_dBA_mlg_strut,OSPL_dBA_nlg_strut, OSPL_dBA_tot=EPNdB_calculations(r2,theta_2,phi_observer,V_approach, S_flap, b_flap,flap_deflection )