# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:41:25 2019

@author: Lisa
"""
from modules.sustainability.noise_defs import *
import inputs.concept_1 as c1

phi_observer=radians(0.1)
flap_deflection=radians(40)

b_slat=c1.b*0.6
b_flap=c1.b/2*0.4
c_flap=0.35*c1.Cr
area_flap=b_flap*c_flap

V_approach=75

r1,r2,r3,theta_1,theta_2,theta_3=simulate_flight_path(V_approach)


OSPL_dBA_tot_2=EPNdB_calculations(r2,theta_2,phi_observer,V_approach, area_flap, b_flap,flap_deflection,b_slat )
#OSPL_dBA_tot_1=EPNdB_calculations(r1,theta_1,phi_observer,V_approach, area_flap, b_flap,flap_deflection,b_slat )
#OSPL_dBA_tot_3=EPNdB_calculations(r3,theta_3,phi_observer,V_approach, area_flap, b_flap,flap_deflection,b_slat )


