# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:41:25 2019

@author: Lisa
"""
from modules.sustainability.noise_defs import *

phi_observer=radians(1)
flap_deflection=np.radians(40)


b_flap=b/2*0.4
c_flap=0.35*Cr
S_flap=b_flap*c_flap

V_approach=64

r1,r2,r3,theta_1,theta_2,theta_3=simulate_flight_path(V_approach)

OSPL_dBA_tot_2=EPNdB_calculations(r2,theta_2,phi_observer,V_approach, S_flap, b_flap,flap_deflection )
OSPL_dBA_tot_1=EPNdB_calculations(r1,theta_1,phi_observer,V_approach, S_flap, b_flap,flap_deflection )
OSPL_dBA_tot_3=EPNdB_calculations(r3,theta_3,phi_observer,V_approach, S_flap, b_flap,flap_deflection )

centrefreq,freq_delta=get_octave_frequency_bands()
