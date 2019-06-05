# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:46:51 2019

@author: Lisa
"""
import numpy as np
from inputs.concept_1 import b, S
from inputs.constants import *
flap_deflection_rad=radians(40)
#noise constants from ANOPP method
#G_wing=0.37*S/b**2*((rho*M_cruise*a*S)/(mu_0*b))**-0.2
#L_wing=G_wing*b
#K_wing=4.464E-5 
#a_wing=5
#
#G_slat=G_wing
#L_slat=L_wing
#K_slat=K_wing
#a_slat=a_wing
#
#G_flap=S_flap/b**2*(sin(flap_deflection_rad))**2
#L_flap=S_flap/b_flap
#K_flap=2.787*10**-4
#a_flap=6
#
#G_lg=N_mw*(d_mw/b)**2
#l_lg=d_mw
#K_lg=3.414E-4  #4 wheels
#K_nlg=4.349E-4 #1 or 2 wheels (nose)
#a_lg=6
#
#K_strut=2.753E-4 
##radius definitions
#r_lateral=450
#r_app=2000
#r_flyover=6500
    
#start with the  1/3 octave bands 
number_of_bands=43
bandnumbers=list(range(1,number_of_bands+1))
centrefreq=[10**(bandnumbers[i]/10) for i in range(len(bandnumbers))]
lowerfreq=[2**(-1/6)*centrefreq[i] for i in range(len(bandnumbers))]
upperfreq=[2**(1/6)*centrefreq[i] for i in range(len(bandnumbers))]
delta_freq=[upperfreq[i]-lowerfreq[i] for i in range(len(bandnumbers))]

