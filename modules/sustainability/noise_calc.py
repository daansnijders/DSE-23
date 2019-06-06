# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:46:51 2019

@author: Lisa
"""
import numpy as np
from inputs.concept_1 import b, S, V_TO, N_mw, N_nw, z_mlg
from inputs.constants import mu_sl, rho_0,M_cruise

d_mw=1.5 #[m]
d_nw=0.7    #[m]
flap_deflection_rad=np.radians(40)
L_strut= -z_mlg
diameter_strut = 0.3

#noise constants from ANOPP method

G_wing=0.37*S/b**2*((rho_0*M_cruise*a*S)/(mu_sl*b))**-0.2
L_wing=G_wing*b
K_wing=4.464E-5 
a_wing=5

G_slat=G_wing
L_slat=L_wing
K_slat=K_wing
a_slat=a_wing

flap_deflection_rad=radians(40)
b_flap=b/2*0.6
c_flap=0.35*Cr
S_flap=b_flap*c_flap

G_flap=S_flap/b**2*(sin(flap_deflection_rad))**2
L_flap=S_flap/b_flap
K_flap=2.787*10**-4
a_flap=6

#G_lg_m=N_mw*(d_mw/b)**2
#G_lg_n=N_nw*(d_nw/b)**2

K_lg=3.414E-4  #4 wheels
K_nlg=4.349E-4 #1 or 2 wheels (nose)
K_strut=2.753E-4 



'CONDITIONS WHERE WE EVALUATE THE NOISE'
#radius definitions
r_lateral=450
r_app=2000
r_flyover=6500
    
#get r and phi and theta for the observers





#start with the  1/3 octave bands 
number_of_bands=43
bandnumbers=list(range(1,number_of_bands+1))
centrefreq=[10**(bandnumbers[i]/10) for i in range(len(bandnumbers))]
lowerfreq=[2**(-1/6)*centrefreq[i] for i in range(len(bandnumbers))]
upperfreq=[2**(1/6)*centrefreq[i] for i in range(len(bandnumbers))]
delta_freq=[upperfreq[i]-lowerfreq[i] for i in range(len(bandnumbers))]

