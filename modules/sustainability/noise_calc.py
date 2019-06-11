# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:46:51 2019

@author: Lisa
"""
import numpy as np
from math import *
from modules.sustainability.noise_defs import *


req_approach_EPNdB=89.9
req_lateral_TO_EPNdB=86
experimental_dBA_to_EPNdB=1.15503

M_TO=V_TO[2]/a_sl


d_mw=1.5 #[m]
d_nw=0.7    #[m]
flap_deflection=np.radians(40)
L_strut_main= -z_mlg-d_mw/2
L_strut_nose=  -z_mlg-d_nw/2

b_flap=b/2*0.4
c_flap=0.35*Cr
S_flap=b_flap*c_flap

#noise constants from ANOPP method

G_wing=0.37*S/b**2*((rho_0*M_TO*a_sl*S)/(mu_sl*b))**-0.2
L_wing=G_wing*b
K_wing=4.464E-5 
a_wing=5

G_slat=G_wing
L_slat=L_wing
K_slat=K_wing
a_slat=a_wing

G_flap=S_flap/b**2*(sin(flap_deflection))**2
L_flap=S_flap/b_flap
K_flap=2.787*10**-4
a_flap=6

G_mlg=N_mw*(d_mw/b)**2
G_nlg=N_nw*(d_nw/b)**2
K_mlg=3.414E-4  #4 wheels
K_nlg=4.349E-4 #1 or 2 wheels (nose)
K_strut=2.753E-4 
a_lg=6





#start with the  1/3 octave bands 
number_of_bands=43
bandnumbers=list(range(1,number_of_bands+1))
centrefreq=[10**(bandnumbers[i]/10) for i in range(len(bandnumbers))]
lowerfreq=[2**(-1/6)*centrefreq[i] for i in range(len(bandnumbers))]
upperfreq=[2**(1/6)*centrefreq[i] for i in range(len(bandnumbers))]
freq_delta=[upperfreq[i]-lowerfreq[i] for i in range(len(bandnumbers))]

