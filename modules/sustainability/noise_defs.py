# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:53:00 2019

@author: Lisa
"""
from math import *
#noise definitions

def get_strouhal_number(f,L,M,theta,a):
    S=f*L*(1-M*cos(theta))/(M*a)
    return S

def get_power(K,M,a,G,rho,a,b):
    P=K*M**a*G*(rho*a**3*b**2)
    return P
    
def get_spectrical_function(S):
    return F

def get_directivity_function(theta,phi):
    return D

def get_effective_pressure(rho,a,P,D,F,r_observer,M,theta):
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_sound_pressure_level(p_e_squared):   #[dB]
    SPL=10*log(p_e_squared/(2E-5)**2)

#noise constants
G_wing=0.37*S/b**2*((rho*M*a*S)/(mu*b))**-0.2
L_wing=G*b
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

G_lg=N_mw*(d_mw/b)**2
l_lg=d_mw
K_lg=3.414E-4  #4 wheels
K_nlg=4.349E-4 #1 or 2 wheels (nose)
a_lg=6

#radius definitions
r_lateral=450
r_app=2000
r_flyover=6500
