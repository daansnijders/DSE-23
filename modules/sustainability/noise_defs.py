# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:53:00 2019

@author: Lisa
"""
from math import *
#noise definitions

#f will be a list which will be an 1./third octave band or complete octave band 
def get_strouhal_number(f,L,M,theta,a):
    S=f*L*(1-M*cos(theta))/(M*a)
    return S

def get_strouhal_number_lg_and_strut(f,d,V):
    S=f*d/V # normal Mach times speed of sound which is V uses both the 
    return S

def get_power(K,M,a_const,G,rho,a,b):
    P=K*M**a_const*G*(rho*a**3*b**2)
    return P
    
def get_spectrical_function_wing(S):
    F=0.613*(10*S)**4*((10*S)**1.5+0.5)**-4
    return F
def get_spectrical_function_flap(S):
    if S<2:
        F=0.048*S
    if S>=2 and S<=20:
        F=0.1406*S**-0.55
    else:
        F=216.49*S**-3
    return F

def get_spectrical_function_slat(S):
    F=0.613*(10*S)**4*((10*S)**1.5+0.5)**-4 + 0.613*(2.19*S)**4*((2.19*S)**1.5+0.5)**-4
    return F

def get_spectrical_function_landinggear(S): #this is for two wheels
    F=2*(13.59*S**2*(S**2+12.5)**-2.25)
    return F

def get_directivity_function_wing(theta,phi):  #same for slat
    D=4*(cos(phi))**2*(cos(theta/2))**2
    return D

def get_directivity_function_flap(flap_deflection_rad,theta,phi):
    D=3*(sin(flap_deflection_rad)*cos(theta)+cos(flap_deflection_rad)*sin(theta)*cos(phi))**2
    return D

def get_directivity_function_landinggear(theta):
    D=3/2*(sin(theta))**2
    return D

def get_effective_pressure(rho,a,P,D,F,r_observer,M,theta):
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_sound_pressure_level(p_e_squared):   #[dB]
    SPL=10*log10(p_e_squared/(2E-5)**2)
    return SPL

def correction_for_thirdbandwidth(SPL,freq_delta):
    SPL_bandwidth=SPL+10*log10(freq_delta)
    return SPL_bandwidth

def get_perceived_noise_level(N):
    L_pn=40+33.3*log10(N)
    return L_pn


def get_effective_perceived_noise_level(L_pn,delta_t):
    L_epn=10*log10(delta_t/10*sum(10**(L_pn/10)))

def get_overall_sound_pressure_level(SPL,centrefreq):
    OSPL=10*log10(sum(10**(SPL/10)))
    
    
    
    
    
    
    