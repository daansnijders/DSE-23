# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:41:25 2019

@author: Lisa
"""
from modules.sustainability.noise_defs import *
from modules.sustainability.noise_calc import *

pe_2_flap =get_effective_pressure_flap(rho_0,a_sl,M,r_observer,theta,phi,L,K,a_const,G,flap_deflection)

pe_2_slat = get_effective_pressure_slat(rho_0,a_sl,M,r_observer,theta,phi,L,K,a_const,G)

pe_2_wing = get_effective_pressure_wing(rho_0,a_sl,M,r_observer,theta,phi,L,K,a_const,G)

pe_2_lg_main =  get_effective_pressure_lg(rho_0,a_sl,M,r_observer,theta,phi,d,K,N_w,d_w)

pe_2_lg_noise =  get_effective_pressure_lg(rho_0,a_sl,M,r_observer,theta,phi,d,K,N_w,d_w)

pe_2_strut= get_effective_pressure_strut(rho_0,a_sl,M,r_observer,theta,phi,d,K,d_strut,L_strut)
