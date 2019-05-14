# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
from modules.initialsizing_weights import *
from modules.initialsizing_planform import *
from inputs.constants import M_cruise, M_x, rho, V_cruise,p,R,T
import numpy as np



T_W    = [0.295,0.295,0.295]                                                    # [-]
W_S    = [4253, 4253, 4253]                                                     # [N/m^2]
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg]
OEW = [27745.73, 27069.62, 38729.81]                                            # [-]

MTOW = get_MTOW(OEW)                                                            # [kg]
M_fuel = get_M_fuel(MTOW,M_ff)                                                  # [kg]
T_req = get_T_req(T_W, MTOW)                                                    # [N]
M_payload = get_M_payload(MTOW,OEW,M_fuel)                                      # [kg]


# Fuselage parameters


# Wing parameters
A = [11,11,11]                                                                  # [-]
S = get_S(MTOW,W_S)                                                             # [m^2]
b = get_b(A,S)                                                                  # [m]
lambda_4_rad = get_lambda_4_rad(M_cruise,M_x)                                   # [rad]
taper_ratio = get_taper_ratio(lambda_4_rad)                                     # [-]
lambda_2_rad = get_lambda_2_rad(lambda_4_rad,A,taper_ratio)                     # [rad]
Cr = get_Cr(S,taper_ratio,b)                                                    # [m]
Ct = get_Ct(Cr, taper_ratio)                                                    # [m]
CL = get_CL(MTOW,rho,V_cruise,S)                                                # [-]
t_c =  get_t_c(lambda_2_rad,M_x, M_cruise,CL)                                   # [-]
MAC = get_MAC(Cr, taper_ratio)                                                  # [m]
y_MAC = get_y_MAC(b, Cr, MAC, Ct)                                               # [m]
dihedral_rad = get_dihedral_rad(lambda_4_rad)                                   # [rad]


# Propulsion parameters
