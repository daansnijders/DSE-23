# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt
from modules.performance import *
from modules.initialsizing_cg import *
from modules.airfoil_calculations import *
from modules.initialsizing_weights import *
from modules.initialsizing_planform import *
from modules.initialsizing_fuselage import *
from modules.initialsizing_empennage import *
from modules.initialsizing_undercarriage import*
from inputs.constants import M_cruise, M_x, rho, V_cruise, N_sa, l_cockpit, inchsq_to_msq

N_pax = [90,120,120]                                                            # [-]
R = [4000E3,2000E3,4000E3]                                                      # [m]

T_W    = [0.295,0.295,0.295]                                                    # [-]
W_S    = [4253, 4253, 4253]                                                     # [N/m^2]
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg]
OEW = [36264.04, 38789.59, 40108.34]                                            # [-]
MTOW = [61064.99,62883.99,70632.89]#get_MTOW(OEW)                   #make this an input or change the definition                                         # [kg]
M_fuel = get_M_fuel(MTOW,M_ff)                                                  # [kg]
T_req = get_T_req(T_W, MTOW)                                                    # [N]
M_payload = get_M_payload_available(MTOW,OEW,M_fuel)                                      # [kg]


# Fuselage parameters
l_cabin = get_l_cabin(N_pax,N_sa)                                               # [m]

d_f_inner = get_d_f_inner(N_sa, seat_width, N_aisle,\
                          armrest, aisle_width, s_clearance)                    # [m]    

d_f_outer = get_d_f_outer(d_f_inner)                                            # [m]
l_nose = get_l_nose(d_f_outer)                                                  # [m]
l_tailcone = get_l_tailcone(d_f_outer)                                          # [m]
l_tail = get_l_tail(d_f_outer)                                                  # [m]
l_f = get_l_fuselage(l_cockpit, l_cabin, l_tail)                                # [m]


R_f = [d_f_outer[i]/2 for i in range(3)]                                        # [m]

V_os = get_overhead_volume(l_cabin)                                             # [m^3]
V_cc = get_cargo_volume(R_f,l_cabin)                                            # [m^3]

Mtot_carry_on, Mtot_check_in, V_carry_on\
, V_check_in = get_masses_volumes(N_pax, V_cc, V_os)                            # [kg,kg,m^3,m^3]

V_cargo_available = get_available_cargo_volume(V_cc,V_os,V_carry_on, V_check_in)# [m^3]
M_cargo_available = get_cargo_mass(N_pax,M_payload)                             # [kg]

#Propulsion
Dfan = 2.006                                                                    # [m]
Dnacel = 1.1*Dfan                                                               # [m]
Lfan = 3.184                                                                    # [m]
Lnacel = 1.1*Lfan                                                               # [m]

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
lambda_le_rad = get_lambda_le_rad(lambda_4_rad, Cr, b, taper_ratio)             # [rad]


# Empennage parameters
V_h = [1.28, 1.28, 1.28]                                                        # [-]
A_h = [4.95, 4.95, 4.95]                                                        # [-]
taper_ratio_h = [0.39, 0.39, 0.39]                                              # [-]
lambda_h_le = [np.deg2rad(34) for i in range(3)]                                # [rad]

V_v = [0.1, 0.1, 0.1]                                                           # [-]
A_v = [1.9, 1.9, 1.9]                                                           # [-]
taper_ratio_v = [0.375, 0.375, 0.375]                                           # [-]
lambda_v_le = [np.deg2rad(40) for i in range(3)]                                # [rad]

x_le_h = get_x_h(l_f)                                                           # [m]
x_le_v = x_le_h                                                                 # [m]

x_cg = get_x_cg(l_f,MTOW, MAC, concept_3 = True)                                # [m]
y_cg = get_y_cg()                                                               # [m]
z_cg = get_z_cg(d_f_outer)                                                      # [m]

S_h = get_S_h(S, MAC, x_cg, V_h, x_le_h)                                        # [m^2]
S_v = get_S_v(S, b, x_cg, V_v, x_le_v)                                          # [m^2]

b_h = get_b_h(S_h, A_h)                                                         # [m]          
b_v = get_b_v(S_v, A_v)                                                         # [m]
Cr_h = get_Cr_h(S_h, taper_ratio_h, b_h)                                        # [m]
Ct_h = get_Ct_h(Cr_h, taper_ratio_h)                                            # [m]
Cr_v = get_Cr_v(S_v, taper_ratio_v, b_v)                                        # [m]
Ct_v = get_Ct_v(Cr_v, taper_ratio_v)                                            # [m]


# Undercarriage
"""Inputs that might be better located in constants"""
N_mw = 4                                                                        # [-] number of wheels mlg
N_nw = 2                                                                        # [-] number of wheels nlg
N_struts = 2                                                                    # [-] number of struts used
stroke = 0.3                                                                    # [m] shock absorber stroke

LCN = 45                                                                        # [-] load classification number
tire_pressure = 430 * np.log(LCN) - 680                                         # [Pa] tire pressure mlg

weight_distribution = 0.08                                                      # [-] weight percentage on nose wheel
y_eng = [0.3*b[i]/2 for i in range(3)]
d_eng = 2.006                                                                   # [m] diameter of the engine
z_eng = -d_eng/2                                                                # [m] z-location of lowest part of the engine

theta = 15                                                                      # [deg] scrape angle
beta = 17                                                                       # [deg] tip-back angle
phi = 5                                                                         # [deg] tip clearance angle
psi = 55                                                                        # [deg] overturn angle
theta_rad = np.deg2rad(theta)                                                   # [rad] scrape angle
beta_rad = np.deg2rad(beta)                                                     # [rad] tip-back angle
phi_rad = np.deg2rad(phi)                                                       # [rad] tip clearance angle
psi_rad = np.deg2rad(psi)                                                       # [rad] overturn angle

P_mw = get_P_mw(MTOW,N_mw,weight_distribution)                                  # [N] static loading on mw
P_nw = get_P_nw(MTOW,N_nw,weight_distribution)                                  # [N] static loading on nw

x_mlg = get_x_mlg(z_cg,theta_rad,beta_rad, x_cg, stroke,l_f)                    # [m] x-location of the mlg
z_mlg = get_z_mlg(x_mlg,beta_rad,x_cg, z_cg, l_f)                               # [m] z-location of the mlg

l_w = get_l_mw(x_mlg,x_cg)                                                      # [m] mlg distance from c.g
l_n = get_l_nw(l_w,P_mw,N_mw,P_nw,N_nw)                                         # [m] nlg distance from c.g

y_mlg = get_y_mlg(b,dihedral_rad,psi_rad,phi_rad,\
                  z_cg,z_mlg,l_n,l_w,y_eng,z_eng,d_eng)                         # [m] y-location of the mlg

x_nlg = get_x_nlg(x_cg,l_n)                                                     # [m] x-location of nlg
y_nlg = [0,0,0]                                                                 # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg


#Airfoil Cl,max from javafoil for Re = [9*10^6, 17*10^6, 20*10^6]
Cl_max = [1.552, 1.582, 1.584]

#airfoil design 
# CLmax: Wing CL max for three Re numbers: [9*10^6, 17*10^6, 20*10^6]
# CL_alpha: Wing CL_alpha for three configurations
Reto1, Reto2, Reto3, CLdes, Cl_des, CL_alpha, CLmax, CLmaxto=airfoil(Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, lambda_le_rad, lambda_2_rad, b, taper_ratio, A, Cl_max)

CD0, CDcruise, LoverD=drag3(A, S, S_h, S_v, l_nose, l_tailcone, l_f, d_f_outer, Dnacel, Lnacel, lambda_le_rad, CLdes)
