# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
from modules.initialsizing_weights import *
from modules.initialsizing_planform import *
from modules.initialsizing_fuselage import *
from modules.initialsizing_empennage import *
from modules.initialsizing_cg import *
from modules.airfoil_calculations import *
from inputs.constants import M_cruise, M_x, rho, V_cruise, N_sa, l_cockpit
import numpy as np

N_pax = [90,120,120]
R = [4000,2000,4000]

T_W    = [0.295,0.295,0.295]                                                    # [-]
W_S    = [4253, 4253, 4253]                                                     # [N/m^2]
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg]
OEW = [27745.73, 27069.62, 38729.81]                                            # [-]

MTOW = get_MTOW(OEW)                                                            # [kg]
M_fuel = get_M_fuel(MTOW,M_ff)                                                  # [kg]
T_req = get_T_req(T_W, MTOW)                                                    # [N]
M_payload = get_M_payload(MTOW,OEW,M_fuel)                                      # [kg]


# Fuselage parameters
l_cabin = get_l_cabin(N_pax,N_sa)                                               # [m]
d_f_inner = get_d_f_inner(N_sa, seat_width, N_aisle,\
                          armrest, aisle_width, s_clearance)                    # [m]
d_f_outer = get_d_f_outer(d_f_inner)                                            # [m]
l_nose = get_l_nose(d_f_outer)                                                  # [m]
l_tailcone = get_l_tailcone(d_f_outer)                                          # [m]
l_tail = get_l_tail(d_f_outer)                                                  # [m]
l_f = get_l_fuselage(l_cockpit, l_cabin, l_tail)                                # [m]

R_f=[d_f_outer[i]/2 for i in range(3)] 

V_os= get_overhead_volume(l_cabin)
V_cc=get_cargo_volume(R_f,l_cabin)
M_cargo_available=get_cargo_mass(N_pax,V_cc, V_os)
M_payload_total=get_payload_mass(M_cargo_available,N_pax,V_cc,V_os)


# Wing parameters
A = [11,11,11]   
e=[0.85,0.85,0.85]                                                              # [-]
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

# Empennage parameters
V_h = [1.28, 1.28, 1.28]                                                        # [-]
A_h = [4.95, 4.95, 4.95]                                                        # [-]
taper_ratio_h = [0.39, 0.39, 0.39]                                              # [-]
lambda_h_le = [np.deg2rad(34) for i in range(3)]                                # [rad]

V_v = [0.1, 0.1, 0.1]                                                           # [-]
A_v = [1.9, 1.9, 1.9]                                                           # [-]
taper_ratio_v = [0.375, 0.375, 0.375]                                           # [-]
lambda_v_le = [np.deg2rad(40) for i in range(3)]                                # [rad]

x_h = get_x_h(l_f)                                                              # [m]
x_v = x_h                                                                       # [m]

x_cg = get_x_cg(l_f,MTOW, MAC)                                                  # [m]
y_cg = get_y_cg()                                                               # [m]
z_cg = get_z_cg(d_f_outer)                                                      # [m]

S_h = get_S_h(S, MAC, x_cg, V_h, x_h)                                           # [m^2]
S_v = get_S_v(S, MAC, x_cg, V_v, x_v)                                           # [m^2]

b_h = get_b_h(S_h, A_h)                                                         # [m]          
b_v = get_b_v(S_v, A_v)                                                         # [m]
Cr_h = get_Cr_h(S_h, taper_ratio_h, b_h)                                        # [m]
Ct_h = get_Ct_h(Cr_h, taper_ratio_h)                                            # [m]
Cr_v = get_Cr_v(S_v, taper_ratio_v, b_v)                                        # [m]
Ct_v = get_Ct_v(Cr_v, taper_ratio_v)                                            # [m]



#airfoil design 

Re1, Re2, Re3, Cl_des=airfoil(Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, lambda_2_rad, b, taper_ratio)