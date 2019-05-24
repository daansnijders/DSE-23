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
from modules.initialsizing_undercarriage import *
from modules.payload_range import *
# from modules.initialsizing_loading import *     # commented out because this import immediately runs the plot......
from inputs.constants import M_cruise, M_x, rho, V_cruise, N_sa, l_cockpit, inchsq_to_msq

N_pax = [90,120,120]                                                            # [-] number of passengers
R = [4000E3,2000E3,4000E3]                                                      # [m] range of the aircraft

T_W    = [0.295,0.295,0.295]                                                    # [-] thrust over weight ratio
W_S    = [4253, 4253, 4253]                                                     # [N/m^2] weight over wing surface area
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg] mass fuel fraction
OEW = [34631.92,38223.31,38729.81]                                              # [kg] operational empty weight
MTOW = [58722.6,67394,68264.27]                                                 # [kg] maximum take-off weight
M_fuel = get_M_fuel(MTOW,M_ff)                                                  # [kg] fuel mass
T_req = get_T_req(T_W, MTOW)                                                    # [N] required thrust
M_payload = get_M_payload_available(MTOW,OEW,M_fuel)                            # [kg] payload mass

d_OEW1,d_OEW2=get_mass_efficiency(OEW)
# Fuselage parameters
l_cabin = get_l_cabin(N_pax,N_sa)                                               # [m] cabin length

d_f_inner = get_d_f_inner(N_sa, seat_width, N_aisle,\
                          armrest, aisle_width, s_clearance)                    # [m] inner diameter fuselage

d_f_outer = get_d_f_outer(d_f_inner)                                            # [m] outer diameter fuselage
l_nose = get_l_nose(d_f_outer)                                                  # [m] nose length
l_tailcone = get_l_tailcone(d_f_outer)                                          # [m] tailcone length
l_tail = get_l_tail(d_f_outer)                                                  # [m] tail length
l_f = get_l_fuselage(l_cockpit, l_cabin, l_tail)                                # [m] length fuselage
R_f = [d_f_outer[i]/2 for i in range(3)]                                        # [m] radius fuselage

V_os = get_overhead_volume(l_cabin)                                             # [m^3] overhead storage volume
V_cc = get_cargo_volume(R_f,l_cabin)                                            # [m^3] total storage volume

Mtot_carry_on, Mtot_check_in, V_carry_on\
, V_check_in = get_masses_volumes(N_pax, V_cc, V_os)                            # [kg,kg,m^3,m^3] mass and volume of check-in/carry-on luggage

V_cargo_available = get_available_cargo_volume(V_cc,V_os,V_carry_on, V_check_in)# [m^3] available cargo volume
M_cargo_available = get_cargo_mass(N_pax,M_payload)                             # [kg] available cargo mass

#Propulsion
d_fan = 2.006                                                                   # [m] diameter of engine fan
d_nacel = 1.1*d_fan                                                             # [m] diameter of engine nacelle
l_eng = 3.184                                                                   # [m] length of the engine
l_nacel = 1.1*l_eng                                                             # [m] length of the engine nacelle

# Wing parameters
A = [11,11,11]                                                                  # [-] aspect ration main wing
e = [0.85,0.85,0.85]                                                            # [-]
S = get_S(MTOW,W_S)                                                             # [m^2] surface area main wing
b = get_b(A,S)                                                                  # [m] span main wing
lambda_4_rad = get_lambda_4_rad(M_cruise,M_x)                                   # [rad] quarter chord sweep angle main wing
taper_ratio = get_taper_ratio(lambda_4_rad)                                     # [-] taper ratio main wing
lambda_2_rad = get_lambda_2_rad(lambda_4_rad,A,taper_ratio)                     # [rad] half chord sweep angle main wing
Cr = get_Cr(S,taper_ratio,b)                                                    # [m] root chord length main wing
Ct = get_Ct(Cr, taper_ratio)                                                    # [m] tip chord length main wing
CL = get_CL(MTOW,rho,V_cruise,S)                                                # [-] lift coefficient aircraft
t_c =  get_t_c(lambda_2_rad,M_x, M_cruise,CL)                                   # [-] thickness over chord main wing
MAC = get_MAC(Cr, taper_ratio)                                                  # [m] mean aerodynamic chord main wing
y_MAC = get_y_MAC(b, Cr, MAC, Ct)                                               # [m] y-location of the MAC of the main wing
dihedral_rad = get_dihedral_rad(lambda_4_rad)                                   # [rad] dihedral angle of the main wing
lambda_le_rad = get_lambda_le_rad(lambda_4_rad, Cr, b, taper_ratio)             # [rad] leading edge sweep angle main wing

# Empennage parameters
V_h = [1.28, 1.28, 1.28]                                                        # [-] volume horizontal tail
A_h = [4.95, 4.95, 4.95]                                                        # [-] aspect ratio horizontal tail
taper_ratio_h = [0.39, 0.39, 0.39]                                              # [-] taper ratio horizontal tail
lambda_h_le = [np.deg2rad(34) for i in range(3)]                                # [rad] leading edge sweep angle horizontal tail
V_v = [0.1, 0.1, 0.1]                                                           # [-] volume vertical tail
A_v = [1.9, 1.9, 1.9]                                                           # [-] aspect ratio vertical tail
taper_ratio_v = [0.375, 0.375, 0.375]
                                      # [-] taper ratio vertical tail
lambda_v_le = [np.deg2rad(40) for i in range(3)]                                # [rad] leading edge sweep angle vertical tail
x_le_h = get_x_h(l_f)                                                           # [m] x-position leading edge horizontal tail
x_le_v = x_le_h                                                                 # [m] x-position leading edge vertical tail
x_cg = get_x_cg(l_f,MTOW, MAC)                                                  # [m] x-location of the centre of mass aircraft
y_cg = get_y_cg()                                                               # [m] y-location of the centre of mass aircraft
z_cg = get_z_cg(d_f_outer)                                                      # [m] z-location of the centre of mass aircraft
S_h = get_S_h(S, MAC, x_cg, V_h, x_le_h)                                        # [m^2] surface area horizontal tail
S_v = get_S_v(S, b, x_cg, V_v, x_le_v)                                          # [m^2] surface area vertical tail
b_h = get_b_h(S_h, A_h)                                                         # [m] span horizontal tail
b_v = get_b_v(S_v, A_v)                                                         # [m] span vertical tail
Cr_h = get_Cr_h(S_h, taper_ratio_h, b_h)                                        # [m] root chord length horizontal tail
Ct_h = get_Ct_h(Cr_h, taper_ratio_h)                                            # [m] tip chord length horizontal tail
Cr_v = get_Cr_v(S_v, taper_ratio_v, b_v)                                        # [m] root chord lengh vertical tail
Ct_v = get_Ct_v(Cr_v, taper_ratio_v)                                            # [m] tip chord length vertical tail


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

l_m = get_l_mw(x_mlg,x_cg)                                                      # [m] mlg distance from c.g
l_n = get_l_nw(l_m,P_mw,N_mw,P_nw,N_nw)                                         # [m] nlg distance from c.g

y_mlg = get_y_mlg(b,dihedral_rad,psi_rad,phi_rad,\
                  z_cg,z_mlg,l_n,l_m,y_eng,z_eng,d_eng)                         # [m] y-location of the mlg

x_nlg = get_x_nlg(x_cg,l_n)                                                     # [m] x-location of nlg
y_nlg = [0,0,0]                                                                 # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg


# Airfoil Cl,max from javafoil for Re = [9*10^6, 17*10^6, 20*10^6]

Cl_max = [1.552, 1.582, 1.584]

# airfoil design
# CLmax: Wing CL max for three Re numbers: [9*10^6, 17*10^6, 20*10^6]
# CL_alpha: Wing CL_alpha for three configurations
Reto1, Reto2, Reto3, CLdes, Cl_des, CL_alpha, CLmax, CLmaxto=airfoil(Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, lambda_le_rad, lambda_2_rad, b, taper_ratio, A, Cl_max)

CD0, CDcruise, LoverD, Wing, Fuselage, Nacelle, Tailplane=drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_f, d_f_outer, d_nacel, l_nacel, lambda_le_rad, CLdes)


# loadingdiagram=plot_loadingdiagram(Sland,CLmaxto,CLmax,CLmaxto,c,f,sigma, TOP, CD0,100,7100,100)

"""
---
~~~ Performance
---
"""
cg_loc = [[12.80353534, 12.76158237], [16.92946042, 17.3060118], [16.93525685, 17.46464004]]  # c.g. location [m],
#  because Daan uses excel

"""
Airport performance
"""

# take-off
take_off_thrust = 2*thrust_max
climb_out_thrust = 2*0.85*thrust_max
take_off_friction_coefficient = [get_friction_coefficient(P_nw[i], MTOW[i], x_mlg[i], x_nlg[i], cg_loc[i][1], z_cg[i] -
                                                          z_mlg[i], g) for i in range(3)]

take_off_field_length = [get_take_off_field_length(rho_0, g, h_screen, MTOW[i], take_off_thrust, climb_out_thrust,
                                                   CDcruise[i], CLmaxto[i], S[i], take_off_friction_coefficient[i])
                         for i in range(3)]

# landing
landing_thrust = 2*thrust_max  # for thrust reversal
landing_mass = [get_m_landing(MTOW[i], 2*thrust_max) for i in range(3)]
landing_friction_coefficient = [get_friction_coefficient(P_nw[i], landing_mass[i], x_mlg[i], x_nlg[i], cg_loc[i][0],
                                                         z_cg[i] - z_mlg[i], g) for i in range(3)]

landing_field_length = [get_landing_field_length(landing_thrust, landing_mass[i], g, h_screen, rho_0, S[i], CLmaxto[i],
                                                 CDcruise[i], landing_friction_coefficient[i])
                        for i in range(3)]

"""
Cruise fuel economy
"""
fuel_cruise = [get_cruise_fuel(get_cruise_thrust(rho_0, V_cruise, S[i], CDcruise[i]), R[i], V_cruise) for i in range(3)]

"""
Climb performance
"""
V_to = [1.05*get_V_min(MTOW[i], g, rho_0, S[i], CLmax[i]) for i in range (3)]  # horizontal velocity during TO climb
V_approach = [1.3*get_V_min(MTOW[i], g, rho_0, S[i], CLmax[i]) for i in range (3)]  # horizontal velocity during landing
CDto = drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_f, d_f_outer, d_nacel, l_nacel, lambda_le_rad, CLmaxto)[1]
climb_gradient = [get_climb_gradient(2*thrust_max, 0.5 * rho_0 * V_approach[i]**2 * CDto[i] * S[i], MTOW[i], g)
                  for i in range(3)]
rate_of_climb = [get_rate_of_climb(V_to[i], climb_gradient[i]) for i in range(3)]

"""
Mass/payload-range diagram
"""
# [generate_payload_range_diagram(M_payload[i], M_fuel[i], MTOW[i], R[i], V_cruise, 0.5*2.832545035E-5, 14, g, OEW[i], i)
#  for i in range(3)]
