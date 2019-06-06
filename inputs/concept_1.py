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
from modules.initialsizing_loading import *     # commented out because this import immediately runs the plot......
from inputs.constants import *


 
#should move to constants
N_pax = [90,120,120]                                                            # [-] number of passengers
R = [4000E3,2000E3,4000E3]                                                      # [m] range of the aircraft
#inputs to this file 
T_W    = [0.29,0.29,0.29]                                                       # [-] thrust over weight ratio
W_S    = [4405, 4405 , 4405]                                                    # [N/m^2] weight over wing surface area
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg] mass fuel fraction
OEW = [34631.92,38223.31-360,38729.81]                                          # [kg] operational empty weight
#MTOW = [58722.6,67394-360,68264.27]                                             # [kg] maximum take-off weight
d_OEW1,d_OEW2=get_mass_efficiency(OEW)




#START SIZING 
# Fuselage parameters]
l_cutout=30/N_sa*seat_pitch + seat_pitch                                        #change with additional safety factors

l_cabin = get_l_cabin(N_pax,N_sa)                                               # [m] cabin length UPDATE THIS TO THE REAL VALUE
l_cabin = [max(l_cabin)-l_cutout, max(l_cabin), max(l_cabin)]
d_f_inner = get_d_f_inner(N_sa, seat_width, N_aisle,\
                          armrest, aisle_width, s_clearance)                    # [m] inner diameter fuselage
d_f_outer = get_d_f_outer(d_f_inner)                                            # [m] outer diameter fuselage
l_nose = get_l_nose(d_f_outer)                                                  # [m] nose length
l_tailcone = get_l_tailcone(d_f_outer)                                          # [m] tailcone length
l_tail = get_l_tail(d_f_outer)                                                  # [m] tail length
l_f = get_l_fuselage(l_cockpit, l_cabin, l_tail)                                # [m] length fuselage
R_f = d_f_outer/2                                                               # [m] radius fuselage
S_fus=[R_f*2*pi*l_f[i] for i in range(3)]                              # [m^2] gross shell fuselage area

V_os = get_overhead_volume(l_cabin)                                             # [m^3] overhead storage volume
V_cc = get_cargo_volume(R_f,l_cabin)                                            # [m^3] total storage volume

Mtot_carry_on, Mtot_check_in, V_carry_on\
, V_check_in = get_masses_volumes(N_pax, V_cc, V_os)                            # [kg,kg,m^3,m^3] mass and volume of check-in/carry-on luggage

V_cargo_available = get_available_cargo_volume(V_cc,V_os,V_carry_on, V_check_in)# [m^3] available cargo volume

#CALCULATE MASSES BASED ON THE FUSALGE LAYOUT
M_pax_and_lugg=get_passenger_luggage_mass(N_pax)
M_cargo_available=[V_cargo_available[i]*rho_cargo for i in range(3)]             # [kg] available cargo mass
M_payload=[M_cargo_available[i]+M_pax_and_lugg[i] for i in range(3)]
MTOW=get_TOW(OEW,M_payload,M_ff)
M_fuel = get_M_fuel(MTOW,M_ff)                                                  # [kg] fuel mass
T_req = get_T_req(T_W, MTOW)

#needed for class2 estimation
M_MZF    = [MTOW[i]-M_fuel[i] for i in range(3)]
M_carried_canard_MZF=[M_MZF[i]-M_MZF[0] for i in range(3)]
M_carried_canard_MTOW=[MTOW[i]-MTOW[0] for i in range(3)]                                                    # [N] required thrust                          

"Change this when correct length of modular part is found, implement the l_cutout here" 
x_cargo = [[Xcargo1, Xcargo2], [Xcargo1+(l_cabin[1]-l_cabin[0]), Xcargo2+(l_cabin[1]-l_cabin[0])]\
           , [Xcargo1+(l_cabin[1]-l_cabin[0]), Xcargo2+(l_cabin[1]-l_cabin[0])]]     



#Propulsion
d_nacel = 1.1*d_fan                                                             # [m] diameter of engine nacelle
l_nacel = 1.1*l_eng                                                             # [m] length of the engine nacelle

# Wing parameters MOVE TO CONSTANTS
A = 9.5                                                                         # [-] aspect ration main wing
e = 0.85                                                                        # [-] 
S = get_S(MTOW,W_S)                                                             # [m^2] surface area main wing
#take canard into account
S_c = [S[0]-S[0],S[1]-S[0],S[2]-S[0]]
S = min(S)                                                                      # [m] update surface area to be the same for all config
b = get_b(A,S)                                                                  # [m] span main wing
lambda_4_rad = get_lambda_4_rad(M_cruise,M_x)                                   # [rad] quarter chord sweep angle main wing
taper_ratio = get_taper_ratio(lambda_4_rad)                                     # [-] taper ratio main wing
lambda_2_rad = get_lambda_2_rad(lambda_4_rad,A,taper_ratio)                     # [rad] half chord sweep angle main wing
Cr = get_Cr(S,taper_ratio,b)                                                    # [m] root chord length main wing
Ct = get_Ct(Cr, taper_ratio)                                                    # [m] tip chord length main wing
CL = get_CL(MTOW,rho,V_cruise,S)                                                # [-] lift coefficient aircraft
CL = CL[0]
t_c =  get_t_c(lambda_2_rad,M_x, M_cruise,CL)                                   # [-] thickness over chord main wing
Cr_t= t_c*Cr                                                                    # [m] thinkness at the root chord
MAC = get_MAC(Cr, taper_ratio)                                                  # [m] mean aerodynamic chord main wing
y_MAC = get_y_MAC(b, Cr, MAC, Ct)                                               # [m] y-location of the MAC of the main wing
dihedral_rad = get_dihedral_rad(lambda_4_rad)                                   # [rad] dihedral angle of the main wing
lambda_le_rad = get_lambda_le_rad(lambda_4_rad, Cr, b, taper_ratio)             # [rad] leading edge sweep angle main wing

#canard parameters
A_c = 6  
b_c = [get_b(A_c,S_c[i]) for i in range(3)]
lambda_c_4_rad = get_lambda_4_rad(M_cruise,M_x)                                 # [rad] quarter chord sweep angle canard
taper_ratio_c = get_taper_ratio(lambda_c_4_rad)                                 # [-] taper ratio canard
lambda_c_2_rad = get_lambda_2_rad_canard(lambda_c_4_rad,A_c,taper_ratio_c)      # [rad] half chord sweep angle canard

Cr_c = [0] + [get_Cr_canard(S_c[i],taper_ratio_c,b_c[i]) for i in range(1,3)]   # [m] root chord length canard
Ct_c = [0] + [get_Ct_canard(Cr_c[i], taper_ratio_c) for i in range(1,3)]        # [m] tip chord length canard
CL_c = [0] + get_CL_canard(M_carried_canard_MTOW,rho,V_cruise,S_c)              # [-] lift coefficient aircraft
t_c_c =  [0]+get_t_c_canard(lambda_c_2_rad,M_x, M_cruise,CL_c)                  # [-] thickness over chord main wing
Cr_t_c= [t_c_c[i]*Cr_c[i] for i in range(3)]                                    # [m] thinkness at the root chord
MAC_c = [0]+get_MAC_canard(Cr_c, taper_ratio_c)                                 # [m] mean aerodynamic chord canard
y_MAC_c = [0]+get_y_MAC_canard(b_c, Cr_c, MAC_c, Ct_c)                          # [m] y-location of the MAC of the canard


#cg and masses of components
#NEED L_H FOR CLASS 2 (DAAAN EN STIJJN)
M_wing, M_eng, M_wing_group=get_mass_winggroup(MTOW)
M_fuselage, x_cg_fuselage=get_mass_fuselage(MTOW,l_f)
M_tail,x_cg_tail=get_mass_tail(MTOW,l_f)
M_fuselage_group, x_cg_fuselage_group=get_mass_fuselagegroup(M_fuselage,M_tail,x_cg_fuselage,x_cg_tail)
x_le_MAC=get_x_le_MAC(l_f,MAC,M_wing_group, M_fuselage_group)

x_cg_wing,x_cg_eng,x_cg_wing_group=get_cg_winggroup(x_le_MAC, MAC,M_wing, M_eng, M_wing_group )

l_h=[x_cg_tail[i]-x_cg_wing[i] for i in range(3)]




x_cg=get_x_cg(M_wing_group, M_fuselage_group,x_cg_wing_group,x_cg_fuselage_group) # [m] x-location of the centre of mass aircraft
y_cg = get_y_cg()                                                               # [m] y-location of the centre of mass aircraft
z_cg = get_z_cg(d_f_outer)                                                      # [m] z-location of the centre of mass aircraft

# Empennage parameters
V_h = 1.28                                                                      # [-] volume horizontal tail
A_h = 4.95                                                                      # [-] aspect ratio horizontal tail
taper_ratio_h = 0.39                                                            # [-] taper ratio horizontal tail
V_v = 0.1                                                                       # [-] volume vertical tail
A_v = 1.9                                                                       # [-] aspect ratio vertical tail
taper_ratio_v = 0.375                                                           # [-] taper ratio vertical tail
x_le_h = get_x_h(l_f)                                                           # [m] x-position leading edge horizontal tail
x_le_v = x_le_h                                                                 # [m] x-position leading edge vertical tail

S_h = get_S_h(S, MAC, x_cg, V_h, x_le_h)                                        # [m^2] surface area horizontal tail
S_v = get_S_v(S, b, x_cg, V_v, x_le_v)                                          # [m^2] surface area vertical tail
b_h = get_b_h(S_h, A_h)                                                         # [m] span horizontal tail
b_v = get_b_v(S_v, A_v)                                                         # [m] span vertical tail
Cr_h = get_Cr_h(S_h, taper_ratio_h, b_h)                                        # [m] root chord length horizontal tail
Ct_h = get_Ct_h(Cr_h, taper_ratio_h)                                            # [m] tip chord length horizontal tail
Cr_v = get_Cr_v(S_v, taper_ratio_v, b_v)                                        # [m] root chord lengh vertical tail
Ct_v = get_Ct_v(Cr_v, taper_ratio_v)                                            # [m] tip chord length vertical tail

lambda_h_le_rad = np.deg2rad(34)                                                # [rad] leading edge sweep angle horizontal tail
lambda_h_4_rad= get_lambda_4_rad_from_lambda_le(lambda_h_le_rad,Cr_h,b_h,taper_ratio_h)
lambda_h_2_rad=get_lambda_2_rad(lambda_h_4_rad,A_h,taper_ratio_h)


lambda_v_le_rad = np.deg2rad(40)                                                # [rad] leading edge sweep angle vertical tail
lambda_v_4_rad= get_lambda_4_rad_from_lambda_le(lambda_v_le_rad,Cr_v,b_v,taper_ratio_v)
lambda_v_2_rad=get_lambda_2_rad(lambda_v_4_rad,A_v,taper_ratio_v)

# engine specifics
x_engine = x_le_MAC                                                             # [m] x-location of the engine
y_engine = 0.3*b/2                                                              # [m] y-location of the engine
z_engine = 0                                                                    # [m] z-location of the engine
i_e_rad = np.deg2rad(i_e)                                                       # [rad] incidence angle of the engine

#undercarriage
tire_pressure = 430 * np.log(LCN) - 680                                         # [Pa] tire pressure mlg

weight_distribution = 0.10                                                      # [-] weight percentage on nose wheel
z_engine_clearance = z_engine - d_eng/2                                         # [m] z-location of lowest part of the engine

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
                  z_cg,z_mlg,l_n,l_m,y_engine,z_engine_clearance,d_eng)                         # [m] y-location of the mlg

x_nlg = get_x_nlg(x_cg,l_n)                                                     # [m] x-location of nlg
y_nlg = [0,0,0]                                                                 # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg


#MAKE THIS ITERABLE WITH THE LOADING DIAGRAM
# Airfoil Cl,max from javafoil for Re = [9*10^6, 17*10^6, 20*10^6]

Cl_max = [1.552, 1.582, 1.584]

# airfoil design
# CLmax: Wing CL max for three Re numbers: [9*10^6, 17*10^6, 20*10^6]
# CL_alpha: Wing CL_alpha for three configurations
Reto1, Re1, Reto3, CLdes, Cl_des, CL_alpha, CLmax, CLmaxto=airfoil(Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, lambda_le_rad, lambda_2_rad, b, taper_ratio, A, Cl_max)
CD0, CDcruise, LoverD, Wing, Fuselage, Nacelle, Tailplane=drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_f, d_f_outer, d_nacel, l_nacel, lambda_le_rad, CLdes)

"""Please add this one below"""
alpha_cruise_rad = np.deg2rad(0)                                                # [rad] angle of attack during cruise


# loadingdiagram=plot_loadingdiagram(Sland,CLmaxto,CLmax,CLmaxto,c,f,sigma, TOP, CD0,100,7100,100)

#create loading diagram with new Cl and Cd
#CD0_roskam, CD0_TO_roskam, CD0_land_roskam=dragcoefficient(Cfe,Swet_S)
#for i in range(3):
#    loadingdiagram=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,Vto1*kts_to_ms,c,f,sigma, TOP, CD0_roskam,100,7100,100)

#V_TO estimation to get Class-II drag calculation going

V_TO = [sqrt(2* x * 9.80665 /(rho_0 * S * Cl_TO)) for x in MTOW]



