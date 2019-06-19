# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt

import math
import modules.initialsizing_cg as initialcg
import modules.airfoil_calculations as airf
import modules.initialsizing_weights as initialw
import modules.initialsizing_planform as initialplanform
import modules.initialsizing_fuselage as initialfus
import modules.initialsizing_empennage as initialemp
import modules.initialsizing_undercarriage as initialunderc
#from modules.initialsizing_loading import *     # commented out because this import immediately runs the plot......
#from modules.performance.payload_range import *
import inputs.performance_inputs as inputperf
#from modules.performance.class2_performance_de5fs import get_thrust_required
import inputs.constants as const

 
#should move to constants
N_pax = [90,120,120]                                                            # [-] number of passengers
R = [4000E3,2000E3,4000E3]                                                      # [m] range of the aircraft
#inputs to this file 
T_W    = [0.29,0.29,0.29]                                                       # [-] thrust over weight ratio
#T_Wnew=[0.287,0.291,0.294]
W_S    = [4405, 4405 , 4405]   
#W_S -[4185,4185,4185]                                                 # [N/m^2] weight over wing surface area
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg] mass fuel fraction
#M_ffnew =[0.8671129569197487,0.9180333011040266,0.8707912564954988]
OEW = [34631.92,38223.31-360,38729.81]                                          # [kg] operational empty weight
#OEWnew=[36819.31992073485,40721.89132627823,40721.89132627823]


d_OEW1,d_OEW2=initialw.get_mass_efficiency(OEW)




#START SIZING 
# Fuselage parameters]
l_cutout=30/const.N_sa*const.seat_pitch + const.seat_pitch                                        #change with additional safety factors

l_cabin = initialfus.get_l_cabin(N_pax,const.N_sa)                                               # [m] cabin length UPDATE THIS TO THE REAL VALUE
l_cabin = [max(l_cabin)-l_cutout, max(l_cabin), max(l_cabin)]
d_f_inner = initialfus.get_d_f_inner(const.N_sa, const.seat_width, const.N_aisle,\
                          const.armrest, const.aisle_width, const.s_clearance)                    # [m] inner diameter fuselage
d_f_outer = initialfus.get_d_f_outer(d_f_inner)                                            # [m] outer diameter fuselage
l_nose = initialfus.get_l_nose(d_f_outer)                                                  # [m] nose length
l_tailcone = initialfus.get_l_tailcone(d_f_outer)                                          # [m] tailcone length
l_tail = initialfus.get_l_tail(d_f_outer)                                                  # [m] tail length
l_f = initialfus.get_l_fuselage(const.l_cockpit, l_cabin, l_tail)                                # [m] length fuselage
R_f = d_f_outer/2                                                               # [m] radius fuselage
S_fus=[R_f*2*math.pi*l_f[i] for i in range(3)]                              # [m^2] gross shell fuselage area

V_os = initialfus.get_overhead_volume(l_cabin)                                             # [m^3] overhead storage volume
V_cc = initialfus.get_cargo_volume(R_f,l_cabin)                                            # [m^3] total storage volume

Mtot_carry_on, Mtot_check_in, V_carry_on\
, V_check_in = initialfus.get_masses_volumes(N_pax, V_cc, V_os)                            # [kg,kg,m^3,m^3] mass and volume of check-in/carry-on luggage

V_cargo_available = initialfus.get_available_cargo_volume(V_cc,V_os,V_carry_on, V_check_in)# [m^3] available cargo volume

#CALCULATE MASSES BASED ON THE FUSALGE LAYOUT
M_pax_and_lugg=initialfus.get_passenger_luggage_mass(N_pax)
M_cargo_available=[V_cargo_available[i]*const.rho_cargo for i in range(3)]             # [kg] available cargo mass
M_payload=[M_cargo_available[i]+M_pax_and_lugg[i] for i in range(3)]
MTOW=initialw.get_TOW(OEW,M_payload,M_ff)
M_fuel = initialw.get_M_fuel(MTOW,M_ff)                                                  # [kg] fuel mass
T_req = initialw.get_T_req(T_W, MTOW)

#needed for class2 estimation
M_MZF    = [MTOW[i]-M_fuel[i] for i in range(3)]
M_carried_canard_MZF=[M_MZF[i]-M_MZF[0] for i in range(3)]              
M_carried_canard_MTOW=[MTOW[i]-MTOW[0] for i in range(3)]                                                    # [N] required thrust                          

"Change this when correct length of modular part is found, implement the l_cutout here" 
x_cargo = [[const.Xcargo1, const.Xcargo2], [const.Xcargo1+(2.), const.Xcargo2+(l_cutout)]\
           , [const.Xcargo1+(2.), const.Xcargo2+(l_cutout)]]     



#Propulsion
d_nacel = 1.1*const.d_fan                                                             # [m] diameter of engine nacelle
l_nacel = 1.1*const.l_eng                                                             # [m] length of the engine nacelle

# Wing parameters MOVE TO CONSTANTS
A = 9.5                                                                         # [-] aspect ration main wing
e = 0.85                                                                        # [-] 
S = initialplanform.get_S(MTOW,W_S)                                                             # [m^2] surface area main wing
#take canard into account
S_c = [S[0]-S[0],S[1]-S[0],S[2]-S[0]]

S = min(S)                                                                      # [m] update surface area to be the same for all config
b = initialplanform.get_b(A,S)                                                                  # [m] span main wing
lambda_4_rad = initialplanform.get_lambda_4_rad(const.M_cruise,const.M_x)                       # [rad] quarter chord sweep angle main wing
taper_ratio = initialplanform.get_taper_ratio(lambda_4_rad)                                     # [-] taper ratio main wing
lambda_2_rad = initialplanform.get_lambda_2_rad(lambda_4_rad,A,taper_ratio)                     # [rad] half chord sweep angle main wing
Cr = initialplanform.get_Cr(S,taper_ratio,b)                                                    # [m] root chord length main wing
Ct = initialplanform.get_Ct(Cr, taper_ratio)                                                    # [m] tip chord length main wing
CL = initialplanform.get_CL(MTOW,const.rho,const.V_cruise,S)                                                # [-] lift coefficient aircraft
CL = CL[0]
t_c =  initialplanform.get_t_c(lambda_2_rad,const.M_x, const.M_cruise,CL)                                   # [-] thickness over chord main wing
Cr_t= t_c*Cr                                                                    # [m] thinkness at the root chord
MAC = initialplanform.get_MAC(Cr, taper_ratio)                                                  # [m] mean aerodynamic chord main wing
y_MAC = initialplanform.get_y_MAC(b, Cr, MAC, Ct)                                               # [m] y-location of the MAC of the main wing
dihedral_rad = initialplanform.get_dihedral_rad(lambda_4_rad)                                   # [rad] dihedral angle of the main wing
lambda_le_rad = initialplanform.get_lambda_le_rad(lambda_4_rad, Cr, b, taper_ratio)             # [rad] leading edge sweep angle main wing

'DELETE THIS'
#canard parameters
A_c = 6
b_c = [initialplanform.get_b(A_c,S_c[i]) for i in range(3)]
lambda_c_4_rad = initialplanform.get_lambda_4_rad(const.M_cruise,const.M_x)                                 # [rad] quarter chord sweep angle canard
taper_ratio_c = initialplanform.get_taper_ratio(lambda_c_4_rad)                                 # [-] taper ratio canard
lambda_c_2_rad = initialplanform.get_lambda_2_rad_canard(lambda_c_4_rad,A_c,taper_ratio_c)      # [rad] half chord sweep angle canard

Cr_c = [0] + [initialplanform.get_Cr_canard(S_c[i],taper_ratio_c,b_c[i]) for i in range(1,3)]   # [m] root chord length canard
Ct_c = [0] + [initialplanform.get_Ct_canard(Cr_c[i], taper_ratio_c) for i in range(1,3)]        # [m] tip chord length canard
CL_c = [0] + initialplanform.get_CL_canard(M_carried_canard_MTOW,const.rho,const.V_cruise,S_c)              # [-] lift coefficient aircraft
t_c_c =  [0]+initialplanform.get_t_c_canard(lambda_c_2_rad,const.M_x, const.M_cruise,CL_c)                  # [-] thickness over chord main wing
Cr_t_c= [t_c_c[i]*Cr_c[i] for i in range(3)]                                    # [m] thinkness at the root chord
MAC_c = [0]+initialplanform.get_MAC_canard(Cr_c, taper_ratio_c)                                 # [m] mean aerodynamic chord canard
y_MAC_c = [0]+initialplanform.get_y_MAC_canard(b_c, Cr_c, MAC_c, Ct_c)                          # [m] y-location of the MAC of the canard


#cg and masses of components Delete this later on 
#NEED L_H FOR CLASS 2 (DAAAN EN STIJJN)
M_wing, M_eng, M_wing_group=initialcg.get_mass_winggroup(MTOW)
M_fuselage, x_cg_fuselage=initialcg.get_mass_fuselage(MTOW,l_f)
M_tail,x_cg_tail=initialcg.get_mass_tail(MTOW,l_f)
M_fuselage_group, x_cg_fuselage_group=initialcg.get_mass_fuselagegroup(M_fuselage,M_tail,x_cg_fuselage,x_cg_tail)
x_le_MAC=initialcg.get_x_le_MAC(l_f,MAC,M_wing_group, M_fuselage_group)
#x_le_MAC=[x_le_MAC[i]+4 for i in range(3)]
x_le_w = initialplanform.get_le_wing(y_MAC,x_le_MAC, lambda_2_rad, MAC, Cr)

x_cg_wing,x_cg_eng,x_cg_wing_group=initialcg.get_cg_winggroup(x_le_MAC, MAC,M_wing, M_eng, M_wing_group )

l_h=[x_cg_tail[i]-x_cg_wing[i] for i in range(3)]




x_cg=initialcg.get_x_cg(M_wing_group, M_fuselage_group,x_cg_wing_group,x_cg_fuselage_group) # [m] x-location of the centre of mass aircraft
y_cg = initialcg.get_y_cg()                                                               # [m] y-location of the centre of mass aircraft
z_cg = initialcg.get_z_cg(d_f_outer)                                                      # [m] z-location of the centre of mass aircraft

# Empennage parameters
V_h = 1.28                                                                      # [-] volume horizontal tail
A_h = 4.95                                                                      # [-] aspect ratio horizontal tail
taper_ratio_h = 0.39                                                            # [-] taper ratio horizontal tail
V_v = 0.1                                                                       # [-] volume vertical tail
A_v = 1.9                                                                       # [-] aspect ratio vertical tail
taper_ratio_v = 0.375                                                           # [-] taper ratio vertical tail
x_le_h = initialemp.get_x_h(l_f)                                                           # [m] x-position leading edge horizontal tail
x_le_v = x_le_h                                                                 # [m] x-position leading edge vertical tail

S_h = initialemp.get_S_h(S, MAC, x_cg, V_h, x_le_h)                                        # [m^2] surface area horizontal tail
S_v = initialemp.get_S_v(S, b, x_cg, V_v, x_le_v)                                          # [m^2] surface area vertical tail
b_h = initialemp.get_b_h(S_h, A_h)                                                         # [m] span horizontal tail
b_v = initialemp.get_b_v(S_v, A_v)                                                         # [m] span vertical tail
Cr_h = initialemp.get_Cr_h(S_h, taper_ratio_h, b_h)                                        # [m] root chord length horizontal tail
Ct_h = initialemp.get_Ct_h(Cr_h, taper_ratio_h)                                            # [m] tip chord length horizontal tail
Cr_v = initialemp.get_Cr_v(S_v, taper_ratio_v, b_v)                                        # [m] root chord lengh vertical tail
Ct_v = initialemp.get_Ct_v(Cr_v, taper_ratio_v)                                            # [m] tip chord length vertical tail

lambda_h_le_rad = np.deg2rad(34)                                                # [rad] leading edge sweep angle horizontal tail
lambda_h_4_rad= initialplanform.get_lambda_4_rad_from_lambda_le(lambda_h_le_rad,Cr_h,b_h,taper_ratio_h) # [rad] quarter chord sweep angle
lambda_h_2_rad= initialplanform.get_lambda_2_rad(lambda_h_4_rad,A_h,taper_ratio_h)               # [rad] half chord sweep angle


lambda_v_le_rad = np.deg2rad(40)                                                # [rad] leading edge sweep angle vertical tail
lambda_v_4_rad= initialplanform.get_lambda_4_rad_from_lambda_le(lambda_v_le_rad,Cr_v,b_v,taper_ratio_v) # [rad] quarter chord sweep angle
lambda_v_2_rad=initialplanform.get_lambda_2_rad(lambda_v_4_rad,A_v,taper_ratio_v)               # [rad] half chord sweep angle

# engine specifics
x_engine = x_le_MAC                                                             # [m] x-location of the engine
y_engine = 0.3*b/2                                                              # [m] y-location of the engine
z_engine = 0                                                                    # [m] z-location of the engine
i_e_rad = np.deg2rad(const.i_e)                                                 # [rad] incidence angle of the engine

#undercarriage
tire_pressure = 430 * np.log(const.LCN) - 680                                   # [Pa] tire pressure mlg

weight_distribution = 0.16                                                      # [-] weight percentage on nose wheel
z_engine_clearance = z_engine - const.d_eng/2                                   # [m] z-location of lowest part of the engine


P_mw = initialunderc.get_P_mw(MTOW,const.N_mw,weight_distribution)                                  # [N] static loading on mw
P_nw = initialunderc.get_P_nw(MTOW,const.N_nw,weight_distribution)                                  # [N] static loading on nw

x_mlg = initialunderc.get_x_mlg(z_cg,const.theta_rad,const.beta_rad, x_cg, const.stroke,l_f)  # [m] x-location of the mlg
x_mlg[1]=min(x_mlg)+l_cutout
x_mlg[2]=min(x_mlg)+l_cutout
z_mlg = initialunderc.get_z_mlg(x_mlg,const.beta_rad,x_cg, z_cg)                              # [m] z-location of the mlg
z_mlg=max(z_mlg)


l_m = initialunderc.get_l_mw(x_mlg,x_cg)                                                      # [m] mlg distance from c.g
l_n = initialunderc.get_l_nw(l_m,P_mw,const.N_mw,P_nw,const.N_nw)                             # [m] nlg distance from c.g

x_nlg = initialunderc.get_x_nlg(x_cg,l_n)                                                     # [m] x-location of nlg
x_nlg=min(x_nlg)
l_n=[x_cg[i]-x_nlg for i in range(3)]


new_beta_rad=initialunderc.get_new_beta_config23(z_mlg,z_cg,l_m[1])
new_beta_deg=new_beta_rad*180/math.pi
#new_theta_rad=check_new_scrape_angle(l_f,x_mlg,d_f_outer,z_mlg)
#new_theta_deg=[new_theta_rad[i]*180/pi for i in range(3)]

y_mlg = initialunderc.get_y_mlg(b,dihedral_rad,const.psi_rad,const.phi_rad,\
               z_cg,z_mlg,l_n,l_m,y_engine,z_engine_clearance,const.d_eng)            # [m] y-location of the mlg
#


y_nlg = [0,0,0]                                                                 # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg

L_strut_mlg= -z_mlg-const.D_mlg/2
L_strut_nlg= -z_nlg-const.D_nlg/2
D_strut_mlg= initialunderc.get_d_lg(max(P_mw)*2,L_strut_mlg)
D_strut_nlg= initialunderc.get_d_lg(max(P_nw)*2,L_strut_nlg)
#MAKE THIS ITERABLE WITH THE LOADING DIAGRAM
# Airfoil Cl,max from javafoil for Re = [9*10^6, 17*10^6, 20*10^6]

Cl_max = [1.552, 1.582, 1.584]

# airfoil design
# CLmax: Wing CL max for three Re numbers: [9*10^6, 17*10^6, 20*10^6]
# CL_alpha: Wing CL_alpha for three configurations
Reto1, Re1, Reto3, CLdes, Cl_des, CL_alpha, CLmax, CLmaxto=airf.airfoil(Ct, Cr, MTOW, inputperf.FF1, inputperf.FF2, inputperf.FF3, inputperf.FF4, inputperf.FF5, S, lambda_le_rad, lambda_2_rad, b, taper_ratio, A, Cl_max)
CD0, CDcruise, LoverD, Wing, Fuselage, Nacelle, Tailplane=airf.drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_f, d_f_outer, d_nacel, l_nacel, lambda_le_rad, CLdes)


"""Please add this one below"""
alpha_cruise_rad = np.deg2rad(0)                                                # [rad] angle of attack during cruise


# loadingdiagram=plot_loadingdiagram(Sland,CLmaxto,CLmax,CLmaxto,c,f,sigma, TOP, CD0,100,7100,100)


#create loading diagram with new Cl and Cd
#CD0_roskam, CD0_TO_roskam, CD0_land_roskam=dragcoefficient(Cfe,Swet_S)
#for i in range(3):
#    loadingdiagram=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,Vto1*kts_to_ms,c,f,sigma, TOP, CD0_roskam,100,7100,100)
# loadingdiagram=plot_loadingdiagram(Sland,CLmaxto,CLmax,CLmaxto,c,f,sigma, TOP, CD0,100,7100,100)

#V_TO estimation to get Class-II drag calculation going

V_TO = [math.sqrt(2* x * 9.80665 /(const.rho_0 * S * inputperf.Cl_TO)) for x in MTOW]

#thrust_cruise = [get_thrust_required(isa(H_m)[2], V_cruise, S, CDcruise[i]) for i in range(3)]