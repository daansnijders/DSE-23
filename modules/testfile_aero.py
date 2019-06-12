# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:32:22 2019

@author: Sybren
"""

from inputs.constants import *
from inputs.performance_inputs import *
from inputs.concept_1 import *
from modules.Aerodynamics import *

""" Yet unknown values needed to start the simulation """
D_nlg = 0.5
b_nlg = 0.25
D_mlg = 1.3
b_mlg = 0.35
D_strutt_nlg = 0.15
D_strutt_mlg = 0.2
CL_alpha_h = 0
de_da_h = 0
i_h = 0
alpha0L_h = 0
CL_alpha_c = 0
de_da_c = 0
i_c = 0
alpha0L_c = 0
l_fueltank = 1.5
d_fueltank = 0.3
S_ef = 0.75 * S_h[0] 
i_n = 0
delta_C_L_h = 0.3
delta_C_L_c = 0.25
Delta_C_L_flap = Cl_land - Cl_clean

#""" HLD design """
#config1_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
#SWF, b_flap, SWF_LE, b_slat = config1_HLD.HLD()
##print(SWF, b_flap, SWF_LE, b_slat)
#
#""" Drag classII estimations """
##Configuration 1
#config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[0],CL_alpha,l_cockpit, l_cabin[0], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, b_nlg, D_strutt_nlg, x_mlg[0], z_mlg, D_mlg, b_mlg, D_strutt_mlg, lambda_h_2_rad[0], lambda_v_2_rad[0], MAC_c[0], Cr_v[0], Ct_v[0], Cr_h[0], Ct_h[0], S_h[0], S_v[0], S_c[0], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#
#CD_w_sub1, CD_w_trans1, CD_0_w1 = config1_Drag.wing_drag()
#CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag()
#CD_h_sub1, CD_v_sub1, CD_c_sub1, CD_h_trans1, CD_v_trans1, CD_c_trans1 = config1_Drag.empennage_drag()
#CD_nacel_sub1, CD_nacel_trans1 = config1_Drag.nacelle_drag()
#CD_flap_TO1, CD_flap_land1, CD_slat1 = config1_Drag.flaps_drag(CD_0_w1)
#CD_gear1 = config1_Drag.landinggear_drag()
#CD_ws1 = config1_Drag.windshield_drag()
#CD_store_sub1, C_D_store_trans1 = config1_Drag.store_drag()
#CD_trim1 = config1_Drag.trim_drag()
#CD_spoiler1 = config1_Drag.spoiler_drag()
#
#CD_cruise1 = CD_w_trans1 + CD_fus_trans1 + CD_h_trans1 + CD_v_trans1 + CD_c_trans1 + CD_nacel_trans1 + CD_ws1 + C_D_store_trans1 + CD_trim1
#CD_TO1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_TO1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
#CD_land1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_land1 + CD_slat1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
#
#
##Configuration 2
#config2_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[1],CL_alpha,l_cockpit, l_cabin[1], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, b_nlg, D_strutt_nlg, x_mlg[1], z_mlg, D_mlg, b_mlg, D_strutt_mlg, lambda_h_2_rad[1], lambda_v_2_rad[1], MAC_c[1], Cr_v[1], Ct_v[1], Cr_h[1], Ct_h[1], S_h[1], S_v[1], S_c[1], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#
#CD_w_sub2, CD_w_trans2, CD_0_w2 = config2_Drag.wing_drag()
#CD_fus_sub2, CD_fus_trans2 = config2_Drag.fuse_drag()
#CD_h_sub2, CD_v_sub2, CD_c_sub2, CD_h_trans2, CD_v_trans2, CD_c_trans2 = config2_Drag.empennage_drag()
#CD_nacel_sub2, CD_nacel_trans2 = config2_Drag.nacelle_drag()
#CD_flap_TO2, CD_flap_land2, CD_slat2 = config2_Drag.flaps_drag(CD_0_w2)
#CD_gear2 = config2_Drag.landinggear_drag()
#CD_ws2 = config2_Drag.windshield_drag()
#CD_store_sub2, C_D_store_trans2 = config2_Drag.store_drag()
#CD_trim2 = config2_Drag.trim_drag()
#CD_spoiler2 = config2_Drag.spoiler_drag()
#
#CD_cruise2 = CD_w_trans2 + CD_fus_trans2 + CD_h_trans2 + CD_v_trans2 + CD_c_trans2 + CD_nacel_trans2 + CD_ws2 + C_D_store_trans2 + CD_trim2
#CD_TO2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_TO2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2
#CD_land2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_land2 + CD_slat2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2
#
##Configuration 3
#config3_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[2],CL_alpha,l_cockpit, l_cabin[2], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, b_nlg, D_strutt_nlg, x_mlg[2], z_mlg, D_mlg, b_mlg, D_strutt_mlg, lambda_h_2_rad[2], lambda_v_2_rad[2], MAC_c[2], Cr_v[2], Ct_v[2], Cr_h[2], Ct_h[2], S_h[2], S_v[2], S_c[2], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#
#CD_w_sub3, CD_w_trans3, CD_0_w3 = config3_Drag.wing_drag()
#CD_fus_sub3, CD_fus_trans3 = config3_Drag.fuse_drag()
#CD_h_sub3, CD_v_sub3, CD_c_sub3, CD_h_trans3, CD_v_trans3, CD_c_trans3 = config3_Drag.empennage_drag()
#CD_nacel_sub3, CD_nacel_trans3 = config3_Drag.nacelle_drag()
#CD_flap_TO3, CD_flap_land3, CD_slat3 = config3_Drag.flaps_drag(CD_0_w3)
#CD_gear3 = config3_Drag.landinggear_drag()
#CD_ws3 = config3_Drag.windshield_drag()
#CD_store_sub3, C_D_store_trans3 = config3_Drag.store_drag()
#CD_trim3 = config3_Drag.trim_drag()
#CD_spoiler3 = config3_Drag.spoiler_drag()
#
#CD_cruise3 = CD_w_trans3 + CD_fus_trans3 + CD_h_trans3 + CD_v_trans3 + CD_c_trans3 + CD_nacel_trans3 + CD_ws3 + C_D_store_trans3 + CD_trim3
#CD_TO3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_TO3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3
#CD_land3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_land3 + CD_slat3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3
#
#print(CD_cruise1, CD_cruise2, CD_cruise3)
#print(CD_TO1, CD_TO2, CD_TO3)
#print(CD_land1, CD_land2, CD_land3)

""" Lift classII estimations """
alpha_0_l = -5.4
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
C_l_alpha_M75 = C_l_alpha / sqrt(1 - M_cruise**2)
i_w= 0
wing_twist = 0


#Configuration 1
config1_Lift = Lift(S,A,rho,rho_0,l_f[0],V_cruise,M_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c, i_h, x_le_MAC)
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1 = config1_Lift.Airfoil_lift_flaps()
CL_alpha1, alpha_0_L1, CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)


#print(delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1)
print(CL_alpha1, alpha_0_L1, CL_max1)

#Configuration 2
config2_Lift = Lift(S,A,rho,rho_0,l_f[1],V_cruise,M_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c, i_h, x_le_MAC)
delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2 = config2_Lift.Airfoil_lift_flaps()
CL_alpha2, alpha_0_L2, CL_max2 = config2_Lift.Airplane_lift(CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2)

#print(delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2)
print(CL_alpha2, alpha_0_L2, CL_max2)

#Configuration 3
config3_Lift = Lift(S,A,rho,rho_0,l_f[2],V_cruise,M_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c, i_h, x_le_MAC)
delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3 = config3_Lift.Airfoil_lift_flaps()
CL_alpha3, alpha_0_L3, CL_max3 = config3_Lift.Airplane_lift(CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3)

#print(delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3)
print(CL_alpha3, alpha_0_L3, CL_max3)
