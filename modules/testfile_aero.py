# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:32:22 2019

@author: Sybren
"""

import numpy as np
import math
import inputs.constants as const
import inputs.performance_inputs as perf
import inputs.concept_1 as conc1
import modules.Aerodynamics as aero

""" Yet unknown values needed to start the simulation """
D_nlg = 0.5
b_nlg = 0.25
D_mlg = 1.3
b_mlg = 0.35
D_strutt_nlg = 0.15
D_strutt_mlg = 0.2
CL_alpha_h = 0
i_h = 0
alpha0L_h = 0
CL_alpha_c = 0
de_da_c = 0
de_da = 0
i_c = 0
alpha0L_c = 0
l_fueltank = 1.5
d_fueltank = 0.3
S_elev = 0.75 * conc1.S_h[0] 
i_n = 0
delta_C_L_h = 0.3
delta_C_L_c = 0.25
Delta_C_L_flap = 2.4 - 1.5          #Cl_land - Cl_clean from performanca_input

""" HLD design """
config1_HLD = aero.HLD_class(perf.Cl_land,perf.Cl_clean,conc1.S,conc1.A,conc1.lambda_4_rad,conc1.taper_ratio,conc1.CL_alpha,conc1.lambda_le_rad,conc1.Cr,conc1.d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config1_HLD.HLD()
#print(SWF, b_flap, SWF_LE, b_slat)

""" Drag classII estimations """
#Configuration 1
config1_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[0], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[0], conc1.lambda_v_2_rad[0], conc1.MAC_c[0], conc1.Cr_v[0], conc1.Ct_v[0], conc1.Cr_h[0], conc1.Ct_h[0], conc1.S_h[0], conc1.S_v[0], conc1.S_c[0], CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

CD_w_sub1, CD_w_trans1, CD0_w1 = config1_Drag.wing_drag()
CD0_fus1, CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag()
CD0_h_tail1, CD0_v_tail1, CD0_c_tail1, CD_h_sub1, CD_v_sub1, CD_c_sub1, CD_h_trans1, CD_v_trans1, CD_c_trans1 = config1_Drag.empennage_drag()
CD0_nacel1, CD_nacel_sub1, CD_nacel_trans1 = config1_Drag.nacelle_drag()
CD_flap_TO1, CD_flap_land1, CD_slat1 = config1_Drag.flaps_drag(CD0_w1)
CD_gear1 = config1_Drag.landinggear_drag()
CD_ws1 = config1_Drag.windshield_drag()
CD0_store1, CD_store_sub1, C_D_store_trans1 = config1_Drag.store_drag()
CD_trim1 = config1_Drag.trim_drag()
CD_spoiler1 = config1_Drag.spoiler_drag()

CD_cruise1 = CD_w_trans1 + CD_fus_trans1 + CD_h_trans1 + CD_v_trans1 + CD_c_trans1 + CD_nacel_trans1 + CD_ws1 + C_D_store_trans1 + CD_trim1
CD_TO1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_TO1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
CD_land1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_land1 + CD_slat1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
CD0_1 = CD0_w1 + CD0_fus1 + CD0_h_tail1 + CD0_v_tail1 + CD0_c_tail1 + CD0_nacel1 + CD0_store1

#Configuration 2
config2_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[1],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[1], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[1], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[1], conc1.lambda_v_2_rad[1], conc1.MAC_c[1], conc1.Cr_v[1], conc1.Ct_v[1], conc1.Cr_h[1], conc1.Ct_h[1], conc1.S_h[1], conc1.S_v[1], conc1.S_c[1], CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

CD_w_sub2, CD_w_trans2, CD0_w2 = config2_Drag.wing_drag()
CD0_fus2, CD_fus_sub2, CD_fus_trans2 = config2_Drag.fuse_drag()
CD0_h_tail2, CD0_v_tail2, CD0_c_tail2, CD_h_sub2, CD_v_sub2, CD_c_sub2, CD_h_trans2, CD_v_trans2, CD_c_trans2 = config2_Drag.empennage_drag()
CD0_nacel2, CD_nacel_sub2, CD_nacel_trans2 = config2_Drag.nacelle_drag()
CD_flap_TO2, CD_flap_land2, CD_slat2 = config2_Drag.flaps_drag(CD0_w2)
CD_gear2 = config2_Drag.landinggear_drag()
CD_ws2 = config2_Drag.windshield_drag()
CD0_store2, CD_store_sub2, C_D_store_trans2 = config2_Drag.store_drag()
CD_trim2 = config2_Drag.trim_drag()
CD_spoiler2 = config2_Drag.spoiler_drag()

CD_cruise2 = CD_w_trans2 + CD_fus_trans2 + CD_h_trans2 + CD_v_trans2 + CD_c_trans2 + CD_nacel_trans2 + CD_ws2 + C_D_store_trans2 + CD_trim2
CD_TO2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_TO2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2
CD_land2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_land2 + CD_slat2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2
CD0_2 = CD0_w2 + CD0_fus2 + CD0_h_tail2 + CD0_v_tail2 + CD0_c_tail2 + CD0_nacel2 + CD0_store2

#Configuration 3
config3_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[2],const.V_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[2],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[2], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[2], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[2], conc1.lambda_v_2_rad[2], conc1.MAC_c[2], conc1.Cr_v[2], conc1.Ct_v[2], conc1.Cr_h[2], conc1.Ct_h[2], conc1.S_h[2], conc1.S_v[2], conc1.S_c[2], CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

CD_w_sub3, CD_w_trans3, CD0_w3 = config3_Drag.wing_drag()
CD0_fus3, CD_fus_sub3, CD_fus_trans3 = config3_Drag.fuse_drag()
CD0_h_tail3, CD0_v_tail3, CD0_c_tail3, CD_h_sub3, CD_v_sub3, CD_c_sub3, CD_h_trans3, CD_v_trans3, CD_c_trans3 = config3_Drag.empennage_drag()
CD0_nacel3, CD_nacel_sub3, CD_nacel_trans3 = config3_Drag.nacelle_drag()
CD_flap_TO3, CD_flap_land3, CD_slat3 = config3_Drag.flaps_drag(CD0_w3)
CD_gear3 = config3_Drag.landinggear_drag()
CD_ws3 = config3_Drag.windshield_drag()
CD0_store3, CD_store_sub3, C_D_store_trans3 = config3_Drag.store_drag()
CD_trim3 = config3_Drag.trim_drag()
CD_spoiler3 = config3_Drag.spoiler_drag()

CD_cruise3 = CD_w_trans3 + CD_fus_trans3 + CD_h_trans3 + CD_v_trans3 + CD_c_trans3 + CD_nacel_trans3 + CD_ws3 + C_D_store_trans3 + CD_trim3
CD_TO3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_TO3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3
CD_land3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_land3 + CD_slat3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3
CD0_3 = CD0_w3 + CD0_fus3 + CD0_h_tail3 + CD0_v_tail3 + CD0_c_tail3 + CD0_nacel3 + CD0_store3

#print(CD_cruise1, CD_cruise2, CD_cruise3)
#print(CD_TO1, CD_TO2, CD_TO3)
#print(CD_land1, CD_land2, CD_land3)

""" Lift classII estimations """
alpha_0_l = -5.4
alpha_star_l = 10        # from -7 to 3 deg
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
C_l_alpha_M75 = C_l_alpha / math.sqrt(1 - const.M_cruise**2)
i_w= 0
wing_twist = -3 

#Configuration 1
config1_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[0], conc1.lambda_c_2_rad, i_c, conc1.S_h[0], conc1.S_c[0], i_h, conc1.x_le_MAC[0], b_flap, SWF)
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1 = config1_Lift.Airfoil_lift_flaps()
C_L_w1, CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1 = config1_Lift.Wing_lift()
delta_CL_w1, delta_CL_alpha_w1, delta_CL_max_w1 = config1_Lift.Wing_lift_flaps(delta_cl_flap1,CL_alpha_w1,C_l_alpha,(delta_clmax_flap1 + delta_clmax_krueger1),b_slat)
CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1, de_da1, de_da_c1, alpha_CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)
delta_CL1, delta_CL_alpha1, delta_CL_max1 = config1_Lift.Airplane_lift_flaps(delta_CL_w1, CL_alpha_h1, CL_alpha_c1, delta_CL_alpha_w1, de_da1, delta_CL_max_w1)

#print(delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1)
#print(CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1)
#print(delta_CL1, delta_CL_alpha1, delta_CL_max1)

#Configuration 2
config2_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[1], conc1.lambda_c_2_rad, i_c, conc1.S_h[1], conc1.S_c[1], i_h, conc1.x_le_MAC[1], b_flap, SWF)
delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2 = config2_Lift.Airfoil_lift_flaps()
C_L_w2, CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2 = config2_Lift.Wing_lift()
delta_CL_w2, delta_CL_alpha_w2, delta_CL_max_w2 = config2_Lift.Wing_lift_flaps(delta_cl_flap2,CL_alpha_w2,C_l_alpha,(delta_clmax_flap2 + delta_clmax_krueger2),b_slat)
CL_alpha_h2, CL_alpha_c2, CL_alpha2, alpha_0_L2, CL_max2, de_da2, de_da_c2, alpha_CL_max2 = config2_Lift.Airplane_lift(CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2)
delta_CL2, delta_CL_alpha2, delta_CL_max2 = config2_Lift.Airplane_lift_flaps(delta_CL_w2, CL_alpha_h2, CL_alpha_c2, delta_CL_alpha_w2, de_da2, delta_CL_max_w2)

#print(delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2)
#print(CL_alpha_h2, CL_alpha_c2, CL_alpha2, alpha_0_L2, CL_max2)

#Configuration 3
config3_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[2],const.V_cruise,const.M_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[2], conc1.lambda_c_2_rad, i_c, conc1.S_h[2], conc1.S_c[2], i_h, conc1.x_le_MAC[2], b_flap, SWF)
delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3 = config3_Lift.Airfoil_lift_flaps()
C_L_w3, CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3 = config3_Lift.Wing_lift()
delta_CL_w3, delta_CL_alpha_w3, delta_CL_max_w3 = config3_Lift.Wing_lift_flaps(delta_cl_flap3,CL_alpha_w3,C_l_alpha,(delta_clmax_flap3 + delta_clmax_krueger3),b_slat)
CL_alpha_h3, CL_alpha_c3, CL_alpha3, alpha_0_L3, CL_max3, de_da3, de_da_c3, alpha_CL_max3 = config3_Lift.Airplane_lift(CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3)
delta_CL3, delta_CL_alpha3, delta_CL_max3 = config3_Lift.Airplane_lift_flaps(delta_CL_w3, CL_alpha_h3, CL_alpha_c3, delta_CL_alpha_w3, de_da3, delta_CL_max_w3)

#print(delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3)
#print(CL_alpha_h3, CL_alpha_c3, CL_alpha3, alpha_0_L3, CL_max3)


plot = config1_Lift.CL_alpha_plot(CL_alpha1, alpha_0_L1, CL_max1, alpha_CL_max1, delta_CL1, delta_CL_alpha1, delta_CL_max1)


""" Moment classII estimations """
'Constants'
x_ref = 0.5                     #Own decision based on nothing
cl_des_airfoil = 0.651          #CL at zero angle of attack from JAVAfoil (RE = 17*10^6, 65-615)

#Configuration 1
config1_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[0], conc1.l_tail)
cm_des_airfoil1, dcm_dcl_airfoil1 = config1_Moment.Airfoil_moment()
delta_cm_flap1, delta_cm_krueger1 = config1_Moment.Airfoil_moment_flaps(cm_des_airfoil1)
Cm0_w_sub1, Cm0_w_trans1, dCm_dCl_w1 = config1_Moment.Wing_moment()
delta_Cm_w_flaps1, delta_Cm_w_krueger1 = config1_Moment.Wing_moment_flaps(Cm0_w_sub1)
#print(cm_des_airfoil1, dcm_dcl_airfoil1)
#print(delta_cm_flap1, delta_cm_krueger1)
#print(Cm0_w_sub1, Cm0_w_trans1, dCm_dCl_w1)
#print(delta_Cm_w_flaps1, delta_Cm_w_krueger1)


#Configuration 2
config2_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[1], conc1.l_tail)
cm_des_airfoil2, dcm_dcl_airfoil2 = config2_Moment.Airfoil_moment()
delta_cm_flap2, delta_cm_krueger2 = config2_Moment.Airfoil_moment_flaps(cm_des_airfoil2)
Cm0_w_sub2, Cm0_w_trans2, dCm_dCl_w2 = config2_Moment.Wing_moment()
delta_Cm_w_flaps2, delta_Cm_w_krueger2 = config2_Moment.Wing_moment_flaps(Cm0_w_sub2)

#Configuration 3
config3_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[1], conc1.l_tail)
cm_des_airfoil3, dcm_dcl_airfoil3 = config3_Moment.Airfoil_moment()
delta_cm_flap3, delta_cm_krueger3 = config3_Moment.Airfoil_moment_flaps(cm_des_airfoil3)
Cm0_w_sub3, Cm0_w_trans3, dCm_dCl_w3 = config3_Moment.Wing_moment()
delta_Cm_w_flaps3, delta_Cm_w_krueger3 = config3_Moment.Wing_moment_flaps(Cm0_w_sub3)

#""" Test """
#config1_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[0], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[0], conc1.lambda_v_2_rad[0], conc1.MAC_c[0], conc1.Cr_v[0], conc1.Ct_v[0], conc1.Cr_h[0], conc1.Ct_h[0], conc1.S_h[0], conc1.S_v[0], conc1.S_c[0], CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#config4_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[0], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[0], conc1.lambda_v_2_rad[0], conc1.MAC_c[0], conc1.Cr_v[0], conc1.Ct_v[0], conc1.Cr_h[0], conc1.Ct_h[0], conc1.S_h[0], conc1.S_v[0], conc1.S_c[0], CL_alpha_h1, de_da1, i_h, alpha_0_L1, conc1.A_h, CL_alpha_c1, de_da_c1, i_c, alpha_0_L1, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#CD_w_sub4, CD_w_trans4, CD_0_w4 = config4_Drag.wing_drag()
#CD_fus_sub4, CD_fus_trans4 = config4_Drag.fuse_drag()
#CD_h_sub4, CD_v_sub4, CD_c_sub4, CD_h_trans4, CD_v_trans4, CD_c_trans4 = config4_Drag.empennage_drag()
#CD_nacel_sub4, CD_nacel_trans4 = config4_Drag.nacelle_drag()
#CD_flap_TO4, CD_flap_land4, CD_slat4 = config4_Drag.flaps_drag(CD_0_w1)
#CD_gear4 = config4_Drag.landinggear_drag()
#CD_ws4 = config4_Drag.windshield_drag()
#CD_store_sub4, C_D_store_trans4 = config4_Drag.store_drag()
#CD_trim4 = config4_Drag.trim_drag()
#CD_spoiler4 = config4_Drag.spoiler_drag()
#
#CD_cruise4 = CD_w_trans4 + CD_fus_trans4 + CD_h_trans4 + CD_v_trans4 + CD_c_trans4 + CD_nacel_trans4 + CD_ws4 + C_D_store_trans4 + CD_trim4
#CD_TO4 = CD_w_sub4 + CD_fus_sub4 + CD_h_sub4 + CD_v_sub4 + CD_c_sub4 + CD_flap_TO4 + CD_gear4 + CD_ws4 + CD_store_sub4 + CD_trim4
#CD_land4 = CD_w_sub4 + CD_fus_sub4 + CD_h_sub4 + CD_v_sub4 + CD_c_sub4 + CD_flap_land4 + CD_slat4 + CD_gear4 + CD_ws4 + CD_store_sub4 + CD_trim4


C_L, C_L_flaps = config1_Lift.get_CL(CL_alpha1, alpha_0_L1, CL_max1, alpha_CL_max1, delta_CL1, delta_CL_alpha1, delta_CL_max1, 0.58)
#print (C_L, C_L_flaps)