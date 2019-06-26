# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 09:08:24 2019

@author: Sybren
"""

""" Sensitivity analysis Aerodynamics """

import numpy as np
import math
import inputs.constants as const
import inputs.performance_inputs as perf
import inputs.concept_1 as conc1
import modules.sensitivity_aero_2 as aero
#import Output.connection_departments as connect

""" Input parameters """
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

""" Lift """
alpha_0_l = -5.4
alpha_star_l = 10        # from -7 to 3 deg
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
C_l_alpha_M75 = C_l_alpha / math.sqrt(1 - const.M_cruise**2)
i_w= 0
wing_twist = -3         #DEG

Kcw = [0.75,0.85,0.95,1.05]
delta_CL = []
delta_CL_max = []

config1_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[0], conc1.lambda_c_2_rad, i_c, conc1.S_h[0], conc1.S_c[0], i_h, conc1.x_le_MAC[0], b_flap, SWF)
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1, delta_cl_flap_TO1, delta_clmax_TO1, clalpha_TO1 = config1_Lift.Airfoil_lift_flaps()
C_L_w1, CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1 = config1_Lift.Wing_lift()
delta_CL_w1, delta_CL_alpha_w1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1 = config1_Lift.Wing_lift_flaps(delta_cl_flap1,CL_alpha_w1,C_l_alpha,(delta_clmax_flap1 + delta_clmax_krueger1),b_slat, delta_cl_flap_TO1, delta_clmax_TO1)
CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1, de_da1, de_da_c1, alpha_CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)

for i in range(len(Kcw)):
    delta_CL1, delta_CL_alpha1, delta_CL_max1, delta_CL_TO1, delta_CL_alpha_TO1, delta_CL_max_TO1 = config1_Lift.Airplane_lift_flaps(delta_CL_w1, CL_alpha_h1, CL_alpha_c1, delta_CL_alpha_w1, de_da1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1, Kcw[i])
    delta_CL.append(delta_CL1)
    delta_CL_max.append(delta_CL_max1)

delta_cl = (np.array(delta_CL) - delta_CL[2])/delta_CL[2] * 100 
delta_clmax = (np.array(delta_CL_max) - delta_CL_max[2])/delta_CL_max[2] * 100

""" Drag """

CD_wave_w = [0.001,0.0015,0.002,0.003,0.004]
CD_wave_fus = [0.0025,0.75*0.005,0.005,0.0075,0.01]
CD_wave_ch = [0.001,0.0015,0.002,0.003,0.004]
CD_wave_nac = [0.001,0.0015,0.002,0.003,0.004]

config1_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[0], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[0], conc1.lambda_v_2_rad[0], conc1.MAC_c[0], conc1.Cr_v[0], conc1.Ct_v[0], conc1.Cr_h[0], conc1.Ct_h[0], conc1.S_h[0], conc1.S_v[0], conc1.S_c[0], CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

CD_w_trans = []
CD_fus_trans = []
CD_h_trans = []
CD_v_trans = []
CD_c_trans = []
CD_nacel_trans = []
CD_cruise = []

CD0_w1 = 0.006694693593553782
CD_flap_TO1, CD_flap_land1, CD_slat1 = config1_Drag.flaps_drag(CD0_w1)
CD_gear1 = config1_Drag.landinggear_drag()
CD_ws1 = config1_Drag.windshield_drag()
CD0_store1, CD_store_sub1, C_D_store_trans1 = config1_Drag.store_drag()
CD_trim1 = config1_Drag.trim_drag()
CD_spoiler1 = config1_Drag.spoiler_drag()

for i in range(len(CD_wave_w)):
    CD_w_sub1, CD_w_trans1, CD0_w1, e1 = config1_Drag.wing_drag(CD_wave_w[i])
    CD0_fus1, CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag(CD_wave_fus[i])
    CD0_h_tail1, CD0_v_tail1, CD0_c_tail1, CD_h_sub1, CD_v_sub1, CD_c_sub1, CD_h_trans1, CD_v_trans1, CD_c_trans1 = config1_Drag.empennage_drag(CD_wave_ch[i])
    CD0_nacel1, CD_nacel_sub1, CD_nacel_trans1 = config1_Drag.nacelle_drag(CD_wave_nac[i])
    
    CD_cruise1 = CD_w_trans1 + CD_fus_trans1 + CD_h_trans1 + CD_v_trans1 + CD_c_trans1 + CD_nacel_trans1 + CD_ws1 + C_D_store_trans1 + CD_trim1
    
    CD_w_trans.append(CD_w_trans1)
    CD_fus_trans.append(CD_fus_trans1)
    CD_h_trans.append(CD_h_trans1)
    CD_v_trans.append(CD_v_trans1)
    CD_c_trans.append(CD_c_trans1)
    CD_nacel_trans.append(CD_nacel_trans1)
    CD_cruise.append(CD_cruise1)
    
CD_TO1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_TO1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
CD_land1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_land1 + CD_slat1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
CD0_1 = CD0_w1 + CD0_fus1 + CD0_h_tail1 + CD0_v_tail1 + CD0_c_tail1 + CD0_nacel1 + CD0_store1

CD_cr = np.array(CD_cruise) / CD_cruise[2] * 100. - 100

""" Pitching Moment """
"""
config1_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[0], conc1.l_tail)
cm_des_airfoil1, dcm_dcl_airfoil1 = config1_Moment.Airfoil_moment()
delta_cm_flap1, delta_cm_krueger1 = config1_Moment.Airfoil_moment_flaps(cm_des_airfoil1)
Cm0_w_sub1, Cm0_w_trans1, dCm_dCl_w1 = config1_Moment.Wing_moment()
delta_Cm_w_flaps1, delta_Cm_w_krueger1 = config1_Moment.Wing_moment_flaps(Cm0_w_sub1)
"""


""" VALIDATION using Boeing 727-100"""

""" Lift """
alpha_0_l = -5.4
alpha_star_l = 10        # from -7 to 3 deg
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
C_l_alpha_M75 = C_l_alpha / math.sqrt(1 - const.M_cruise**2)
i_w= 0
wing_twist = -3         #DEG

Kcw = 0.95

S = 153
A = 7.67
l_f = 35.41
MAC = 1.75 * 40.60 / 15.45
Cr = 3.7 * 40.60 / 15.45
Ct = 1 * 40.60 / 15.45
b = math.sqrt(S*A)
taper_ratio = Ct / Cr
lambda_le_rad = 36 * math.pi / 180
lambda_4_rad = 32 * math.pi / 180
lambda_2_rad = 
A_h = 
A_c = 5
lambda_h_2_rad = 30 * math.pi / 180
lambda_c_2_rad = 0
i_c = 0 
S_h = 
S_c = 0
i_h =
x_le_MAC = 
b_flap =
SWF =
d_f_outer = 3.76
x_nlg = 
z_nlg = 
D_nlg = 
w_nlg = 
D_strut_nlg = 
x_mlg = 
z_mlg = 
D_mlg = 
w_mlg = 
D_strut_mlg = 
lambda_h_2_rad = 
lambda_v_2_rad = 
MAC_c = 0
Cr_v = 
Ct_v = 
Cr_h = 
Ct_h = 
S_v = 
CL_alpha_h = 
de_da = 
i_h = 
alpha0L_h = 
A_h = 
CL_alpha_c = 
de_da_c = 
i_c = 
alpha0L_c = 
A_c = 
l_fueltank = 0
d_fueltank = 0
delta_C_L_h = 
delta_C_L_c = 
S_elev = 
l_nacel = 
d_nacel = 
i_n = 
SWF = 
SWF_LE = 
Delta_C_L_flap = 
b_slat = 
b_flap = 
CLdes = 
CL_alpha = 
l_cockpit = 5.03 - 0.5*2.05
l_cabin = 
l_tail = 3 * 40.60 / 15.45


config1_Lift = aero.Lift(S,A,const.rho,const.rho_0,l_f,const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,MAC,Cr,Ct,b,taper_ratio,conc1.d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c, i_h, x_le_MAC, b_flap, SWF)
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1, delta_cl_flap_TO1, delta_clmax_TO1, clalpha_TO1 = config1_Lift.Airfoil_lift_flaps()
C_L_w1, CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1 = config1_Lift.Wing_lift()
delta_CL_w1, delta_CL_alpha_w1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1 = config1_Lift.Wing_lift_flaps(delta_cl_flap1,CL_alpha_w1,C_l_alpha,(delta_clmax_flap1 + delta_clmax_krueger1),b_slat, delta_cl_flap_TO1, delta_clmax_TO1)
CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1, de_da1, de_da_c1, alpha_CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)
delta_CL1, delta_CL_alpha1, delta_CL_max1, delta_CL_TO1, delta_CL_alpha_TO1, delta_CL_max_TO1 = config1_Lift.Airplane_lift_flaps(delta_CL_w1, CL_alpha_h1, CL_alpha_c1, delta_CL_alpha_w1, de_da1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1, Kcw)

""" Drag """

CD_wave_w = 0.002
CD_wave_fus = 0.005
CD_wave_ch = 0.002
CD_wave_nac = 0.002

config1_Drag = aero.Drag(S,A,const.rho,const.rho_0,l_f,const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes,CL_alpha,l_cockpit, l_cabin, l_tail, lambda_2_rad, lambda_4_rad, x_nlg, z_nlg, D_nlg, w_nlg, D_strut_nlg, x_mlg, z_mlg, D_mlg, w_mlg, D_strut_mlg, lambda_h_2_rad, lambda_v_2_rad, MAC_c, Cr_v, Ct_v, Cr_h, Ct_h, S_h, S_v, S_c, CL_alpha_h, de_da, i_h, alpha0L_h, conc1.A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
CD_w_sub1, CD_w_trans1, CD0_w1, e1 = config1_Drag.wing_drag(CD_wave_w)
CD0_fus1, CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag(CD_wave_fus)
CD0_h_tail1, CD0_v_tail1, CD0_c_tail1, CD_h_sub1, CD_v_sub1, CD_c_sub1, CD_h_trans1, CD_v_trans1, CD_c_trans1 = config1_Drag.empennage_drag(CD_wave_ch)
CD0_nacel1, CD_nacel_sub1, CD_nacel_trans1 = config1_Drag.nacelle_drag(CD_wave_nac)
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

