'inputs'
# from inputs.constants import *
# from inputs.performance_inputs import *
import math
import numpy as np
import inputs.constants as const
import inputs.performance_inputs as perf
import inputs.concept_1 as conc1
import matplotlib.pyplot as plt

'department modules'
import modules.Aerodynamics as aero
import modules.sustainability.greenhousegasemissions as gasemissions
import modules.sustainability.noise_defs as noise
import modules.initialsizing_loading as loadingdiagram
import modules.performance.class2_performance as class2performance
from Structure.Wing.wing_canard_iteration import wing_struc_analysis
import modules.initialsizing_undercarriage as initialunderc

'Iteration values'
import Output.read_load_variables as variab

print('check')



"""
AERODYNAMICS
"""

# initial guess for some values or needed from other departments, must be updated after iteration
i_w = 0                 #Angle of incidence wing in DEG
i_h = 0                 #Angle of incidence horinzontal tail in DEG
i_c = 0                 #Angle of incidence canard in DEG
i_n = 0                 #Angle of incidence nacelle in DEG
wing_twist = -3         #Wing twist in DEG
S_elev = 9.987716166863338
Delta_C_L_flap = 2.4 - 1.5          #Cl_land - Cl_clean from performanca_input

print('check')
'HLD'
#config1_HLD = aero.HLD_class(perf.Cl_land,perf.Cl_clean,conc1.S,conc1.A,conc1.lambda_4_rad,conc1.taper_ratio,conc1.CL_alpha,conc1.lambda_le_rad,conc1.Cr,conc1.d_f_outer)
config1_HLD = aero.HLD_class(variab.CL_land,variab.CL_clean_max,conc1.S,conc1.A,conc1.lambda_4_rad,conc1.taper_ratio,conc1.CL_alpha,conc1.lambda_le_rad,conc1.Cr,conc1.d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config1_HLD.HLD()


'Lift'
# initial values for lift
alpha_0_l = -5.4
alpha_star_l = 10        # from -7 to 3 deg
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
C_l_alpha_M75 = C_l_alpha / math.sqrt(1 - const.M_cruise**2)


#Set up different configurations
#config1_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[0], conc1.lambda_c_2_rad, i_c, conc1.S_h[0], conc1.S_c[0], i_h, conc1.x_le_MAC[0], b_flap, SWF)
#config2_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[1], conc1.lambda_c_2_rad, i_c, conc1.S_h[1], conc1.S_c[1], i_h, conc1.x_le_MAC[1], b_flap, SWF)
#config3_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, conc1.A_h, conc1.A_c,conc1.lambda_h_2_rad[2], conc1.lambda_c_2_rad, i_c, conc1.S_h[2], conc1.S_c[2], i_h, conc1.x_le_MAC[2], b_flap, SWF)


config1_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, variab.A_h, 0,variab.lambda_h_2_rad, 0.0, i_c, variab.S_h, 0., i_h, variab.x_le_MAC[0], b_flap, SWF)
config2_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, variab.A_h, variab.A_c2,variab.lambda_h_2_rad, variab.lambda_c_2_rad2, i_c, variab.S_h, variab.S_c2, i_h, variab.x_le_MAC[1], b_flap, SWF)
config3_Lift = aero.Lift(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, variab.A_h, variab.A_c3,variab.lambda_h_2_rad, variab.lambda_c_2_rad3, i_c, variab.S_h, variab.S_c3, i_h, variab.x_le_MAC[2], b_flap, SWF)




#Lift outputs
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1, delta_cl_flap_TO1, delta_clmax_TO1, clalpha_TO1 = config1_Lift.Airfoil_lift_flaps()
C_L_w1, CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1 = config1_Lift.Wing_lift()
delta_CL_w1, delta_CL_alpha_w1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1 = config1_Lift.Wing_lift_flaps(delta_cl_flap1,CL_alpha_w1,C_l_alpha,(delta_clmax_flap1 + delta_clmax_krueger1),b_slat, delta_cl_flap_TO1, delta_clmax_TO1)
CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1, de_da1, de_da_c1, alpha_CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)
delta_CL1, delta_CL_alpha1, delta_CL_max1, delta_CL_TO1, delta_CL_alpha_TO1, delta_CL_max_TO1 = config1_Lift.Airplane_lift_flaps(delta_CL_w1, CL_alpha_h1, CL_alpha_c1, delta_CL_alpha_w1, de_da1, delta_CL_max_w1, delta_CL_w_TO1, delta_CL_alpha_w_TO1, delta_CL_max_w_TO1)

delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2, delta_cl_flap_TO2, delta_clmax_TO2, clalpha_TO2 = config2_Lift.Airfoil_lift_flaps()
C_L_w2, CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2 = config2_Lift.Wing_lift()
delta_CL_w2, delta_CL_alpha_w2, delta_CL_max_w2, delta_CL_w_TO2, delta_CL_alpha_w_TO2, delta_CL_max_w_TO2 = config2_Lift.Wing_lift_flaps(delta_cl_flap2,CL_alpha_w2,C_l_alpha,(delta_clmax_flap2 + delta_clmax_krueger2),b_slat, delta_cl_flap_TO2, delta_clmax_TO2)
CL_alpha_h2, CL_alpha_c2, CL_alpha2, alpha_0_L2, CL_max2, de_da2, de_da_c2, alpha_CL_max2 = config2_Lift.Airplane_lift(CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2)
delta_CL2, delta_CL_alpha2, delta_CL_max2, delta_CL_TO2, delta_CL_alpha_TO2, delta_CL_max_TO2 = config2_Lift.Airplane_lift_flaps(delta_CL_w2, CL_alpha_h2, CL_alpha_c2, delta_CL_alpha_w2, de_da2, delta_CL_max_w2, delta_CL_w_TO2, delta_CL_alpha_w_TO2, delta_CL_max_w_TO2)

delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3, delta_cl_flap_TO3, delta_clmax_TO3, clalpha_TO3 = config3_Lift.Airfoil_lift_flaps()
C_L_w3, CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3  = config3_Lift.Wing_lift()
delta_CL_w3, delta_CL_alpha_w3, delta_CL_max_w3, delta_CL_w_TO3, delta_CL_alpha_w_TO3, delta_CL_max_w_TO3 = config3_Lift.Wing_lift_flaps(delta_cl_flap3,CL_alpha_w3,C_l_alpha,(delta_clmax_flap3 + delta_clmax_krueger3),b_slat, delta_cl_flap_TO3, delta_clmax_TO3)
CL_alpha_h3, CL_alpha_c3, CL_alpha3, alpha_0_L3, CL_max3, de_da3, de_da_c3, alpha_CL_max3 = config3_Lift.Airplane_lift(CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3)
delta_CL3, delta_CL_alpha3, delta_CL_max3, delta_CL_TO3, delta_CL_alpha_TO3, delta_CL_max_TO3 = config3_Lift.Airplane_lift_flaps(delta_CL_w3, CL_alpha_h3, CL_alpha_c3, delta_CL_alpha_w3, de_da3, delta_CL_max_w3, delta_CL_w_TO3, delta_CL_alpha_w_TO3, delta_CL_max_w_TO3)





'Drag'
#Values that must come from other departments

#D_strutt_nlg = conc1.D_strut_mlg                                                #UPDATE
#D_strutt_mlg = conc1.D_strut_nlg                                                #UPDATE
D_strutt_nlg = variab.D_strut_mlg                                                #UPDATE
D_strutt_mlg = variab.D_strut_nlg                                                #UPDATE
l_fueltank = 1.5                                                                #UPDATE
d_fueltank = 0.3                                                                #UPDATE


#Hard to calculate, not important for first iteration, maybe from Daan and Stijn
delta_C_L_h = 0                 #When i_h = 0
delta_C_L_c = 0                 #When i_c = 0

#Set up different configurations
#config1_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[0], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[0], conc1.lambda_v_2_rad[0], conc1.MAC_c[0], conc1.Cr_v[0], conc1.Ct_v[0], conc1.Cr_h[0], conc1.Ct_h[0], conc1.S_h[0], conc1.S_v[0], conc1.S_c[0], CL_alpha_h1, de_da1, i_h, alpha_0_L1, conc1.A_h, CL_alpha_c1, de_da_c1, i_c, alpha_0_L1, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#config2_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[1],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[1], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[1], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[1], conc1.lambda_v_2_rad[1], conc1.MAC_c[1], conc1.Cr_v[1], conc1.Ct_v[1], conc1.Cr_h[1], conc1.Ct_h[1], conc1.S_h[1], conc1.S_v[1], conc1.S_c[1], CL_alpha_h2, de_da2, i_h, alpha_0_L2, conc1.A_h, CL_alpha_c2, de_da_c2, i_c, alpha_0_L2, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#config3_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[2],const.V_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[2],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[2], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,conc1.x_nlg, conc1.z_nlg, const.D_nlg, const.w_nlg, conc1.D_strut_nlg, conc1.x_mlg[2], conc1.z_mlg, const.D_mlg, const.w_mlg, conc1.D_strut_mlg, conc1.lambda_h_2_rad[2], conc1.lambda_v_2_rad[2], conc1.MAC_c[2], conc1.Cr_v[2], conc1.Ct_v[2], conc1.Cr_h[2], conc1.Ct_h[2], conc1.S_h[2], conc1.S_v[2], conc1.S_c[2], CL_alpha_h3, de_da3, i_h, alpha_0_L3, conc1.A_h, CL_alpha_c3, de_da_c3, i_c, alpha_0_L3, conc1.A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

config1_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[0],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[0], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,variab.x_nlg, variab.z_nlg, const.D_nlg, const.w_nlg, variab.D_strut_nlg, variab.x_mlg[0], variab.z_mlg, const.D_mlg, const.w_mlg, variab.D_strut_mlg, variab.lambda_h_2_rad, variab.lambda_v_2_rad, 0.0, variab.Cr_v, variab.Ct_v, variab.Cr_h, variab.Ct_h, variab.S_h, variab.S_v, 0., CL_alpha_h1, de_da1, i_h, alpha_0_L1, variab.A_h, CL_alpha_c1, de_da_c1, i_c, alpha_0_L1, 0, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
config2_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[1],const.V_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[1],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[1], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,variab.x_nlg, variab.z_nlg, const.D_nlg, const.w_nlg, variab.D_strut_nlg, variab.x_mlg[1], variab.z_mlg, const.D_mlg, const.w_mlg, variab.D_strut_mlg, variab.lambda_h_2_rad, variab.lambda_v_2_rad, variab.MAC_c2, variab.Cr_v, variab.Ct_v,variab.Cr_h, variab.Ct_h, variab.S_h, variab.S_v, variab.S_c2, CL_alpha_h2, de_da2, i_h, alpha_0_L2, variab.A_h, CL_alpha_c2, de_da_c2, i_c, alpha_0_L2, variab.A_c2, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
config3_Drag = aero.Drag(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[2],const.V_cruise,conc1.V_TO[2],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.CLdes[2],conc1.CL_alpha,const.l_cockpit, conc1.l_cabin[2], conc1.l_tail, conc1.lambda_2_rad, conc1.lambda_4_rad,variab.x_nlg, variab.z_nlg, const.D_nlg, const.w_nlg, variab.D_strut_nlg, variab.x_mlg[2], variab.z_mlg, const.D_mlg, const.w_mlg, variab.D_strut_mlg, variab.lambda_h_2_rad, variab.lambda_v_2_rad, variab.MAC_c3, variab.Cr_v, variab.Ct_v, variab.Cr_h, variab.Ct_h, variab.S_h, variab.S_v, variab.S_c3, CL_alpha_h3, de_da3, i_h, alpha_0_L3, variab.A_h, CL_alpha_c3, de_da_c3, i_c, alpha_0_L3, variab.A_c3, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, conc1.l_nacel, conc1.d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

#Drag outputs
CD_w_sub1, CD_w_trans1, CD0_w1, e1 = config1_Drag.wing_drag()
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
CD_TO1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_TO1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1 + CD_nacel_sub1
CD_land1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_land1 + CD_slat1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1 + CD_nacel_sub1
CD0_1 = CD0_w1 + CD0_fus1 + CD0_h_tail1 + CD0_v_tail1 + CD0_c_tail1 + CD0_nacel1 + CD0_store1

CD_w_sub2, CD_w_trans2, CD0_w2, e2 = config2_Drag.wing_drag()
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
CD_TO2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_TO2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2 + CD_nacel_sub2
CD_land2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_land2 + CD_slat2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2 + CD_nacel_sub2
CD0_2 = CD0_w2 + CD0_fus2 + CD0_h_tail2 + CD0_v_tail2 + CD0_c_tail2 + CD0_nacel2 + CD0_store2

CD_w_sub3, CD_w_trans3, CD0_w3, e3 = config3_Drag.wing_drag()
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
CD_TO3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_TO3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3 + CD_nacel_sub3
CD_land3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_land3 + CD_slat3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3 + CD_nacel_sub3
CD0_3 = CD0_w3 + CD0_fus3 + CD0_h_tail3 + CD0_v_tail3 + CD0_c_tail3 + CD0_nacel3 + CD0_store3

'Moments'
#Initial values
x_ref = 0.5                     #Own decision based on nothing
cl_des_airfoil = 0.651          #CL at zero angle of attack from JAVAfoil (RE = 17*10^6, 65-615)

#Set up different configurations
config1_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[0],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[0], conc1.l_tail)
config2_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[1], conc1.l_tail)
config3_Moment = aero.Moment(conc1.S,conc1.A,const.rho,const.rho_0,conc1.l_f[0],const.V_cruise,const.M_cruise,conc1.V_TO[1],const.mu_37,const.mu_sl,conc1.MAC,conc1.Cr,conc1.Ct,conc1.b,conc1.taper_ratio,conc1.d_f_outer,conc1.lambda_le_rad,conc1.lambda_4_rad,conc1.lambda_2_rad, conc1.t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, conc1.y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat, const.l_cockpit, conc1.l_cabin[1], conc1.l_tail)




cm_des_airfoil1, dcm_dcl_airfoil1 = config1_Moment.Airfoil_moment()
delta_cm_flap1, delta_cm_krueger1 = config1_Moment.Airfoil_moment_flaps(cm_des_airfoil1)
Cm0_w_sub1, Cm0_w_trans1, dCm_dCl_w1 = config1_Moment.Wing_moment()
delta_Cm_w_flaps1, delta_Cm_w_krueger1 = config1_Moment.Wing_moment_flaps(Cm0_w_sub1)

cm_des_airfoil2, dcm_dcl_airfoil2 = config2_Moment.Airfoil_moment()
delta_cm_flap2, delta_cm_krueger2 = config2_Moment.Airfoil_moment_flaps(cm_des_airfoil2)
Cm0_w_sub2, Cm0_w_trans2, dCm_dCl_w2 = config2_Moment.Wing_moment()
delta_Cm_w_flaps2, delta_Cm_w_krueger2 = config2_Moment.Wing_moment_flaps(Cm0_w_sub2)

cm_des_airfoil3, dcm_dcl_airfoil3 = config3_Moment.Airfoil_moment()
delta_cm_flap3, delta_cm_krueger3 = config3_Moment.Airfoil_moment_flaps(cm_des_airfoil3)
Cm0_w_sub3, Cm0_w_trans3, dCm_dCl_w3 = config3_Moment.Wing_moment()
delta_Cm_w_flaps3, delta_Cm_w_krueger3 = config3_Moment.Wing_moment_flaps(Cm0_w_sub3)


CL_clean_max, C_L_flaps_unoptimized = config1_Lift.get_CL(CL_alpha1, alpha_0_L1, CL_max1, alpha_CL_max1, delta_CL1, delta_CL_alpha1, delta_CL_max1, 9.75)
CL_unoptiized, CL_flaps_max = config1_Lift.get_CL(CL_alpha1, alpha_0_L1, CL_max1, alpha_CL_max1, delta_CL1, delta_CL_alpha1, delta_CL_max1, 9.5)

#CL_flaps_max = CL_max1 + delta_CL_max1

print('AERO DONE')


CL_TO= (C_L_flaps_unoptimized+CL_clean_max)/2 # at an angle of 9.75 deg
CL_cruise1=conc1.CLdes[0]

CL_cruise2=conc1.CLdes[1]
CL_cruise3=conc1.CLdes[2]
CL_land= CL_flaps_max # at an angle of 9.5 deg


"""
PERFORMANCE & PROPULSION
"""
'simulation accuracy inputs'
max_airport_altitude = 2000.  # [m]
altitude_resolution = 5  # number of different altitudes considered
mass_resolution = 20  # resolution of plotting mass vs take-off field length

'general inputs'
engine_failure = False
show_performance_plots = True
show_airport_plots = False
show_rate_of_climb_plots = False

# temporary values, will be removed as soon as aerodynamics inputs are ready


#CL_TO = 1.9
#CL_land = 2.3
#C_D_0 = 0.018117539865047032
#CL_cruise = 0.8


lift_over_drag1 = CL_cruise1/CD_cruise1
lift_over_drag2 = CL_cruise2/CD_cruise2
lift_over_drag3 = CL_cruise3/CD_cruise3
aspect_ratio = conc1.A
oswald_efficiency_number = conc1.e

'analysis'

config1_Performance = class2performance.Performance(CL_TO, CL_land, CL_cruise1, CD0_1, CD_TO1, CD_land1, CD_cruise1, conc1.S, conc1.OEW[0],
                                  conc1.MTOW[0], const.g, perf.screen_height_to, perf.screen_height_la, perf.thrust_max,
                                  perf.friction_coefficient_to, perf.friction_coefficient_la,
                                  perf.reverse_thrust_factor, engine_failure, perf.thrust_setting_climb_out,
                                  perf.thrust_setting_transition, conc1.M_payload[0], conc1.M_fuel[0],
                                  max_airport_altitude, altitude_resolution, mass_resolution, perf.thrust_setting_climb,
                                  const.H_m, const.V_cruise, conc1.R[0], lift_over_drag1, aspect_ratio,
                                  oswald_efficiency_number, perf.correction_factor_to, show_performance_plots,
                                  show_airport_plots, perf.thrust_setting_descent, show_rate_of_climb_plots, 1)


print('check between conc')

config2_Performance =class2performance.Performance(CL_TO, CL_land, CL_cruise2, CD0_2, CD_TO2, CD_land2, CD_cruise2, conc1.S, conc1.OEW[1],
                                  conc1.MTOW[1], const.g, perf.screen_height_to, perf.screen_height_la, perf.thrust_max,
                                  perf.friction_coefficient_to, perf.friction_coefficient_la,
                                  perf.reverse_thrust_factor, engine_failure, perf.thrust_setting_climb_out,
                                  perf.thrust_setting_transition, conc1.M_payload[1], conc1.M_fuel[1],
                                  max_airport_altitude, altitude_resolution, mass_resolution, perf.thrust_setting_climb,
                                  const.H_m, const.V_cruise, conc1.R[1], lift_over_drag2, aspect_ratio,
                                  oswald_efficiency_number, perf.correction_factor_to, show_performance_plots,
                                  show_airport_plots, perf.thrust_setting_descent, show_rate_of_climb_plots, 2)

config3_Performance = class2performance.Performance(CL_TO, CL_land, CL_cruise3, CD0_3, CD_TO3, CD_land3, CD_cruise3, conc1.S, conc1.OEW[2],
                                  conc1.MTOW[2], const.g, perf.screen_height_to, perf.screen_height_la, perf.thrust_max,
                                  perf.friction_coefficient_to, perf.friction_coefficient_la,
                                  perf.reverse_thrust_factor, engine_failure, perf.thrust_setting_climb_out,
                                  perf.thrust_setting_transition, conc1.M_payload[2], conc1.M_fuel[2],
                                  max_airport_altitude, altitude_resolution, mass_resolution, perf.thrust_setting_climb,
                                  const.H_m, const.V_cruise, conc1.R[2], lift_over_drag3, aspect_ratio,
                                  oswald_efficiency_number, perf.correction_factor_to, show_performance_plots,
                                  show_airport_plots, perf.thrust_setting_descent, show_rate_of_climb_plots, 3)

'needed in further programs of iteration'
'fuel fractions'
Mff1 = config1_Performance.fuel_fraction_total
Mff2 = config2_Performance.fuel_fraction_total
Mff3 = config3_Performance.fuel_fraction_total

f1= config1_Performance.fuel_fraction_cruise_breguet
f2= config1_Performance.fuel_fraction_cruise_breguet
f3= config1_Performance.fuel_fraction_cruise_breguet
'Climb gradient'



'Climb velocity'
V_climb1=config1_Performance.take_off_velocity*1.3
V_climb2=config2_Performance.take_off_velocity*1.3
V_climb3=config3_Performance.take_off_velocity*1.3
print('P&P DONE')



"""
CONTROL AND STABILITY
"""
from Output.read_load_variables import x_le_MAC,l_h, x_mlg, l_n, l_m
from inputs.concept_1 import  MAC, S, y_MAC, lambda_2_rad, Cr, Ct, b, l_f,l_cutout
#from inputs.concept_1 import x_le_MAC,l_h, MAC, S, y_MAC, lambda_2_rad, Cr, Ct, b, l_f, x_mlg, l_cutout, l_n, l_m
import modules.initialsizing_planform as initialplanform

from modules.Stability.cg_weight_loadingdiagram import cg1_pass, cg2_pass, cg1_fuel, cg2_fuel, weight_fuel
from modules.Stability.control_surf_func import get_c_elev, get_S_elev, get_b_elev, get_c_rud, get_S_rud, get_b_rud, get_c_ail, get_S_ail, get_b_ail, get_c_splr, get_b_splr
from modules.Stability.check_ground import update_x_mlg, update_z_mlg, update_y_mlg, check_ground
from modules.Stability.cg_weight_loadingdiagram import  weight_pass, x_cg_min_flight1, x_cg_max_flight1, x_cg_max_flight2, x_cg_max_flight3
from modules.Stability.empennage import empennage
#from modules.testfile_aero import CL_alpha_h1, CL_alpha_w1, de_da, CL_max_w1, CL_alpha_c2, alpha_0_l
#from Output.class2_integration import config1_cg, config2_cg, config3_cg
from Output.class2_integration import config1_cg, config2_cg, config3_cg

'DELETE CONC1 BEFORE X LE MAC'
"""NEED FROM OTHER FILES"""
V_critical = min(config1_Performance.decision_speed,config1_Performance.approach_velocity )
#V_critical = 1.2                                                                        # [m/s] V1 speed/V_app
etah       = 0.9                                                                # [-] eta_h of aerodynamics
x_ac      = [x_le_MAC[0]+0.25*MAC, x_le_MAC[1]+0.25*MAC, x_le_MAC[2]+0.25*MAC]   # [m] x-location of the main wing ac
#CL_a_h = 2*np.pi
CL_a_h    = CL_alpha_h1                                                         # [-] CL_alpha_h
CL_a_ah   = CL_alpha_w1                                                         # [-] CL_alpha_(A-h)
de_da     = de_da1                                                              # [-] downwash
Vh_V      = 1.                                                                  # [-] V_h/V velocity factors
Cm_ac     = Cm0_w_trans1                                                        # [-] moment coefficient of main wing ac
#Cm_ac       = 0.2
CL_ah     = CL_max_w1                                                           # [-] CL_(A-h)
x_cg      = x_cg_max_flight1                                                    # [m] x-location of the most aft cg location for configuration 1 during flight
CL_h      = -0.8                                                                # [-] lift coefficient htail

CL_c      = 1.4                                                                 # [-] lift coefficient canard

CL_a_c    = CL_alpha_c2                                                         # [-] CL_alpha_canard



a_0       = alpha_0_l*np.pi/180                                                 # [rad] zero lift angle of attack
CN_h_a    = CL_a_h                                                              # [-] C_N_h_alpha htail
CN_w_a    = CL_alpha_w1                                                         # [-] C_N_w_alpha main wing
CN_c_a    = CL_a_c                                                              # [-] C_N_c_alpha canard
CN_h_def  = 0.5                                        #zelf                    # [-] C_N_h_de elevator deflection
Vc_V      = 1.                                          #zelf                   # [-] V_c/V velocity factors
"""====================="""


# initialize class:
empennage1 = empennage(2, x_ac, CL_a_h, CL_a_ah, de_da, l_h, conc1.S, conc1.MAC, Vh_V, x_le_MAC[0], Cm_ac, CL_ah, x_cg, CL_h, CL_c, CL_a_c, a_0, i_h, i_c, CN_h_a, CN_w_a, CN_c_a, CN_h_def, Vc_V, V_critical)
empennage2 = empennage(3, x_ac, CL_a_h, CL_a_ah, de_da, l_h, conc1.S, conc1.MAC, Vh_V, x_le_MAC[0], Cm_ac, CL_ah, x_cg, CL_h, CL_c, CL_a_c, a_0, i_h, i_c, CN_h_a, CN_w_a, CN_c_a, CN_h_def, Vc_V, V_critical)
if empennage2.Sc_S >= empennage1.Sc_S:
    empennage1.Sc_S = empennage2.Sc_S
else:
    empennage2.Sc_S = empennage1.Sc_S

# outputs:
x_le_MAC        = empennage1.x_le_MAC_out                                       # [m] x-location of MAC main wing
x_ac            = empennage1.x_ac                                                  # [m] x-location of the aerodynamic centre of the main wing
x_le_MAC_l_f    = empennage1.x_le_MAC_l_f                                       # [-] xlemac over fuselage length
x_le_w          = initialplanform.get_le_wing(conc1.y_MAC,x_le_MAC, conc1.lambda_2_rad, conc1.MAC, conc1.Cr)            # [m] x-location of le main wing with updated lemac

S_h             = empennage1.S_h                                                # [m^2] surface area of htail
A_h             = empennage1.A_h                                                # [-] aspect ratio htail
b_h             = empennage1.b_h                                                # [m] span htail
Cr_h            = empennage1.Cr_h                                               # [m] root chord htail
Ct_h            = empennage1.Ct_h                                               # [m] tip chord htail
taper_ratio_h   = empennage1.taper_ratio_h                                      # [-] taper ratio htail
lambda_h_le_rad = empennage1.lambda_h_le_rad                                    # [rad] leading edge sweep htail
lambda_h_2_rad  = empennage1.lambda_h_2_rad                                     # [rad] half chord sweep htail
lambda_h_4_rad  = empennage1.lambda_h_4_rad                                     # [rad] quarter chord sweep htail
x_h             = empennage1.x_h                                                # [m] x-location of ac of the htail?
l_h             = empennage1.l_h                                                # [m] distance between c/4 on MAC of the main wing and horizontal tail

S_v             = empennage1.S_v                                                # [m^2] surface area of vtail
A_v             = empennage1.A_v                                                # [-] aspect ratio vtail
b_v             = empennage1.b_v                                                # [m] span vtail
Cr_v            = empennage1.Cr_v                                               # [m] root chord vtail
Ct_v            = empennage1.Ct_v                                               # [m] tip chord vtail
taper_ratio_v   = empennage1.taper_ratio_v                                      # [-] taper ratio vtail
lambda_v_le_rad = empennage1.lambda_v_le_rad                                    # [rad] leading edge sweep vtail
lambda_v_2_rad  = empennage1.lambda_v_2_rad                                     # [rad] half chord sweep vtail
lambda_v_4_rad  = empennage1.lambda_v_4_rad                                     # [rad] quarter chord sweep vtail
x_v             = empennage1.x_v                                                # [m] x-location of ac of the vtail?

taper_ratio_c2  = empennage1.taper_ratio_c                                      # [-] taper ratio canard
lambda_c_le_rad2= empennage1.lambda_c_le_rad                                    # [rad] leading edge sweep angle canard

t_c_c2          = empennage1.t_c_c                                              # [-] tickness over chord ratio canard
Sc_S2           = empennage1.Sc_S                                               # [-] Ratio area canard (assumed for now)
S_c2            = empennage1.S_c                                                # [m^2] Surface area of the canard
A_c2            = empennage1.A_c                                                # [-] Aspect ratio of the canard
b_c2            = empennage1.b_c                                                # [m] span canard
Cr_c2           = empennage1.Cr_c                                               # [m] root chord length canard
Ct_c2           = empennage1.Ct_c                                               # [m] tip chord length canard
Cr_t_c2         = Cr_c2*t_c_c2                                                  #[m] thickness at the chord canard
z_c2            = empennage1.z_c                                                # [m] veritcal height of the canard
l_c2            = empennage1.l_c                                                # [m] distance 0.25mac-wing to 0.25MAC canard
lambda_c_2_rad2 = empennage1.lambda_c_2_rad
MAC_c2          = empennage1.MAC_c

taper_ratio_c3  = empennage2.taper_ratio_c                                      # [-] taper ratio canard
lambda_c_le_rad3= empennage2.lambda_c_le_rad                                    # [rad] leading edge sweep angle canard
t_c_c3          = empennage2.t_c_c                                              # [-] tickness over chord ratio canard
Sc_S3           = empennage2.Sc_S                                               # [-] Ratio area canard (assumed for now)
S_c3            = empennage2.S_c                                                # [m^2] Surface area of the canard
A_c3            = empennage2.A_c                                                # [-] Aspect ratio of the canard
b_c3            = empennage2.b_c                                                # [m] span canard
Cr_c3           = empennage2.Cr_c                                               # [m] root chord length canard
Cr_t_c3         = Cr_c3*t_c_c3                                                  #[m] thickness at the chord canard
Ct_c3           = empennage2.Ct_c                                               # [m] tip chord length canard
z_c3            = empennage2.z_c                                                # [m] veritcal height of the canard
l_c3            = empennage2.l_c                                                # [m] distance 0.25mac-wing to 0.25MAC canard
lambda_c_2_rad3 = empennage2.lambda_c_2_rad
MAC_c3          = empennage2.MAC_c

# control surfaces: (inputs still need to be worked on...)
c_elev = get_c_elev(Cr_h, Ct_h, b_h)                                            # [m] chord length elevator
S_elev = get_S_elev(S_h)                                                        # [m^2] surface area elevator
b_elev = get_b_elev(S_elev,c_elev)                                              # [m] span elevator

c_rud = get_c_rud(Cr_v, Ct_v, b_v)                                              # [m] chord length rudder
S_rud = get_S_rud(S_v)                                                          # [m^2] surface area rudder
b_rud = get_b_rud(S_rud,c_rud)                                                  # [m] span rudder

c_ail = get_c_ail(conc1.Cr,conc1.Ct,conc1.b)                                    # [m] chord length aileron
S_ail = get_S_ail(conc1.S)                                                      # [m^2] surface area aileron
b_ail = get_b_ail(conc1.b)                                                      # [m] span aileron

c_splr = get_c_splr(conc1.Cr, conc1.Ct, conc1.b)                                # [m] chord length spoiler
b_splr = get_b_splr(conc1.b)                                                    # [m] span spoiler


# Update cg's DIFFERENT CONFIG'S
# landing gear placement
x_mlg[0] = update_x_mlg(config1_cg.calc_z_cg(),const.theta_rad,const.beta_rad, x_cg_max_flight1, const.stroke,conc1.l_f[0]) # [m] x-location of the mlg
x_mlg[1] = max([x_mlg[0] + conc1.l_cutout, update_x_mlg(config2_cg.calc_z_cg(),const.theta_rad,const.beta_rad, x_cg_max_flight2, const.stroke,conc1.l_f[1])])
x_mlg[2] = max([x_mlg[0] + conc1.l_cutout, update_x_mlg(config3_cg.calc_z_cg(),const.theta_rad,const.beta_rad, x_cg_max_flight3, const.stroke,conc1.l_f[2])])

z_mlg = update_z_mlg(x_mlg[0],const.beta_rad,x_cg_max_flight1, config1_cg.calc_z_cg()) # [m] z-location of the mlg
z_nlg = z_mlg

x_nlg = 2

l_m1 = x_cg_max_flight1 - x_mlg[0]
l_m2 = x_cg_max_flight2 - x_mlg[1]
l_m3 = x_cg_max_flight3 - x_mlg[2]


l_n1 = x_cg_max_flight1 - x_nlg
l_n2 = x_cg_max_flight2 - x_nlg
l_n3 = x_cg_max_flight3 - x_nlg

moment = x_cg_min_flight1 * conc1.MTOW[0]
F_h_req = moment / (x_h-x_cg_min_flight1) 
assert S_h > F_h_req / (0.5*1.225*58.61409173313211**2 * 0.860)

y_mlg = update_y_mlg(config1_cg.calc_z_cg(),z_mlg,l_n,l_m)            # [m] y-location of the mlg

config1_ground      = check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])
config2_ground      = check_ground(cg1_pass[1], cg2_pass[1], weight_pass[1], cg1_fuel[1], cg2_fuel[1], weight_fuel[1], x_nlg, x_mlg[1])
config3_ground      = check_ground(cg1_pass[2], cg2_pass[2], weight_pass[2], cg1_fuel[2], cg2_fuel[2], weight_fuel[2], x_nlg, x_mlg[2])


frac = np.ones((3,2))
frac[0,0], frac[0,1], frac1 = config1_ground.check_equilibrium()
frac[1,0], frac[1,1], frac2 = config2_ground.check_equilibrium()
frac[2,0], frac[2,1], frac2 = config3_ground.check_equilibrium()
load_nlg_strut=[frac[0,1]*conc1.MTOW[0]*9.81, frac[1,1]*conc1.MTOW[1]*9.81, frac[2,1]*conc1.MTOW[2]*9.81]
load_nlg_strut=max(load_nlg_strut)
load_mlg_strut=[(1-frac[0,0])/2*conc1.MTOW[0]*9.81, (1-frac[1,0])/2*conc1.MTOW[1]*9.81, (1-frac[1,0])/2*conc1.MTOW[2]*9.81]
load_mlg_strut=max(load_mlg_strut)
if load_mlg_strut/2>const.max_P_mlg:
    print('landing gear main tires cannot handle the pressure')
if load_nlg_strut/2> const.max_P_nlg:
    print('landing gear nose tires cannot handle the pressure')


L_strut_mlg= -z_mlg-const.D_mlg/2
L_strut_nlg= -z_nlg-const.D_nlg/2
D_strut_mlg= initialunderc.get_d_lg(load_mlg_strut,L_strut_mlg)
D_strut_nlg= initialunderc.get_d_lg(load_nlg_strut,L_strut_nlg)

print('CS DONE')

"""
STRUCTURES
"""

    #Calculate fuel mass available for storage in wings (0.75 of one side of the wing)
fuel_mass_available = wing_struc_analysis(CL_flaps_max,const.V_cruise,conc1.b,conc1.Cr,conc1.Ct,conc1.S)


"""
SUSTAINABILITY
"""

'Greenhouse Gas Emissions'


config1_emissions   =gasemissions.greenhousegas_emissions(config1_Performance,1)
config1_NOx         =config1_emissions.get_NOx_mass()
config1_CO2         =config1_emissions.get_CO2_per_passenger_per_km()


config2_emissions   =gasemissions.greenhousegas_emissions(config2_Performance,2)
config2_NOx         =config2_emissions.get_NOx_mass()
config2_CO2         =config2_emissions.get_CO2_per_passenger_per_km()

config3_emissions   =gasemissions.greenhousegas_emissions(config3_Performance,3)
config3_NOx         =config3_emissions.get_NOx_mass()
config3_CO2         =config3_emissions.get_CO2_per_passenger_per_km()



'Noise'
phi_observer=np.radians(1)
flap_deflection=np.radians(50)
c_flap_end=2*conc1.S/((1+conc1.taper_ratio)*conc1.b)*(1-(1-conc1.taper_ratio)/conc1.b*2*b_flap)
area_flap=b_flap*2*(0.35*conc1.Cr+0.35*c_flap_end)/2

r1,r2,r3,theta_1,theta_2,theta_3=noise.simulate_flight_path(config1_Performance.approach_velocity)

OSPL_dBA_tot_straight=noise.EPNdB_calculations(r2,theta_2,phi_observer,config1_Performance.approach_velocity, area_flap, b_flap,flap_deflection,b_slat )
OSPL_dBA_tot_up=noise.EPNdB_calculations(r1,theta_1,phi_observer,config1_Performance.approach_velocity, area_flap, b_flap,flap_deflection,b_slat )
OSPL_dBA_tot_down=noise.EPNdB_calculations(r3,theta_3,phi_observer,config1_Performance.approach_velocity, area_flap, b_flap,flap_deflection, b_slat )



config1_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,CL_TO,CL_cruise1,CL_land,V_climb1,perf.c,f1,perf.sigma, perf.TOP, CD0_1,conc1.A,conc1.e,1000,7000,100)
config2_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,CL_TO,CL_cruise2,CL_land,V_climb2,perf.c,f2,perf.sigma, perf.TOP, CD0_2,conc1.A,conc1.e,1000,7000,100)
config3_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,CL_TO,CL_cruise3,CL_land,V_climb3,perf.c,f3,perf.sigma, perf.TOP, CD0_3,conc1.A,conc1.e,1000,7000,100)

#
##'WRITE THE CSV FILE'
#output_file = open('output_detailedsizing.dat' ,  'w')
##    output_file.write('V_h =' + str(V_h) + '\n')
#output_file.write('A_h  =' + str(A_h) + '\n')
#output_file.write('taper_ratio_h =' + str(taper_ratio_h) + '\n')
#output_file.write('lambda_h_2_rad =' + str(lambda_h_2_rad) + '\n')
##    output_file.write('x_le_h =' + str(x_le_h) + '\n')
#output_file.write('S_h=' + str(S_h) + '\n')
#output_file.write('b_h =' + str(b_h) + '\n')
#output_file.write('Cr_h =' + str(Cr_h) + '\n')
#output_file.write('Ct_h =' + str(Ct_h) + '\n')
#
#
##    output_file.write('V_v =' + str(V_v) + '\n')
#output_file.write('A_v =' + str(A_v) + '\n')
#output_file.write('lambda_v_2_rad =' + str(lambda_v_2_rad) + '\n')
##   output_file.write('x_le_v =' + str(x_le_v) + '\n')
#output_file.write('S_v  =' + str(S_v) + '\n')
#output_file.write('b_v =' + str(b_v) + '\n')
#output_file.write('Cr_v =' + str(Cr_v) + '\n')
#output_file.write('Ct_v =' + str(Ct_v) + '\n')
#
#
##    output_file.write('V_v =' + str(V_v) + '\n')
#output_file.write('A_c2 =' + str(A_c2) + '\n')
#output_file.write('lambda_c_2_rad2 =' + str(lambda_c_2_rad2) + '\n')
##    output_file.write('x_le_v =' + str(x_le_v) + '\n')
#output_file.write('S_c2  =' + str(S_c2) + '\n')
#output_file.write('b_c2 =' + str(b_c2) + '\n')
#output_file.write('Cr_c2 =' + str(Cr_c2) + '\n')
#output_file.write('Ct_c2 =' + str(Ct_c2) + '\n')
#output_file.write('Cr_t_c2=' + str(Cr_t_c2) + '\n')
#
#
##    output_file.write('V_v =' + str(V_v) + '\n')
#output_file.write('A_c3 =' + str(A_c3) + '\n')
#output_file.write('lambda_c_2_rad3 =' + str(lambda_c_2_rad3) + '\n')
##    output_file.write('x_le_v =' + str(x_le_v) + '\n')
#output_file.write('S_c3  =' + str(S_c3) + '\n')
#output_file.write('b_c3 =' + str(b_c3) + '\n')
#output_file.write('Cr_c3 =' + str(Cr_c3) + '\n')
#output_file.write('Ct_c3 =' + str(Ct_c3) + '\n')
#output_file.write('Cr_t_c3=' + str(Cr_t_c3) + '\n')
#
#output_file.write('lh=' + str(l_h)+ '\n')
#output_file.write('xleMAC1=' + str(x_le_MAC[0])+ '\n')
#output_file.write('xleMAC2=' + str(x_le_MAC[1])+ '\n')
#output_file.write('xleMAC3=' + str(x_le_MAC[2])+ '\n')
#
#output_file.write('MAC canard 2=' + str(MAC_c2) + '\n')
#output_file.write('MAC canard 3=' + str(MAC_c3) + '\n')
#
#
#output_file.write('lc_conf2=' + str(l_c2)+ '\n')
#output_file.write('lc_conf3=' + str(l_c3)+ '\n')
#
#output_file.write('x loc main 1 =' + str(x_mlg[0]) + '\n')
#output_file.write('x loc main 2 =' + str(x_mlg[1]) + '\n')
#output_file.write('x loc main 3 =' + str(x_mlg[2]) + '\n')
#output_file.write('x loc nose ='   + str(x_nlg) + '\n')
#output_file.write('l m 1 =' + str(l_m1) + '\n')
#output_file.write('l m 2 =' + str(l_m2) + '\n')
#output_file.write('l m 3 =' + str(l_m3) + '\n')
#output_file.write('l n 1 =' + str(l_n1) + '\n')
#output_file.write('l n 2 =' + str(l_n2) + '\n')
#output_file.write('l n 3 =' + str(l_n3) + '\n')
#output_file.write('L main landing gear strut ='+ str(L_strut_mlg)+ '\n')
#output_file.write('L nose landing gear strut ='+ str(L_strut_nlg)+ '\n')
#output_file.write('D main landing gear strut ='+ str(D_strut_mlg)+ '\n')
#output_file.write('D nose landing gear strut ='+ str(D_strut_nlg)+ '\n')
#
#output_file.write('CL_clean_max=' + str(CL_clean_max)+'\n')
#output_file.write('CL_flaps_max=' + str(CL_flaps_max)+'\n')
#output_file.write('CL_TO=' + str(CL_TO)+'\n')
#output_file.write('CL_land=' + str(CL_land)+'\n')
#output_file.write('CL_cruise1=' + str(CL_cruise1)+'\n')
#output_file.write('CL_cruise2=' + str(CL_cruise2)+'\n')
#output_file.write('CL_cruise3=' + str(CL_cruise3)+'\n')
#output_file.write('CD_cruise1=' + str(CD_cruise1)+'\n')
#output_file.write('CD_cruise2=' + str(CD_cruise2)+'\n')
#output_file.write('CD_cruise3=' + str(CD_cruise3)+'\n')
#
#output_file.write('z_nlg='+str(z_nlg)+'\n')
#output_file.write('z_mlg='+str(z_mlg)+'\n')
#
#
##output_file.write('A=' + str(conc1.A)+'\n')
##output_file.write('S=' + str(conc1.S)+'\n')
##output_file.write('b=' + str(conc1.b)+'\n')
##output_file.write('l_cutout=' + str(conc1.l_cutout)+'\n')
##output_file.write('l fuselage 1=' + str(conc1.l_f[0])+'\n')
##output_file.write('l cabin 1=' + str(conc1.l_cabin[0])+'\n')
##output_file.write('l fuselage extended=' + str(conc1.l_f[1])+'\n')
##output_file.write('l cabin extended=' + str(conc1.l_cabin[1])+'\n')
#output_file.write('mff1=' + str(Mff1)+'\n')
#output_file.write('mff2=' + str(Mff2)+'\n')
#output_file.write('mff3=' + str(Mff3)+'\n')
#
#output_file.write('co2emissions 1 per pax=' + str(config1_CO2[0])+'\n')
#output_file.write('co2emissions 1 requirement=' + str(config1_CO2[1])+'\n')
#output_file.write('co2emissions 2 per pax=' + str(config2_CO2[0])+'\n')
#output_file.write('co2emissions 2 requirement=' + str(config2_CO2[1])+'\n')
#output_file.write('co2emissions 3 per pax=' + str(config3_CO2[0])+'\n')
#output_file.write('co2emissions 3 requirement=' + str(config3_CO2[1])+'\n')
#output_file.write('nox emissions 1=' + str(config1_NOx) + '\n')
#output_file.write('nox emissions 2=' + str(config2_NOx) + '\n')
#output_file.write('nox emissions 3=' + str(config3_NOx) + '\n')
#output_file.close()













