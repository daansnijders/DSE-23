# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:01:30 2019

@author: Lisa
"""
'inputs'
from inputs.constants import *
from inputs.performance_inputs import *
from inputs.concept_1 import *

'department modules'
from modules.Aerodynamics import *
from modules.sustainability.greenhousegasemissions import *



'AERODYNAMICS'
'Initial guess for some values, must be updated after iteration' 
i_w = 0                 #Angle of incidence wing in DEG
i_h = 0                 #Angle of incidence horinzontal tail in DEG
i_c = 0                 #Angle of incidence canard in DEG
wing_twist = -3         #Wing twist in DEG

'HLD'
config_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config_HLD.HLD()


'Lift'
#config1_Lift = Lift(S,A,rho,rho_0,l_f[0],V_cruise,M_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[0], lambda_c_2_rad, i_c, S_h[0], S_c[0], i_h, x_le_MAC[0], b_flap, SWF)
#config2_Lift = Lift(S,A,rho,rho_0,l_f[1],V_cruise,M_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[1], lambda_c_2_rad, i_c, S_h[1], S_c[1], i_h, x_le_MAC[1], b_flap, SWF)
#config3_Lift = Lift(S,A,rho,rho_0,l_f[2],V_cruise,M_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[2], lambda_c_2_rad, i_c, S_h[2], S_c[2], i_h, x_le_MAC[2], b_flap, SWF)


'Drag'
'Values that must come from other departments'
D_strutt_nlg = 0.15
D_strutt_mlg = 0.2
#config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[0],CL_alpha,l_cockpit, l_cabin[0], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[0], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[0], lambda_v_2_rad[0], MAC_c[0], Cr_v[0], Ct_v[0], Cr_h[0], Ct_h[0], S_h[0], S_v[0], S_c[0], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#config2_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[1],CL_alpha,l_cockpit, l_cabin[1], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[1], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[1], lambda_v_2_rad[1], MAC_c[1], Cr_v[1], Ct_v[1], Cr_h[1], Ct_h[1], S_h[1], S_v[1], S_c[1], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
#config3_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[2],CL_alpha,l_cockpit, l_cabin[2], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[2], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[2], lambda_v_2_rad[2], MAC_c[2], Cr_v[2], Ct_v[2], Cr_h[2], Ct_h[2], S_h[2], S_v[2], S_c[2], CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)









'Performance and propulsion'









'Control and Stability'










'sustainability '
'greenhouse gas emissions'
#CO_2_emissions=get_CO2_emissions(M_fuel_burnt)
#
#EI_NOx_phase=get_EI_NOX_fuelflow(fuel_flow)
#    
#  
#NOx_phase = get_NOx_emissions_total(fuel_flow_phase,time_in_flight_phase,EI_NOx_phase)
# 
#    
#Dp_Foo_NOx=get_Dp_Foo_NOx_specific(NOx_total,T_phase)               
#    
#NO_x_req=get_NOx_reduction_CAEP(Dp_Foo_NOx_caep,Dp_Foo_NOx)  #should come to 45%
#     


'Structures'