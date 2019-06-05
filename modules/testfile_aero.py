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


""" HLD design """
config1_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config1_HLD.HLD()
print(SWF, b_flap, SWF_LE, b_slat)

""" Drag classII estimations """
config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[0],CL_alpha,l_cockpit, l_cabin[0], l_tail, lambda_2_rad, x_nlg[0], z_nlg[0], D_nlg, b_nlg, x_mlg[0], z_mlg[0], D_mlg, b_mlg, lambda_h_2_rad[0], lambda_v_2_rad[0], MAC_c[0], Cr_v[0], Ct_v[0], Cr_h[0], Ct_h[0])
CD0_fus1, CDL_fus1, CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag()
print(CD0_fus1, CDL_fus1,CD_fus_sub1, CD_fus_trans1)

config2_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[1],CL_alpha,l_cockpit, l_cabin[1], l_tail, lambda_2_rad, x_nlg[1], z_nlg[1], D_nlg, b_nlg, x_mlg[1], z_mlg[1], D_mlg, b_mlg, lambda_h_2_rad[1], lambda_v_2_rad[1], MAC_c[1], Cr_v[1], Ct_v[1], Cr_h[1], Ct_h[1])
CD0_fus2, CDL_fus2,  CD_fus_sub2, CD_fus_trans2 = config2_Drag.fuse_drag()
print(CD0_fus2, CDL_fus2,CD_fus_sub2, CD_fus_trans2)

config3_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[2],CL_alpha,l_cockpit, l_cabin[2], l_tail, lambda_2_rad, x_nlg[2], z_nlg[2], D_nlg, b_nlg, x_mlg[2], z_mlg[2], D_mlg, b_mlg, lambda_h_2_rad[2], lambda_v_2_rad[2], MAC_c[2], Cr_v[2], Ct_v[2], Cr_h[2], Ct_h[2])
CD0_fus3, CDL_fus3,  CD_fus_sub3, CD_fus_trans3 = config3_Drag.fuse_drag()
print(CD0_fus3, CDL_fus3,CD_fus_sub3, CD_fus_trans3)

