# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:32:22 2019

@author: Sybren
"""

from inputs.constants import *
from inputs.performance_inputs import *
from inputs.concept_1 import *
from modules.Aerodynamics import *


""" HLD design """
config1_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config1_HLD.HLD()
print(SWF, b_flap, SWF_LE, b_slat)

""" Drag classII estimations """
config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_to[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[0],CL_alpha,l_cockpit, l_cabin[0], l_tail)
CD0_fus1, CDL_fus1, CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag()
print(CD0_fus1, CDL_fus1,CD_fus_sub1, CD_fus_trans1)

config2_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_to[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[1],CL_alpha,l_cockpit, l_cabin[1], l_tail)
CD0_fus2, CDL_fus2,  CD_fus_sub2, CD_fus_trans2 = config2_Drag.fuse_drag()
print(CD0_fus2, CDL_fus2,CD_fus_sub2, CD_fus_trans2)

config3_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_to[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[2],CL_alpha,l_cockpit, l_cabin[2], l_tail)
CD0_fus3, CDL_fus3,  CD_fus_sub3, CD_fus_trans3 = config3_Drag.fuse_drag()
print(CD0_fus3, CDL_fus3,CD_fus_sub3, CD_fus_trans3)
>>>>>>> d74e35ad67a005be595c0f0b904a6718a5fec142
