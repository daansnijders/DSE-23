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
config1_HLD = HLD_class(Cl_land, Cl_clean, S, A, lambda_4_rad, taper_ratio, CL_alpha, lambda_le_rad)
SWF, SWF_LE = config1_HLD.HLD()
print(SWF, SWF_LE)


""" Drag classII estimations """
config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_to[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad)
CD0_fus1 = config1_Drag.fuse_drag()
print(CD0_fus1)

config1_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_to[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad)
CD0_fus2 = config1_Drag.fuse_drag()
print(CD0_fus2)

config1_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_to[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad)
CD0_fus3 = config1_Drag.fuse_drag()
print(CD0_fus3)