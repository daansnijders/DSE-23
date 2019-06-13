# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:01:30 2019

@author: Lisa
"""
from inputs.constants import *
from inputs.performance_inputs import *
from inputs.concept_1 import *
from modules.Aerodynamics import *



'AERODYNAMICS'
'HLD'
config_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config_HLD.HLD()













'Performance and propulsion'









'Control and Stability'










'Structures'