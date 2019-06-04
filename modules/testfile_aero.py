# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:32:22 2019

@author: Sybren
"""

from inputs.constants import *
from inputs.performance_inputs import *
from inputs.concept_1 import *


""" HLD design """
config1_HLD = HLD_class(Cl_land, Cl_clean, S, A, lambda_4_rad, taper_ratio,CL_alpha)
SWF = HLD_class.HLD(config1_HLD)
print(SWF)

