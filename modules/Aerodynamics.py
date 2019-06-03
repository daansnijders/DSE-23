# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD:
    def __HLD__(self, Cl_TO, CLmax, S, A, lambda_4_rad, taper_ratio):
        Delta_CLmax = Cl_TO - CLmax
        lambda_hl_rad = np.arctan(np.tan(lambda_4_rad)-4/A*(13/20-1/4)(1-taper_ratio)/(1+taper_ratio))
        
        