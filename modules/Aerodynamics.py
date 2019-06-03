# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD:
    def __init__(self, Cl_TO, CLmax, S, A, lambda_4_rad, taper_ratio):
        self.Cl_TO          = Cl_TO
        self.CLmax          = CLmax
        self.S              = S
        self.A              = A
        self.lambda_4_rad   = lambda_4_rad
        self.taper_ratio    = taper_ratio
        
    
    def HLD(self):
        self.Delta_CLmax = Cl_TO - CLmax
        self.hl = 13/20              # Location hinge line on chord
        self.lambda_hl_rad = np.arctan(np.tan(lambda_4_rad)-4/A*(hl-1/4)(1-taper_ratio)/(1+taper_ratio))
        self.c_prime = 1 + 0.57*(1-hl)
        self.Delta_cl = 1.3*(c_prime)
        
        self.SWF = (Delta_CLmax*S)/(0.9*Delta_cl*np.cos(lambda_hl_rad))
        
        return(SWF)