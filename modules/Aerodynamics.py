# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD_class:
    def __init__(self, Cl_TO, CLmax, S, A, lambda_4_rad, taper_ratio):
        self.Cl_TO          = Cl_TO
        self.CLmax          = CLmax
        self.S              = S[0]
        self.A              = A[0]
        self.lambda_4_rad   = lambda_4_rad[0]
        self.taper_ratio    = taper_ratio[0]
        
    
    def HLD(self):
        Delta_CLmax = [-i + self.Cl_TO for i in self.CLmax]
        hl = 0.65              # Location hinge line on chord
        lambda_hl_rad = np.arctan(np.tan(self.lambda_4_rad)-(4/self.A)*(hl-1/4)*(1-self.taper_ratio)/(1+self.taper_ratio))
        c_prime = 1 + 0.57*(1-hl)
        Delta_cl = 1.6*(c_prime)
        SWF = [(x *self.S)/(0.9*Delta_cl*np.cos(lambda_hl_rad)) for x in Delta_CLmax]
        SWF_S = SWF[0]/self.S
        delta_alpha = deg2rad(-15) * SWF_S * np.cos(lambda_hl_rad)
        
        HLD_clearance = 0.5     #Clearance between fuselage and the HLD's 
        
        return(SWF)