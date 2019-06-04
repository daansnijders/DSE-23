# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD_class:
    def __init__(self, Cl_land, Cl_clean, S, A, lambda_4_rad, taper_ratio):
        self.Cl_land        = Cl_land
        self.Cl_clean       = Cl_clean
        self.S              = S[0]
        self.A              = A[0]
        self.lambda_4_rad   = lambda_4_rad[0]
        self.taper_ratio    = taper_ratio[0]
        
    
    def HLD(self):
        Delta_CLmax = self.Cl_land - self.Cl_clean    #[-i + self.Cl_land for i in self.Cl_clean]
        hl = 0.65              # Location hinge line on chord
        lambda_hl_rad = np.arctan(np.tan(self.lambda_4_rad)-(4/self.A)*(hl-1/4)*(1-self.taper_ratio)/(1+self.taper_ratio))
        c_prime = 1 + 0.57*(1-hl)
        Delta_cl = 1.6*(c_prime)
        SWF = (Delta_CLmax *self.S)/(0.9*Delta_cl*np.cos(lambda_hl_rad))
        #[(x *self.S)/(0.9*Delta_cl*np.cos(lambda_hl_rad)) for x in Delta_CLmax]
        SWF_S = SWF/self.S
        delta_alpha = np.deg2rad(-15) * SWF_S * np.cos(lambda_hl_rad)
        print(SWF_S)
        
        HLD_clearance = 0.5     #Clearance between fuselage and the HLD's 
        
        return(SWF)