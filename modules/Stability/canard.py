# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:13:41 2019

@author: daansnijders
"""

from modules.Stability.Stability_runfile import *

class canard:
    def __init__ (self, weight_pass, config):
        self.weight_pass = weight_pass
        self.config = config
        self.additional_mass = weight_pass[config][-1] - weight_pass[0][-1]
        self.additional_weight = self.additional_mass * g
        
    def size_canard(self):
        # determine airfoil/angle of attack during cruise
        
        
        # determine force required/surface area
        S_c = self.additional_weight / (C_L_canard * 0.5 * rho * V_cruise**2)
        
        
        
        
        # determine location by use of the moment caused by the aditional module
        
        
        
        
        
        # determine aspect ratio/ taper ratio/ sweep/ ect.
        
        
    
