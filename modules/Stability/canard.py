# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:13:41 2019

@author: daansnijders
"""

from modules.Stability.Stability_runfile import *
from modules.EXECUTE_FILE import * 

class canard():
    def __init__ (self, weight_pass, config):
        self.config = config - 1                                                    # [-] configuration selection
        self.weight_pass = weight_pass                                          # [kg] mass increase per passenger
        self.additional_mass = weight_pass[self.config][-1] - weight_pass[0][-1] # [kg] mass difference between config 1 and 2/3
        self.additional_weight = self.additional_mass * g                       # [N] weight difference between config 1 and 2/3
        self.weight = self.weight_pass[self.config][-1] * g                     # [N] MTOW config
        
      
        config1_cg.calc_x_cg()
        self.size_canard()
        
    def size_canard(self):
        # determine airfoil/angle of attack during cruise
        
        C_L_canard = 0.5
        # determine force required/surface area
        S_c = self.additional_weight / (C_L_canard * 0.5 * rho * V_cruise**2)   # [m^2] surface area required by canard
        
        self.L_canard = self.additional_weight
        
        cg = [config1_cg.calc_x_cg(), config2_cg.calc_x_cg(), config3_cg.calc_x_cg()] # [m] x-location of the c.g.
        moment_arm = cg[self.config] - cg[0]                                    # [m] distance between canard ac and c.g.
        
        moment = moment_arm * self.weight                                       # [N*m] moment caused by the c.g. shift
        
        self.l_c = moment / self.L_canard                                       # [m] required distance between canard and c.g.
        self.x_c = cg[self.config] - self.l_c                                   # [m] x-location of canard ac
        
        # determine location by use of the moment caused by the aditional module
        
        
        
        
        # determine aspect ratio/ taper ratio/ sweep/ ect.
        
        
c = canard(weight_pass,2)

