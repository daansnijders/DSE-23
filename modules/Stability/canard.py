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
        
        # determine location by use of the moment caused by the aditional module
        cg_x = [config1_cg.calc_x_cg(), config2_cg.calc_x_cg(), config3_cg.calc_x_cg()] # [m] x-location of the c.g.
        cg_z = [config1_cg.calc_z_cg(), config2_cg.calc_z_cg(), config3_cg.calc_z_cg()] # [m] x-location of the c.g.

        moment_arm = (cg_x[self.config] - cg_x[0]  - l_cutout)*-1               # [m] distance between canard ac and c.g.
        moment = moment_arm * self.weight                                       # [N*m] moment caused by the c.g. shift
    
        self.l_c = moment / self.L_canard                                       # [m] required distance between canard and c.g.of config 1
        self.x_c = cg_x[0] + l_cutout - self.l_c                                # [m] x-location of canard ac
    
        
        x_c = 5                                                                 # [m] x-location of canard ac
        l_c = cg_x[self.config] - x_c
        l_h = x_le_h[self.config] - cg_x[self.config]
        F_h = 0
        l_cg = (x_le_MAC[self.config] + 0.25*MAC) - cg_x[self.config]
        z_e = cg_z[self.config] - z_engine
        F_e = thrust_max
        
        
        F_w = (l_c * (self.weight + F_h) + l_h * F_h + F_e * z_e) / (l_cg + l_c)
        F_c = -F_w + self.weight + F_h

        
        margin = 1E-8
        print(F_c, F_w, F_h)
        assert -margin <= F_c + F_w - self.weight - F_h <= margin
        assert -margin <= l_c * F_c - l_cg * F_w + l_h * F_h + F_e * z_e <= margin

        # determine location by use of the moment caused by the aditional module
    
        w = self.weight_pass[0][-1] * g
        x_h = 0.9* l_f[0]
        F_h = (-w*(config1_cg.x_cg_wing - config1_cg_x) + thrust_max * (config1_cg_z - config1_cg.z_cg_engines))/-(-config1_cg.x_cg_wing + config1_cg_x + x_h - config1_cg_x)
        F_w = w + F_h 

        
        # determine aspect ratio/ taper ratio/ sweep/ ect.
        
        
c = canard(weight_pass,2)

