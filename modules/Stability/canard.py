# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:13:41 2019

@author: daansnijders
"""

from modules.Stability.cg_weight_loadingdiagram import *
from modules.main_class2 import * 
from modules.Stability.empennage import *

class canard():
    def __init__ (self, weight_pass, config, CL_c):
        self.config = config - 1                                                    # [-] configuration selection
        self.weight_pass = weight_pass                                          # [kg] mass increase per passenger
        self.additional_mass = weight_pass[self.config][-1] - weight_pass[0][-1] # [kg] mass difference between config 1 and 2/3
        self.additional_weight = self.additional_mass * g                       # [N] weight difference between config 1 and 2/3
        self.weight = self.weight_pass[self.config][-1] * g                     # [N] MTOW config
        self.CL_c = CL_c                                                        # [-] CL canard
        self.Vc_V = 1
      
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

        x_c = 7.5                                                               # [m] x-location of canard ac
        self.l_c = cg_x[self.config] - x_c                                      # [m] distance between canard ac and c.g.
        l_h = x_le_h[0] + l_cutout - cg_x[self.config]                          # [m] distance between htail ac and c.g.
        l_cg = (x_le_MAC[self.config] + 0.25*MAC) - cg_x[self.config]           # [m] distance between wing ac and c.g.
        z_e = cg_z[self.config] - z_engine                                      # [m] distance between engine and c.g.
        F_e = thrust_max                                                        # [N] maximum thrust produced by the engine
        
        w = self.weight_pass[0][-1] * g
        x_h = x_le_h[0]
        
        self.F_h = (-w*((x_le_MAC[0] + 0.25*MAC) - cg_x[0]) + thrust_max * (cg_z[0] - config1_cg.z_cg_engines))/-((x_le_MAC[0] + 0.25*MAC) - cg_x[0] - x_h + config1_cg_x)
        self.F_w = w - self.F_h 
        
        margin = 1E-8
        assert -margin <= -self.F_w * ((x_le_MAC[0] + 0.25*MAC) - cg_x[0]) - self.F_h * (x_h - cg_x[0]) + F_e * (cg_z[0] - z_engine) <= margin
        assert -margin <= self.F_w + self.F_h - w <= margin
        
        self.F_w = (self.l_c * (self.weight - self.F_h) - l_h * self.F_h + F_e * z_e) / (l_cg + self.l_c)
        #self.F_h = (self.l_c * (self.weight - self.F_w) - l_cg * self.F_w + F_e * z_e) / (l_h + self.l_c)

        self.F_c = -self.F_w + self.weight - self.F_h
        
        margin = 1E-8
        assert -margin <= self.F_c + self.F_w + self.F_h - self.weight <= margin
        assert -margin <= self.l_c * self.F_c - l_cg * self.F_w - l_h * self.F_h + F_e * z_e <= margin
        
        CL_h2 = self.F_h / (0.5*rho*V_cruise**2*e2.S_h)
        S_c = self.F_c / (0.5*rho*V_cruise**2*self.CL_c)

        # determine location by use of the moment caused by the aditional module
    
        
        # determine aspect ratio/ taper ratio/ sweep/ ect.


    def plot_stability_canard(self, plot = True):
        aa = 1/(self.CL_a_c / e2.CL_a_ah * self.l_c * self.Vc_V**2)
        bb = -e2.x_ac - e2.CL_a_h / e2.CL_a_ah * (1-e2.de_da) * e2.S_h_S * e2.l_h