# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD_class:
    def __init__(self, Cl_land, Cl_clean, S, A, lambda_4_rad, taper_ratio, CL_alpha, lambda_le_rad):
        self.Cl_land        = Cl_land
        self.Cl_clean       = Cl_clean
        self.S              = S
        self.A              = A
        self.lambda_4_rad   = lambda_4_rad
        self.taper_ratio    = taper_ratio
        self.CL_alpha_clean = CL_alpha
        self.lambda_le_rad  = lambda_le_rad
    
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
        
        Sprime_S = 1 + SWF_S * (c_prime - 1)
        
        CL_alpha_flapped = Sprime_S * self.CL_alpha_clean
        
        
        HLD_clearance = 0.5     #Clearance between fuselage and the HLD's 
        SWF_LE = (0.1*self.S)/(0.9*0.3*self.lambda_le_rad)
        return(SWF, SWF_LE)
        
class Drag:
    def __init__(self,S,A,rho,l_f,V_cruise,mu_37):
        self.S          = S
        self.A          = A
        self.rho        = rho
        self.l_f        = l_f
        self.V_cruise   = V_cruise
        self.mu_37      = mu_37
        
        
    def wing_drag(self):
        Re_f  = rho   * V_cruise * l_f / mu_37
        Re_f0 = rho_0 * V_TO     * l_f / mu_sl
        #These result in
        Rwf =           #(Figure 4.1)
        
        cos_lambda_c_2_rad = cos(lambda_c_2_rad)
        #This results in
        R_LS =          #(Figure 4.2)
        
        
        
        

