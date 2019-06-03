# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:09:22 2019

@author: Stijn
"""
from inputs.constants import *
from inputs.concept_1 import *
from modules.CG.CG_func import *

class get_cg(object):
    def __init__(self,b,Cr,Ct,t_c,lambda_le_rad,y_MAC,x_le_MAC):
        self.x_cg          = None                                               # [m]
        self.y_cg          = None                                               # [m]
        self.z_cg          = None                                               # [m]
        self.b             = b                                                  # [m]
        self.Cr            = Cr                                                 # [m]
        self.Ct            = Ct                                                 # [m]
        self.t_c           = t_c                                                # [-]
        self.lambda_le_rad = lambda_le_rad                                      # [rad]
        self.y_MAC         = y_MAC                                              # [m]
        self.x_le_MAC      = x_le_MAC                                           # [m]
        
    def __str__(self):
        return 'ok'
        
    def calc_fuselage(self):
        pass
    
    def calc_wing(self):
        return get_cg_wing(self.b,self.Cr,self.Ct,self.t_c,self.lambda_le_rad,self.y_MAC,self.x_le_MAC)
    
    def calc_engine(self):
        pass
    
    def calc_empennage(self):
        pass
        
    def get_x_cg(self):
        return self.x_cg
    
    def get_y_cg(self):
        return self.y_cg
    
    def get_z_cg(self):
        return self.z_cg

cg = get_cg(b[0],Cr[0],Ct[0],t_c[0],lambda_le_rad[0],y_MAC[0],x_le_MAC[0])

print(cg.calc_wing())