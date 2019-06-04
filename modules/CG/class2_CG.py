# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:09:22 2019

@author: Stijn
"""
from inputs.constants import *
from inputs.concept_1 import *
from modules.CG.CG_func import *

class get_cg(object):
    def __init__(self,b,Cr,Ct,t_c,lambda_le_rad,y_MAC,x_le_MAC,weights):
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
        self.weights       = weights
        
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
    
from modules.class2 import *
weights     = Class2_weight(1,N_pax[0],MTOW[0],M_carried_canard_MZF[0],min(M_MZF), n_max[0],V_dive[0],M_fuel[0], max(T_req), l_f[0],d_f_inner,d_f_outer,l_cabin[0], l_h[0], S, S_c[0], b, b_c[0], S_v[0],S_h[0],Cr_t,Cr_t_c[0],lambda_2_rad,lambda_h_2_rad[0], lambda_v_2_rad[0],lambda_c_2_rad, S_fus[0])     
    
    
a = get_cg(b,Cr,Ct,t_c,lambda_le_rad, y_MAC,x_le_MAC,weights)
