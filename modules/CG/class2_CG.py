# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:09:22 2019

@author: Stijn
"""
from inputs.constants import *
from inputs.concept_1 import *
from modules.CG.CG_func import *

class get_cg(object):
    def __init__(self,x_le_MAC,weights):
        self.x_le_MAC      = x_le_MAC                                           # [m]
        self.weights       = weights                                            # [kg] class
        self.config        = self.weights.config-1                              # [-]
    
    def calc_x_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        self.cg_wing = get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,self.x_le_MAC[self.config])[0]
        wing = self.weights.M_wing * self.cg_wing
        fuselage = self.weights.M_fuselage * get_cg_fuselage(l_f[self.config],d_f_outer)[0]
        h_tail = self.weights.M_horizontaltail * get_cg_hwing(b_h[self.config],Cr_h[self.config],Ct_h[self.config],lambda_h_le_rad,x_le_h[self.config],d_f_outer)[0]
        v_tail = self.weights.M_verticaltail * get_cg_vwing(b_v[self.config],Cr_v[self.config],Ct_v[self.config],lambda_v_le_rad,x_le_v[self.config],d_f_outer)[0]
        engine = (self.weights.M_nacelle + M_engine) * get_cg_engines(x_cg_eng[self.config])[0]
        
        # [m] returning location of cg 
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + M_engine]))
    
    def calc_y_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,x_le_MAC[self.config])[1]
        fuselage = self.weights.M_fuselage * get_cg_fuselage(l_f[self.config],d_f_outer)[1]
        h_tail = self.weights.M_horizontaltail * get_cg_hwing(b_h[self.config],Cr_h[self.config],Ct_h[self.config],lambda_h_le_rad,x_le_h[self.config],d_f_outer)[1]
        v_tail = self.weights.M_verticaltail * get_cg_vwing(b_v[self.config],Cr_v[self.config],Ct_v[self.config],lambda_v_le_rad,x_le_v[self.config],d_f_outer)[1]
        engine = (self.weights.M_nacelle + M_engine) * get_cg_engines(x_cg_eng[self.config])[1]
        
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + M_engine]))
    
    def calc_z_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,self.x_le_MAC[self.config])[2]
        fuselage = self.weights.M_fuselage * get_cg_fuselage(l_f[self.config],d_f_outer)[2]
        h_tail = self.weights.M_horizontaltail * get_cg_hwing(b_h[self.config],Cr_h[self.config],Ct_h[self.config],lambda_h_le_rad,x_le_h[self.config],d_f_outer)[2]
        v_tail = self.weights.M_verticaltail * get_cg_vwing(b_v[self.config],Cr_v[self.config],Ct_v[self.config],lambda_v_le_rad,x_le_v[self.config],d_f_outer)[2]
        engine = (self.weights.M_nacelle + M_engine) * get_cg_engines(x_cg_eng[self.config])[2]
        
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + M_engine]))