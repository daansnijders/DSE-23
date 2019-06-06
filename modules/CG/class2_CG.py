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
        self.x_cg_wing,self.y_cg_wing,self.z_cg_wing=get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,self.x_le_MAC[self.config])
        self.x_cg_fuselage,self.y_cg_fuselage,self.z_cg_fuselage=get_cg_fuselage(l_f[self.config],d_f_outer)
        self.x_cg_htail,self.y_cg_htail,self.z_cg_htail=get_cg_hwing(b_h[self.config],Cr_h[self.config],Ct_h[self.config],lambda_h_le_rad,x_le_h[self.config],d_f_outer)
        self.x_cg_vtail,self.y_cg_vtail,self.z_cg_vtail= get_cg_vwing(b_v[self.config],Cr_v[self.config],Ct_v[self.config],lambda_v_le_rad,x_le_v[self.config],d_f_outer)
        self.x_cg_engines,self.y_cg_engines,self.z_cg_engines=get_cg_engines(x_cg_eng[self.config])
     
        
    def calc_x_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing *  self.x_cg_wing
        fuselage = self.weights.M_fuselage *  self.x_cg_fuselage
        h_tail = self.weights.M_horizontaltail *  self.x_cg_htail
        v_tail = self.weights.M_verticaltail *  self.x_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) * self.x_cg_engines
        self.x_cg_winggroup=sum([wing,engine])/(sum([self.weights.M_wing,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.x_cg_fuselagegroup=sum([fuselage,h_tail,v_tail])/sum([self.weights.M_fuselage,self.weights.M_horizontaltail,self.weights.M_verticaltail])
        # [m] returning location of cg 
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + self.weights.M_engines_total]))
    
    def calc_y_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * self.y_cg_wing
        fuselage = self.weights.M_fuselage *  self.y_cg_fuselage
        h_tail = self.weights.M_horizontaltail * self.y_cg_htail
        v_tail = self.weights.M_verticaltail *self.y_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) *self.y_cg_engines
        self.y_cg_winggroup=sum([wing,engine])/(sum([self.weights.M_wing,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.y_cg_fuselagegroup=sum([fuselage,h_tail,v_tail])/sum([self.weights.M_fuselage,self.weights.M_horizontaltail,self.weights.M_verticaltail])
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle +self.weights.M_engines_total]))
    
    def calc_z_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * self.z_cg_wing
        fuselage = self.weights.M_fuselage * self.z_cg_fuselage
        h_tail = self.weights.M_horizontaltail * self.z_cg_htail
        v_tail = self.weights.M_verticaltail * self.z_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) * self.z_cg_engines
        self.z_cg_winggroup=sum([wing,engine])/(sum([self.weights.M_wing,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.z_cg_fuselagegroup=sum([fuselage,h_tail,v_tail])/sum([self.weights.M_fuselage,self.weights.M_horizontaltail,self.weights.M_verticaltail])
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine])/(sum([self.weights.M_wing, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + self.weights.M_engines_total]))
    
    
   