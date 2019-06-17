# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:09:22 2019

@author: Stijn
"""
import inputs.constants as const
import inputs.concept_1 as c1
import modules.CG.CG_func as cgcomp

class get_cg(object):
    def __init__(self,x_le_MAC,weights):
        self.x_le_MAC      = x_le_MAC                                           # [m]
        self.weights       = weights                                            # [kg] class
        self.config        = self.weights.config-1                              # [-] 
        self.x_cg_wing,self.y_cg_wing,self.z_cg_wing=cgcomp.get_cg_wing( c1.b, c1.Cr, c1.Ct, c1.t_c, c1.lambda_le_rad, c1.y_MAC,self.x_le_MAC[self.config]) # [m] x,y,z-location of the main wing
        self.x_cg_fuselage,self.y_cg_fuselage,self.z_cg_fuselage=cgcomp.get_cg_fuselage( c1.l_f[self.config], c1.d_f_outer) # [m] x,y,z-location of the fuselage
        self.x_cg_htail,self.y_cg_htail,self.z_cg_htail=cgcomp.get_cg_hwing( c1.b_h[self.config], c1.Cr_h[self.config], c1.Ct_h[self.config], c1.lambda_h_le_rad, c1.x_le_h[self.config],c1.d_f_outer) # [m] x,y,z-location of the htail
        self.x_cg_vtail,self.y_cg_vtail,self.z_cg_vtail= cgcomp.get_cg_vwing( c1.b_v[self.config], c1.Cr_v[self.config], c1.Ct_v[self.config], c1.lambda_v_le_rad, c1.x_le_v[self.config], c1.d_f_outer) # [m] x,y,z-location of the vtail
        self.x_cg_engines,self.y_cg_engines,self.z_cg_engines=cgcomp.get_cg_engines( c1.x_cg_eng[self.config]) # [m] x,y,z-location of the engines
        self.x_cg_canard,self.y_cg_canard,self.z_cg_canard=cgcomp.get_cg_canard( c1.Cr_c[self.config], c1.t_c_c[self.config], c1.l_cutout, const.l_cockpit) # [m] x,y,z-location of the canard
        self.x_cg_landinggear_main, self.y_cg_landinggear_main, self.z_cg_landinggear_main = cgcomp.get_cg_landinggear_main( c1.z_mlg, c1.x_mlg[self.config]) # [m] x,y,z-location of the mlg
        self.x_cg_landinggear_nose, self.y_cg_landinggear_nose, self.z_cg_landinggear_nose = cgcomp.get_cg_landinggear_nose( c1.z_nlg, c1.x_nlg) # [m] x,y,z-location of the nlg
        
    def calc_x_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing *  self.x_cg_wing
        fuselage = self.weights.M_fuselage *  self.x_cg_fuselage
        h_tail = self.weights.M_horizontaltail *  self.x_cg_htail
        v_tail = self.weights.M_verticaltail *  self.x_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) * self.x_cg_engines
        landinggear_main= self.weights.M_landinggear_main*self.x_cg_landinggear_main
        landinggear_nose=self.weights.M_landinggear_nose*self.x_cg_landinggear_nose
        if self.weights.config==1:
            canard=0
        else:
            canard=self.weights.M_canard*self.x_cg_canard
        self.x_cg_winggroup=sum([wing,engine,landinggear_main])/(sum([self.weights.M_wing,self.weights.M_landinggear_main,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.x_cg_fuselagegroup=sum([fuselage,h_tail,v_tail,canard,landinggear_nose])/sum([self.weights.M_fuselage,self.weights.M_landinggear_nose,self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_canard])
        # [m] returning location of cg 
        return sum([wing,fuselage,h_tail,v_tail,engine,canard,landinggear_main,landinggear_nose])/(sum([self.weights.M_wing, self.weights.M_canard,self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + self.weights.M_engines_total,self.weights.M_landinggear_nose,self.weights.M_landinggear_main]))
    
    def calc_y_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * self.y_cg_wing
        fuselage = self.weights.M_fuselage *  self.y_cg_fuselage
        h_tail = self.weights.M_horizontaltail * self.y_cg_htail
        v_tail = self.weights.M_verticaltail *self.y_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) *self.y_cg_engines
        landinggear_main= self.weights.M_landinggear_main*self.y_cg_landinggear_main
        landinggear_nose=self.weights.M_landinggear_nose*self.y_cg_landinggear_nose
        if self.weights.config==1:
            canard=0
        else:
            canard=self.weights.M_canard*self.y_cg_canard
        self.y_cg_winggroup=sum([wing,engine,landinggear_main])/(sum([self.weights.M_wing,self.weights.M_landinggear_main,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.y_cg_fuselagegroup=sum([fuselage,h_tail,v_tail,canard,landinggear_nose])/sum([self.weights.M_fuselage,self.weights.M_landinggear_nose,self.weights.M_canard,self.weights.M_horizontaltail,self.weights.M_verticaltail])
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine,canard,landinggear_main,landinggear_nose])/(sum([self.weights.M_wing,self.weights.M_canard, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle +self.weights.M_engines_total,self.weights.M_landinggear_nose,self.weights.M_landinggear_main]))
    
    def calc_z_cg(self):                                                        # [kg*m] mass times c.g distance of the different groups. 
        wing = self.weights.M_wing * self.z_cg_wing
        fuselage = self.weights.M_fuselage * self.z_cg_fuselage
        h_tail = self.weights.M_horizontaltail * self.z_cg_htail
        v_tail = self.weights.M_verticaltail * self.z_cg_vtail
        engine = (self.weights.M_nacelle + self.weights.M_engines_total) * self.z_cg_engines
        landinggear_main= self.weights.M_landinggear_main*self.z_cg_landinggear_main
        landinggear_nose=self.weights.M_landinggear_nose*self.z_cg_landinggear_nose
        if self.weights.config==1:
            canard=0
        else:
            canard=self.weights.M_canard*self.z_cg_canard
        self.z_cg_winggroup=sum([wing,engine,landinggear_main])/(sum([self.weights.M_wing,self.weights.M_landinggear_main,self.weights.M_nacelle + self.weights.M_engines_total]))
        self.z_cg_fuselagegroup=sum([fuselage,h_tail,v_tail,canard,landinggear_nose])/sum([self.weights.M_fuselage,self.weights.M_landinggear_nose,self.weights.M_canard,self.weights.M_horizontaltail,self.weights.M_verticaltail])
        # [m] returning location of cg
        return sum([wing,fuselage,h_tail,v_tail,engine,canard,landinggear_main,landinggear_nose ])/(sum([self.weights.M_wing,self.weights.M_canard, self.weights.M_fuselage, self.weights.M_horizontaltail,self.weights.M_verticaltail,self.weights.M_nacelle + self.weights.M_engines_total,self.weights.M_landinggear_nose,self.weights.M_landinggear_main]))
    
    
   