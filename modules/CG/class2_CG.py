# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:09:22 2019

@author: Stijn
"""
from inputs.constants import *

class get_cg(object):
    def __init__(self,x_cg = None,y_cg = None,z_cg = None):
        self.x_cg = x_cg                                                        # [m]
        self.y_cg = y_cg                                                        # [m]
        self.z_cg = z_cg                                                        # [m]
        
    def calc_fuselage(self):
        pass
    
    def calc_wing(self):
        pass
    
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

cg = get_cg()