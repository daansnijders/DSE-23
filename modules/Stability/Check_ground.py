# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:17:48 2019

@author: daansnijders
"""
from modules.loaddiagram_detailed import *
from modules.Stability.Stability_runfile import *

class Stability_check_ground:
    def __init__ (self, cg_min, cg_max, weight_max):
        self.cg_min=cg_min
        self.cg_max=cg_max
        self.weight_max=weight_max
        
               
    def check_equilibrium(self):
        
        
        
        return check