# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:13:41 2019

@author: daansnijders
"""

from modules.Stability.Stability_runfile import *

class canard:
    def __init__ (self, weight_pass, config):
        self.weight_pass = weight_pass
        self.config = config
        self.additional_weight = weight_pass[config][-1] - weight_pass[0][-1]
        
    def size_canard(self):
        
        
    
