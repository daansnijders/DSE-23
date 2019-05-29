# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:44:46 2019

@author: Niels
"""
from inputs.constants import *
from inputs.concept_1 import *

def get_xcg_wing(b,Ct,Cr):
    y_loc = 0.35
    dis = b/2*0.35*sin(lambda_le_rad)
    location = 0.5*chordlength(y_loc)*0.7 + 0.25 *chordlength(y_loc)
    location_from_LE = location+dis
    

def chordlength(y_loc):
    decrease = (Cr-Ct)/(b/2)
    return Cr-decrease*(b/2)*y_loc
    
#print (get_xcg_wing(10,1,2))