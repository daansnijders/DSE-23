# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:44:46 2019

@author: Niels
"""
from inputs.constants import *
from inputs.concept_1 import *
print (config1.lambda_le_rad)

def get_cg_wing(b,Cr,Ct):
    y_loc = 0.35
    dis = b/2*y_loc*sin(lambda_le_rad)
    location = 0.5*chordlength(y_loc, Cr, Ct)*0.7 + 0.25 *chordlength(y_loc, Cr, Ct)
    x_cg_wing = location+dis
    y_cg_wing = 0
    z_cg_wing = t_c*Cr/2
    return x_cg_wing, y_cg_wing, z_cg_wing

def get_cg_hwing(b_h,Ct_h,Cr_h):
    y_loc = 0.38
    dis = b_h/2*y_loc*sin(lambda_h_le_rad)
    location = 0.42*chordlength(y_loc, Cr_h, Ct_h)
    x_cg_hwing = location+dis
    y_cg_hwing = 0
    z_cg_hwing = 0.65*d_f_outer
    return x_cg_hwing, y_cg_hwing, z_cg_hwing

def get_cg_vwing(b_v,Ct_v,Cr_v):
    z_loc  = 0.38
    dis = b_v/2*z_loc*sin(lambda_v_le_rad)\
    location = 0.42*chordlength(z_loc, Cr_v, Ct_v)
    x_cg_hwing = location+dis
    y_cg_hwing = 0
    z_cg_hwing = d_f_outer+0.38*(b_v/2)
    return x_cg_vwing, y_cg_vwing, z_cg_vwing
    

def chordlength(y_loc, Cr, Ct):
    decrease = (Cr-Ct)/(b/2)
    return Cr-decrease*(b/2)*y_loc
    
#print (get_xcg_wing(10,1,2))