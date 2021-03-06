# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:44:46 2019

@author: Niels
"""
#from inputs.constants import *
#import inputs.concept_1 import 
import numpy as np

def chordlength(y_loc, Cr, Ct,b):
    decrease = (Cr-Ct)/(b/2)                                                    # [-] slope of decrease of the chord along the length
    return Cr-decrease*(b/2)*y_loc                                              # [m] length of the chord at a certain y-location

def get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,x_le_MAC):
    y_loc = 0.35                                                                # [-] percentage of half-span
    dis = b/2*y_loc*np.sin(lambda_le_rad)                                       # [m] y-loation of the c.g. of the wing
    location = 0.5*chordlength(y_loc, Cr, Ct,b)*0.7 + 0.25 * chordlength(y_loc, Cr, Ct,b) # [m] x-location of the c.g. wrt the wing
    x_cg_wingstart = x_le_MAC - y_MAC*np.sin(lambda_le_rad)                     # [m] x-location of the start of the wing at the c.g. y-location
    x_cg_wing = x_cg_wingstart + (location+dis)                                 # [m] x-location of the c.g. of the wing
    y_cg_wing = 0                                                               # [m] y-location of the c.g. of the wing
    z_cg_wing = t_c*Cr/2                                                        # [m] z-location fo the c.g. of the wing
    return x_cg_wing, y_cg_wing, z_cg_wing


def get_cg_hwing(b_h,Cr_h,Ct_h,lambda_h_le_rad,x_le_h,d_f_outer):
    y_loc = 0.38                                                                # [-] percentage of half-span
    dis = b_h/2*y_loc*np.sin(lambda_h_le_rad)                                   # [m] y-loation of the c.g. of the horizontal wing
    location = 0.42*chordlength(y_loc, Cr_h, Ct_h,b_h)                          # [m] x-location of the c.g. wrt the horizontal wing
    x_cg_hwing = location+dis+x_le_h                                            # [m] x-location of the c.g. of the horizontal wing
    y_cg_hwing = 0                                                              # [m] y-location of the c.g. of the horizontal wing
    z_cg_hwing = 0.65*d_f_outer                                                 # [m] z-location of the c.g. of the horizontal wing
    return x_cg_hwing, y_cg_hwing, z_cg_hwing


def get_cg_vwing(b_v,Cr_v,Ct_v,lambda_v_le_rad,x_le_v,d_f_outer):
    z_loc  = 0.38                                                               # [-] percentage of half-span
    dis = b_v/2*z_loc*np.sin(lambda_v_le_rad)                                   # [m] y-loation of the c.g. of the vertical wing
    location = 0.42*chordlength(z_loc, Cr_v, Ct_v,b_v)                          # [m] x-location of the c.g. wrt the vertical wing
    x_cg_vwing = location+dis+x_le_v                                            # [m] x-location of the c.g. of the vertical wing
    y_cg_vwing = 0                                                              # [m] y-location of the c.g. of the vertical wing
    z_cg_vwing = d_f_outer+0.38*(b_v/2)                                         # [m] z-location of the c.g. of the vertical wing
    return x_cg_vwing, y_cg_vwing, z_cg_vwing
    
def get_cg_fuselage(l_f,d_f_outer):
    x_cg_fuselage = 0.435 * l_f                                                 # [m] x-location of the c.g. of the fuselage
    y_cg_fuselage = 0                                                           # [m] y-location of the c.g. of the fuselage
    z_cg_fuselage = 0.5*d_f_outer                                               # [m] z-location of the c.g. of the fuselage
    return x_cg_fuselage, y_cg_fuselage, z_cg_fuselage

def get_cg_canard(Cr_c,t_c_c,l_cutout,l_cockpit):
    x_cg_canard=l_cockpit+l_cutout/2 # add the gallery length                   # [m] x-location of the c.g. of the canard
    y_cg_canard=0                                                               # [m] y-location of the c.g. of the canard
    z_cg_canard=Cr_c*t_c_c/2                                                    # [m] z-location of the c.g. of the canard
    return x_cg_canard, y_cg_canard, z_cg_canard

def get_cg_landinggear_main(z_mlg,x_mlg):
    x_cg_landinggear_main=x_mlg                                                 # [m] x-location of the c.g. of the mlg
    y_cg_landinggear_main=0                                                     # [m] y-location of the c.g. of the mlg
    z_cg_landinggear_main=z_mlg/2                                               # [m] z-location of the c.g. of the mlg
    
    return x_cg_landinggear_main, y_cg_landinggear_main, z_cg_landinggear_main

def get_cg_landinggear_nose(z_nlg,x_nlg):
    x_cg_landinggear_nose=x_nlg                                                 # [m] x-location of the c.g. of the nlg
    y_cg_landinggear_nose=0                                                     # [m] y-location of the c.g. of the nlg
    z_cg_landinggear_nose=z_nlg/2                                               # [m] z-location of the c.g. of the nlg
    
    return x_cg_landinggear_nose, y_cg_landinggear_nose, z_cg_landinggear_nose
"""
NOT APPLICABLE YET

def get_cg_ngear():
    x_cg_ngear = x_nlg
    y_cg_ngear = 0
    z_cg_ngear = 0.5*z_nlg
    return x_cg_ngear, y_cg_ngear, z_cg_ngear

def get_cg_mgear():
    x_cg_mgear = x_mlg
    y_cg_mgear = 0
    z_cg_mgear = 0.5*z_mlg
    return x_cg_mgear, y_cg_mgear, z_cg_mgear
""" 
def get_cg_engines(x_cg_eng):       
    x_cg_engines = x_cg_eng                                                     # [m] x-location of the c.g. of the two engines
    y_cg_engines = 0                                                            # [m] y-location of the c.g. of the two engines
    z_cg_engines = 0                                                            # [m] z-location of the c.g. of the two engines
    return x_cg_engines, y_cg_engines, z_cg_engines
  
