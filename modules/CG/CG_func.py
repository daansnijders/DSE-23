# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:44:46 2019

@author: Niels
"""
#from inputs.constants import *
from inputs.concept_1 import *
import numpy as np

def chordlength(y_loc, Cr, Ct,b):
    decrease = (Cr-Ct)/(b/2)
    return Cr-decrease*(b/2)*y_loc

def get_cg_wing(b,Cr,Ct,t_c,lambda_le_rad,y_MAC,x_le_MAC):
    y_loc = 0.35
    dis = b/2*y_loc*np.sin(lambda_le_rad)
    location = 0.5*chordlength(y_loc, Cr, Ct,b)*0.7 + 0.25 * chordlength(y_loc, Cr, Ct,b)
    x_cg_wingstart = x_le_MAC - y_MAC*np.sin(lambda_le_rad)
    x_cg_wing = x_cg_wingstart + (location+dis)
    y_cg_wing = 0
    z_cg_wing = t_c*Cr/2
    return x_cg_wing, y_cg_wing, z_cg_wing

def get_cg_hwing(b_h,Cr_h,Ct_h,lambda_h_le_rad,x_le_h,d_f_outer):
    y_loc = 0.38
    dis = b_h/2*y_loc*np.sin(lambda_h_le_rad)
    location = 0.42*chordlength(y_loc, Cr_h, Ct_h,b_h)
    x_cg_hwing = location+dis+x_le_h
    y_cg_hwing = 0
    z_cg_hwing = 0.65*d_f_outer
    return x_cg_hwing, y_cg_hwing, z_cg_hwing


def get_cg_vwing(b_v,Cr_v,Ct_v):
    z_loc  = 0.38
    dis = b_v/2*z_loc*np.sin(lambda_v_le_rad)
    location = 0.42*chordlength(z_loc, Cr_v, Ct_v,b_v)
    x_cg_hwing = location+dis+x_le_v
    y_cg_hwing = 0
    z_cg_hwing = d_f_outer+0.38*(b_v/2)
    return x_cg_vwing, y_cg_vwing, z_cg_vwing
    
def get_cg_fuselage():
    x_cg_fuselage = 0.435 * l_fuselage
    y_cg_fuselage = 0
    z_cg_fuselage = 0.5*d_f_outer
    return x_cg_fuselage, y_cg_fuselage, z_cg_fuselage

def get_cg_nacelle():
    x_cg_nacelle = 0.4*l_nacel+x_cg_eng                                                 
    y_cg_nacelle = 0
    z_cg_nacelle = 0
    return x_cg_nacelle, y_cg_nacelle, z_cg_nacelle
    
def get_cg_ngear():
    x_cg_ngear = x_nlg
    y_cg_ngear = 0
    z_cg_ngear = 0.5*z_nlg
    return x_cg_ngear, y_cg_ngear, z_cg_ngear

def get_cg_mgear():
    x_cg_mgear = x_mlg
    y_cg_mgear = 0
    z_cg_mgear = 0.5*m_nlg
    return x_cg_mgear, y_cg_mgear, z_cg_mgear
    
  
