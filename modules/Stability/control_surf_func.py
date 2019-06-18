# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:57:56 2019

@author: Stijn
"""
import numpy as np

def chordlength(perc, Cr, Ct, b):
    decrease = (Cr-Ct)/(b/2)                                                    # [-] slope of decrease of the chord along the length
    return Cr-decrease*b/2*perc                                                 # [m] length of the chord at a certain y-location

# elevator

def get_S_elev(S_h):
    return 0.25 * S_h

def get_c_elev(Cr_h, Ct_h, b_h):
    c_elev_r = chordlength(0, Cr_h, Ct_h, b_h)*0.3
    c_elev_t = chordlength(1, Cr_h, Ct_h, b_h)*0.3
    
    return np.array([c_elev_r, c_elev_t])

def get_b_elev(S_e,c_e):
    return S_e / c_e

# rudder
    
def get_S_rud(S_v):
    return 0.325 * S_v

def get_c_rud(Cr_v, Ct_v, b_v):
    c_rud_r = chordlength(0, Cr_v ,Ct_v, b_v)*0.3
    c_rud_t = chordlength(1, Cr_v, Ct_v, b_v)*0.3
    
    return c_rud_r, c_rud_t

def get_b_rud(S_r,c_r):
    return S_r / c_r

# aileron
    
def get_S_ail(S):
    return 0.0325 * S

def get_c_ail(Cr, Ct, b):
    c_ail_r = chordlength(0.74, Cr ,Ct, b)*0.22
    c_ail_t = chordlength(0.95, Cr, Ct, b)*0.22
    
    return c_ail_r, c_ail_t

def get_b_ail(b):
    return 0.74 * b/2 , 0.95 * b/2

# spoiler
    
def get_c_splr(Cr, Ct, b):
    c_splr_r = chordlength(0.4, Cr ,Ct, b)*0.14
    c_splr_t = chordlength(0.65, Cr, Ct, b)*0.14
    
    return c_splr_r, c_splr_t

def get_b_splr(b):
    return 0.4 * b/2 , 0.65 * b/2
