# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:57:56 2019

@author: Stijn
"""

# elevator

def get_c_elev(c_h):
    return 0.3 * c_h

def get_S_elev(S_h):
    return 0.25 * S_h

def get_b_elev(S_e,c_e):
    return S_e / c_e

# rudder
    
def get_c_rud(c_v):
    return 0.3 * c_v

def get_S_rud(S_v):
    return 0.325 * S_v

def get_b_rud(S_r,c_r):
    return S_r / c_r

# aileron
    
def get_c_ail(c):
    return 0.22 * c

def get_S_ail(S):
    return 0.0325 * S

def get_b_ail(b):
    return 0.74 * b/2 , 0.95 * b/2

# spoiler
def get_c_splr(c):
    return 0.14 * c

def get_b_splr(b):
    return 0.4 * b/2 , 0.65 * b/2
