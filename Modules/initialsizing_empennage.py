# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:30:21 2019

@author: Stijn
"""
import numpy as np

def get_x_h(l_f):
    return [0.9* l_f[i] for i in range(3)]

def get_S_h(S, MAC, x_cg, V_h, x_h):
    return [V_h[i] * S[i] * MAC[i] / (x_h[i] - x_cg[i]) for i in range(3)]

def get_b_h(S_h, A_h):
    return [np.sqrt(S_h[i]*A_h[i]) for i in range(3)]

def get_Cr_h(S_h, taper_ratio_h, b_h):
    return [2*S_h[i]/((1+taper_ratio_h[i])*b_h[i]) for i in range(3)]

def get_Ct_h(Cr_h, taper_ratio_h):
    return [Cr_h[i] * taper_ratio_h[i] for i in range(3)]

def get_S_v(S, MAC, x_cg, V_v, x_v):
     return [V_v[i]*S[i]* MAC[i] / (x_v[i] - x_cg[i]) for i in range(3)]
 
def get_b_v(S_v, A_v):
    return [np.sqrt(S_v[i]*A_v[i]) for i in range(3)]

def get_Cr_v(S_v, taper_ratio_v, b_v):
    return [2*S_v[i]/((1+taper_ratio_v[i])*b_v[i]) for i in range(3)]

def get_Ct_v(Cr_v, taper_ratio_v):
    return [Cr_v[i] * taper_ratio_v[i] for i in range(3)]