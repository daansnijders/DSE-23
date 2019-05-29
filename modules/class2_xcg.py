# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:44:46 2019

@author: Niels
"""

def get_xcg_wing(b,Ct,Cr):
    dis = b/2*0.35*sin(lamda_le_rad)
    spar_front = 0.2* chordlength

def chordlength(x_loc):
    decrease = (Cr-Ct)/(b/2)
    return Cr-decrease*(b/2)*x_loc
    
def spar_front(x_loc):
    chord=chordlength(x_loc) 
    return chordlength
    